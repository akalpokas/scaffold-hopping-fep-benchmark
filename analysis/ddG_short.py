"""
ddG_short
=========

Time-truncated free energy analysis module.

Motivation
----------
`BioSimSpace.FreeEnergy.Relative.analyse()` consumes every sample present in the
parquet energy trajectories and offers no way to restrict the analysis to a time
window. To answer "what would this edge have given me if I had only run 5 ns?"
we therefore have to bypass `Relative.analyse` and rebuild the MBAR estimate
ourselves from the same underlying `u_nk` dataframes that the convergence module
already loads via `Relative._somd2_extract`.

Everything downstream of the per-replicate free energy (mean/std over
replicates, bound-free subtraction, error propagation, rich table rendering,
PMF plotting) is inherited unchanged from `ddGAnalyzer`, so a 5 ns number and a
10 ns number are aggregated and reported through exactly the same code path.

Usage
-----
    from .ddg_short import ddGShortAnalyzer

    modules = [
        ddGAnalyzer(),                              # full 10 ns
        ddGShortAnalyzer(t_max_ps=5000),            # first 5 ns
    ]

Point the orchestrator at a *separate* out_dir per module, or rely on the
``file_tag`` suffix applied to this module's outputs.

Validation
----------
Before trusting the comparison, run this module once with ``t_max_ps=None``
(no truncation) and confirm it reproduces your `ddG` module to within ~0.05
kcal/mol. Any residual difference is the gap between this MBAR call and
whatever preprocessing `Relative.analyse` applies internally, and it will be
common to both the 5 ns and 10 ns numbers rather than an artefact of
truncation. `check_against_full_analysis()` at the bottom of this file does
that check for a single directory.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .base import AnalysisModule  # noqa: F401  (kept for registry consistency)
from .ddG import ddGAnalyzer

logger = logging.getLogger(__name__)

# Gas constant in kcal/(mol K); converts reduced potentials (kT) -> kcal/mol
_R_KCAL = 0.0019872041

_DEFAULT_EXCLUDE_DIRS = ("OptimizeLambdaProbabilities",)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _EnergyShim:
    """
    Minimal stand-in for a BioSimSpace energy type.

    The parent class calls ``pmf[i][1].value()`` when plotting, so the PMF we
    hand back must expose ``.value()``. We prefer real BSS units when they are
    importable and fall back to this shim otherwise.
    """

    __slots__ = ("_value", "_unit")

    def __init__(self, value, unit="kcal/mol"):
        self._value = float(value)
        self._unit = unit

    def value(self):
        return self._value

    def unit(self):
        return self._unit

    def __float__(self):
        return self._value

    def __repr__(self):
        return f"{self._value:.4f} {self._unit}"


def _as_energy(value):
    """Wrap a float in a BSS energy object, falling back to the local shim."""
    try:
        from BioSimSpace.Units.Energy import kcal_per_mol

        return float(value) * kcal_per_mol
    except Exception:
        return _EnergyShim(value)


def _extract_u_nk(run_dir, T=300.0, exclude_dirs=_DEFAULT_EXCLUDE_DIRS, **kwargs):
    """
    Load one `u_nk` dataframe per lambda window from a SOMD2 run directory.

    Unlike the convergence module's `_analyse_somd2`, this does **not** delete
    the OptimizeLambdaProbabilities folder -- it simply refuses to glob into it.
    Deleting is destructive and would prevent re-running the full-length ddG
    analysis afterwards. Files already renamed to `.bak` by
    `ddGAnalyzer._mask_optimization_files` are likewise ignored, since we only
    match `*.parquet`.
    """
    import BioSimSpace as BSS

    run_dir = Path(run_dir)
    files = sorted(run_dir.glob("**/energy_traj_*.parquet"))
    files = [f for f in files if not any(d in f.parts for d in exclude_dirs)]

    if not files:
        raise FileNotFoundError(
            f"No 'energy_traj_*.parquet' files found under {run_dir}."
        )

    df_list = [BSS.FreeEnergy.Relative._somd2_extract(f, T=T, **kwargs) for f in files]

    lengths = {len(df) for df in df_list}
    if len(lengths) != 1:
        logger.warning(
            "Lambda windows in %s have unequal sample counts: %s", run_dir, lengths
        )

    return df_list


def _truncate(df, t_min_ps, t_max_ps):
    """Mask a u_nk dataframe on its 'time' index level (units: ps)."""
    time_vals = df.index.get_level_values("time")
    mask = time_vals >= t_min_ps
    if t_max_ps is not None:
        mask = mask & (time_vals <= t_max_ps)
    return df[mask]


def _lambda_axis(index):
    """
    Coerce an MBAR state index into floats for plotting.

    Handles both scalar lambda schedules and multi-dimensional (tuple) lambdas,
    falling back to evenly spaced indices if the labels are not numeric.
    """
    values = []
    for entry in index:
        if isinstance(entry, tuple):
            entry = entry[0]
        try:
            values.append(float(entry))
        except (TypeError, ValueError):
            values = None
            break

    if values is None or len(set(values)) != len(values):
        n = len(index)
        return list(np.linspace(0.0, 1.0, n)) if n > 1 else [0.0]

    return values


def _mbar_pmf(df_list, T=300.0, decorrelate=False):
    """
    Run MBAR over a list of (already truncated) u_nk dataframes.

    Returns
    -------
    pmf : list of (lambda, energy) tuples, referenced to lambda = 0
    dg : float, total free energy difference in kcal/mol
    dg_err : float, MBAR statistical uncertainty in kcal/mol
    n_samples : int, total number of samples entering the estimate
    """
    from alchemlyb import concat
    from alchemlyb.estimators import MBAR

    if decorrelate:
        from alchemlyb.preprocessing import decorrelate_u_nk

        processed = []
        for i, df in enumerate(df_list):
            try:
                processed.append(decorrelate_u_nk(df, remove_burnin=False))
            except Exception as exc:
                logger.warning(
                    "Decorrelation failed for window %d (using raw data): %s", i, exc
                )
                processed.append(df)
        df_list = processed

    u_nk = concat(df_list)
    n_samples = len(u_nk)

    mbar = MBAR().fit(u_nk)

    kT = _R_KCAL * T
    deltas = mbar.delta_f_.iloc[0] * kT
    errors = mbar.d_delta_f_.iloc[0] * kT

    lam_vals = _lambda_axis(mbar.delta_f_.index)
    pmf = [(lam, _as_energy(val)) for lam, val in zip(lam_vals, deltas.values)]

    return pmf, float(deltas.values[-1]), float(errors.values[-1]), n_samples


# ---------------------------------------------------------------------------
# Plugin class
# ---------------------------------------------------------------------------


class ddGShortAnalyzer(ddGAnalyzer):
    """
    Free energy analysis restricted to a fixed simulation-time window.

    Parameters
    ----------
    t_max_ps : float or None
        Upper bound of the analysis window in ps. Default 5000 (5 ns).
        ``None`` disables the upper bound, which is the mode to use when
        validating this module against the full `ddG` analyzer.
    t_min_ps : float
        Lower bound in ps. Default 0. Raise this to discard equilibration and
        keep the *amount* of sampling fixed (e.g. 1000 -> 6000 gives 5 ns of
        production data after a 1 ns burn-in).
    T : float
        Temperature in K, used both for the reduced-potential extraction and
        the kT -> kcal/mol conversion.
    decorrelate : bool
        Apply `decorrelate_u_nk` (without burn-in removal) before estimating.
        Off by default so that the short and full analyses differ only in the
        length of data consumed.
    strict : bool
        If True, raise when a replicate contains less data than the requested
        window. If False (default), warn and analyse what is available -- but
        the shortfall is recorded so it shows up in the summary CSV.
    file_tag : str
        Suffix appended to this module's output files to avoid collisions with
        the full-length ddG module if both write to the same directory.
    """

    def __init__(
        self,
        t_max_ps=5000.0,
        t_min_ps=0.0,
        T=300.0,
        decorrelate=False,
        strict=False,
        file_tag=None,
    ):
        super().__init__()

        if t_max_ps is not None and t_max_ps <= t_min_ps:
            raise ValueError(
                f"t_max_ps ({t_max_ps}) must be greater than t_min_ps ({t_min_ps})."
            )

        self.t_min_ps = float(t_min_ps)
        self.t_max_ps = None if t_max_ps is None else float(t_max_ps)
        self.T = float(T)
        self.decorrelate = bool(decorrelate)
        self.strict = bool(strict)

        if file_tag is not None:
            self.file_tag = file_tag
        elif self.t_max_ps is None:
            self.file_tag = "full"
        else:
            self.file_tag = f"{self._window_ns():g}ns"

        # rep-level provenance: how much data each replicate actually supplied
        self.meta = {}

    # -- identity ----------------------------------------------------------

    @property
    def name(self) -> str:
        return "ddG_short"

    def _window_ns(self):
        if self.t_max_ps is None:
            return np.nan
        return (self.t_max_ps - self.t_min_ps) / 1000.0

    def _window_label(self):
        if self.t_max_ps is None:
            return f"{self.t_min_ps / 1000:g} ns - end"
        return f"{self.t_min_ps / 1000:g} - {self.t_max_ps / 1000:g} ns"

    # -- per replicate -----------------------------------------------------

    def analyze_replicate(
        self,
        edge_id: str,
        run_dir: Path,
        out_dir: Path,
        leg: str,
        rep: int,
        T: float = None,
        **kwargs,
    ):
        """
        Estimate dG for one replicate using only the configured time window.

        Signature matches the parent so the orchestrator can call either module
        interchangeably. Note that ``_mask_optimization_files`` is deliberately
        *not* called here -- we never invoke `Relative.analyse`, so there is
        nothing to hide from it.
        """
        temperature = self.T if T is None else float(T)

        self.data.setdefault(edge_id, {}).setdefault(leg, {})
        self.meta.setdefault(edge_id, {}).setdefault(leg, {})

        dg = np.nan
        dg_err = np.nan
        pmf = None
        record = {
            "requested_window_ns": self._window_ns(),
            "available_ns": np.nan,
            "used_ns": np.nan,
            "n_samples": 0,
            "truncated": False,
            "short": False,
            "status": "ok",
        }

        try:
            df_list = _extract_u_nk(run_dir, T=temperature, **kwargs)

            available_ps = float(
                min(df.index.get_level_values("time").max() for df in df_list)
            )
            record["available_ns"] = available_ps / 1000.0

            if self.t_max_ps is not None:
                if available_ps + 1e-6 < self.t_max_ps:
                    record["short"] = True
                    msg = (
                        f"{edge_id} | {leg} | rep {rep}: only "
                        f"{available_ps / 1000:.2f} ns available, "
                        f"but window extends to {self.t_max_ps / 1000:.2f} ns."
                    )
                    if self.strict:
                        raise ValueError(msg)
                    logger.warning(msg)
                    self.console.print(f"[yellow]Warning: {msg}[/yellow]")
                else:
                    record["truncated"] = True

            df_list = [_truncate(df, self.t_min_ps, self.t_max_ps) for df in df_list]

            if any(len(df) == 0 for df in df_list):
                raise ValueError(
                    "At least one lambda window has no samples inside "
                    f"[{self.t_min_ps}, {self.t_max_ps}] ps."
                )

            used_ps = float(
                min(df.index.get_level_values("time").max() for df in df_list)
            ) - float(max(df.index.get_level_values("time").min() for df in df_list))
            record["used_ns"] = used_ps / 1000.0

            pmf, dg, dg_err, n_samples = _mbar_pmf(
                df_list, T=temperature, decorrelate=self.decorrelate
            )
            record["n_samples"] = n_samples

        except Exception as exc:
            record["status"] = f"failed: {exc}"
            self.console.print(
                f"[red]Error analyzing edge '{edge_id}', replicate {rep} "
                f"in leg '{leg}' ({self._window_label()}): {exc}[/red]"
            )

        # Parent's aggregate_leg reads "dg" and "pmf"; extra keys are harmless.
        self.data[edge_id][leg][rep] = {
            "dg": dg,
            "pmf": pmf,
            "dg_err": dg_err,
            "n_samples": record["n_samples"],
        }
        self.meta[edge_id][leg][rep] = record

    # -- per leg -----------------------------------------------------------

    def aggregate_leg(self, edge_id: str, out_dir: Path, leg: str):
        """
        Reuse the parent aggregation/plotting, then tag the outputs and write a
        provenance CSV recording how much data each replicate contributed.
        """
        super().aggregate_leg(edge_id, out_dir, leg)

        out_dir = Path(out_dir)
        prefix = f"{edge_id}_{leg}"

        # Tag parent-generated files so both modules can share an out_dir.
        for path in list(out_dir.glob(f"{prefix}_rep_*_pmf.csv")) + list(
            out_dir.glob(f"{prefix}_pmf.png")
        ):
            if f"_{self.file_tag}" in path.stem:
                continue
            tagged = path.with_name(f"{path.stem}_{self.file_tag}{path.suffix}")
            path.replace(tagged)

        records = self.meta.get(edge_id, {}).get(leg, {})
        if records:
            meta_df = pd.DataFrame(
                [{"replicate": rep, **info} for rep, info in sorted(records.items())]
            )
            meta_df.insert(0, "leg", leg)
            meta_df.insert(0, "edge", edge_id)
            meta_df.to_csv(
                out_dir / f"{prefix}_{self.file_tag}_sampling.csv",
                index=False,
                float_format="%.3f",
            )

            if meta_df["short"].any():
                self.console.print(
                    f"[yellow]{edge_id} | {leg}: {int(meta_df['short'].sum())} "
                    f"replicate(s) had less data than the requested "
                    f"{self._window_ns():g} ns window.[/yellow]"
                )

    # -- per edge ----------------------------------------------------------

    def compare_edge(self, edge_id: str, out_dir: Path):
        """
        Render the parent's summary table, then append a machine-readable row so
        that windows can be compared across runs (see `compare_windows`).
        """
        self.console.print(
            f"[dim]Analysis window: {self._window_label()} "
            f"({'truncated' if self.t_max_ps is not None else 'full length'})[/dim]"
        )

        super().compare_edge(edge_id, out_dir)

        if edge_id not in self.data:
            return

        edge_info = self.data[edge_id]
        legs = [k for k in edge_info.keys() if isinstance(edge_info[k], dict)]

        if "ddG" in edge_info:
            value = edge_info["ddG"]
            error = edge_info.get("ddG_std", np.nan)
            quantity = "ddG"
        elif len(legs) == 1:
            value = edge_info.get(f"{legs[0]}_mean", np.nan)
            error = edge_info.get(f"{legs[0]}_std", np.nan)
            quantity = "dG"
        else:
            value, error, quantity = np.nan, np.nan, "dG"

        row = {
            "edge": edge_id,
            "quantity": quantity,
            "window_ns": self._window_ns(),
            "t_min_ps": self.t_min_ps,
            "t_max_ps": self.t_max_ps if self.t_max_ps is not None else np.nan,
            "value_kcal_mol": value,
            "std_kcal_mol": error,
            "n_replicates": len(edge_info.get(legs[0], {})) if legs else 0,
        }
        for leg in sorted(legs):
            row[f"{leg}_mean"] = edge_info.get(f"{leg}_mean", np.nan)
            row[f"{leg}_std"] = edge_info.get(f"{leg}_std", np.nan)

        summary_path = Path(out_dir) / f"ddG_summary_{self.file_tag}.csv"
        frame = pd.DataFrame([row])
        if summary_path.exists():
            existing = pd.read_csv(summary_path)
            existing = existing[existing["edge"] != edge_id]
            frame = pd.concat([existing, frame], ignore_index=True)
        frame.to_csv(summary_path, index=False, float_format="%.4f")


# ---------------------------------------------------------------------------
# Cross-window comparison
# ---------------------------------------------------------------------------


def compare_windows(
    short_csv, full_csv, out_csv=None, label_short=None, label_full=None
):
    """
    Merge two `ddG_summary_*.csv` files and report per-edge differences.

    Returns a dataframe with the two estimates side by side, their difference,
    and whether that difference is covered by the combined replicate spread
    (a crude but useful "did 5 ns already get me the answer?" flag).
    """
    short_df = pd.read_csv(short_csv)
    full_df = pd.read_csv(full_csv)

    label_short = label_short or f"{short_df['window_ns'].iloc[0]:g}ns"
    label_full = label_full or (
        "full"
        if np.isnan(full_df["window_ns"].iloc[0])
        else f"{full_df['window_ns'].iloc[0]:g}ns"
    )

    keep = ["edge", "quantity", "value_kcal_mol", "std_kcal_mol"]
    merged = short_df[keep].merge(
        full_df[keep],
        on=["edge", "quantity"],
        suffixes=(f"_{label_short}", f"_{label_full}"),
    )

    v_s = merged[f"value_kcal_mol_{label_short}"]
    v_f = merged[f"value_kcal_mol_{label_full}"]
    e_s = merged[f"std_kcal_mol_{label_short}"]
    e_f = merged[f"std_kcal_mol_{label_full}"]

    merged["difference"] = v_f - v_s
    merged["combined_std"] = np.sqrt(e_s**2 + e_f**2)
    merged["within_error"] = merged["difference"].abs() <= merged["combined_std"]

    finite = merged["difference"].dropna()
    if len(finite):
        print(
            f"MUE  ({label_full} vs {label_short}): {finite.abs().mean():.3f} kcal/mol"
        )
        print(
            f"RMSE ({label_full} vs {label_short}): {np.sqrt((finite**2).mean()):.3f} kcal/mol"
        )
        print(f"Max drift: {finite.abs().max():.3f} kcal/mol")
        print(
            f"Edges agreeing within combined error: "
            f"{int(merged['within_error'].sum())}/{len(merged)}"
        )

    if out_csv is not None:
        merged.to_csv(out_csv, index=False, float_format="%.4f")

    return merged


def check_against_full_analysis(run_dir, T=300.0, decorrelate=False):
    """
    Sanity check: compare this module's untruncated MBAR estimate against
    `Relative.analyse` on the same directory. Run once before trusting any
    5 ns vs 10 ns comparison.
    """
    from BioSimSpace.FreeEnergy import Relative

    df_list = _extract_u_nk(run_dir, T=T)
    _, dg_local, err_local, n = _mbar_pmf(df_list, T=T, decorrelate=decorrelate)

    pmf, _ = Relative.analyse(str(run_dir))
    dg_bss = Relative.difference(pmf)[0].value()

    print(
        f"ddG_short (untruncated) : {dg_local:8.4f} ± {err_local:.4f} kcal/mol "
        f"({n} samples)"
    )
    print(f"Relative.analyse        : {dg_bss:8.4f} kcal/mol")
    print(f"Difference              : {dg_local - dg_bss:8.4f} kcal/mol")

    return dg_local, dg_bss
