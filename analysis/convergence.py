import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import glob
import logging

logger = logging.getLogger(__name__)

# HELPER FUNCTIONS


def _analyse_somd2(work_dir, T=300.0, use_checkpoint_files=False, **kwargs):
    """
    Analyzes the results of SOMD2 simulations.
    """
    import pathlib as _pathlib

    files = glob.glob(f"{work_dir}/energy_traj_*.parquet")
    glob_path = _pathlib.Path(work_dir)

    analysed_df_list = []

    if use_checkpoint_files:
        files = sorted(glob_path.glob("**/*.s3"))
        logger.debug(f"Found {len(files)} checkpoint files in {work_dir}")
        for f in files:
            logger.debug(f"Loading checkpoint file: {f}")
            import sire as sr

            system = sr.stream.load(str(f))
            analysed_df_list.append(system.energy_trajectory().to_alchemlyb())
    else:
        files = sorted(glob_path.glob("**/energy_traj_*.parquet"))
        for f in files:
            path = Path(f)
            import BioSimSpace as BSS

            analysed_df = BSS.FreeEnergy.Relative._somd2_extract(path, T=T, **kwargs)
            analysed_df_list.append(analysed_df)

    # check if all dataframes contain the same number of entries
    lengths = [len(df) for df in analysed_df_list]
    if len(set(lengths)) != 1:
        logger.warning(
            f"Not all dataframes contain the same number of entries: {lengths}"
        )

    return analysed_df_list


def _plot_alchemlyb_convergence(
    analysed_df_list, estimator="MBAR", units="kcal/mol", number_of_points=10
):
    """
    Computes the forward-backward convergence for a list of alchemlyb dataframes.
    """
    from alchemlyb.convergence import forward_backward_convergence
    from alchemlyb.postprocessors.units import to_kcalmol
    from alchemlyb.visualisation import plot_convergence

    convergence_df = forward_backward_convergence(
        analysed_df_list, estimator=estimator, num=number_of_points
    )

    if units == "kcal/mol":
        convergence_df = to_kcalmol(convergence_df)

    ax = plot_convergence(convergence_df)

    if units == "kcal/mol":
        ax.set_ylabel(r"$\Delta G$ ({})".format("kcal/mol"))

    return convergence_df


# PLUGIN CLASS


class ConvergenceAnalyzer:
    """Plugin to analyze and plot alchemical free energy convergence."""

    def __init__(self, base_strategies=None, override_mismatch=False):
        self._name = "convergence"
        self.override_mismatch = override_mismatch

        # 1. DEFINE STRATEGIES
        self.base_strategies = base_strategies or {
            "Full": {"t_min": 0, "t_max": None},
            "Manual Discard": {"t_min": 4999, "t_max": None},
            "Rapid": {"t_min": 1001, "t_max": 5001},
        }

        # Dictionary to track expected sampling times across replicates
        self.expected_max_ps = {}

    @property
    def name(self) -> str:
        return self._name

    def analyze_replicate(
        self, edge_id: str, run_dir: Path, out_dir: Path, leg: str, rep: int
    ):
        """1. Analyzes a single trajectory and caches the convergence dataframe."""
        print(f"[{self.name}] Processing {edge_id} | {leg} | replicate {rep}...")

        # Load raw data using the helper function
        try:
            analysed_df_list = _analyse_somd2(str(run_dir))
        except Exception as e:
            print(f"    [!] Failed to load or analyze {run_dir}: {e}")
            return

        print(
            f"    Analysed dataframes for Replicate {rep}: {[df.shape for df in analysed_df_list]}"
        )

        # detect global maximum time
        global_max_ps = int(analysed_df_list[0][0].index.get_level_values("time").max())

        # Validation Check against Replicate 1 (or whichever runs first)
        leg_key = f"{edge_id}_{leg}"
        if leg_key not in self.expected_max_ps:
            self.expected_max_ps[leg_key] = global_max_ps
        elif (
            global_max_ps != self.expected_max_ps[leg_key]
            and not self.override_mismatch
        ):
            raise ValueError(
                f"Sampling time mismatch for {leg_key}! Expected max time of {self.expected_max_ps[leg_key]/1000} ns, "
                f"but Replicate {rep} has a max time of {global_max_ps/1000} ns. "
                "All replicates must be strictly equal for accurate convergence analysis."
            )

        global_max_ns = round(global_max_ps / 1000)
        convergence_dfs = []

        for base_name, bounds in self.base_strategies.items():
            trimmed_analysed_df_list = []

            # Apply the time mask for this strategy
            for df in analysed_df_list:
                time_vals = df.index.get_level_values("time")
                mask = time_vals >= bounds["t_min"]
                if bounds["t_max"] is not None:
                    mask = mask & (time_vals <= bounds["t_max"])
                trimmed_analysed_df_list.append(df[mask])

            # Calculate accurate times for the trimmed block
            start_time = int(
                trimmed_analysed_df_list[0][0].index.get_level_values("time").min()
            )
            start_ns = start_time / 1000

            # Create title label
            end_ns_label = (
                round(bounds["t_max"] / 1000)
                if bounds["t_max"] is not None
                else global_max_ns
            )
            strategy_name = f"{base_name} ({round(start_ns)} - {end_ns_label} ns)"

            # Run alchemlyb calculation using the helper function
            convergence_df = _plot_alchemlyb_convergence(
                trimmed_analysed_df_list, number_of_points=10
            )

            # Normalize data_fraction column
            convergence_df["data_fraction"] = (
                convergence_df["data_fraction"] / convergence_df["data_fraction"].max()
            )

            plt.close(
                "all"
            )  # Silently close the unwanted figure(s) generated by alchemlyb

            convergence_df["replicate"] = rep
            convergence_df["strategy"] = strategy_name

            convergence_dfs.append(convergence_df)

        # Cache this replicate's data so aggregate_leg can pick it up
        combined_df = pd.concat(convergence_dfs, ignore_index=True)
        csv_path = out_dir / f"{edge_id}_{leg}_rep_{rep}_convergence.csv"
        combined_df.to_csv(csv_path, index=False)

    def aggregate_leg(self, edge_id: str, out_dir: Path, leg: str):
        """2. Aggregates the replicates for a single leg and plots the results."""
        print(f"[{self.name}] Aggregating {edge_id} | {leg}...")

        # Read the cached CSVs for this leg
        csv_files = list(out_dir.glob(f"{edge_id}_{leg}_rep_*_convergence.csv"))
        if not csv_files:
            print(
                f"    [!] No convergence data found to aggregate for {edge_id} {leg}."
            )
            return

        combined_dfs = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)

        # 3. SUMMARY STATISTICS
        final_rows = (
            combined_dfs.sort_values("data_fraction")
            .groupby(["strategy", "replicate"])
            .tail(1)
        )
        summary_df = (
            final_rows.groupby("strategy")
            .agg(
                Mean_Final_Forward=("Forward", "mean"),
                Std_Final_Forward=("Forward", "std"),
                Mean_Final_Backward=("Backward", "mean"),
                Std_Final_Backward=("Backward", "std"),
            )
            .reset_index()
        )

        print(
            f"\nFinal Free Energy Estimates Summary ({edge_id} - {leg.capitalize()}):"
        )
        print(summary_df.to_string(index=False))

        # 4. PLOTTING
        sns.set_theme(style="ticks", context="paper", font_scale=1.3)

        strategy_labels = combined_dfs["strategy"].unique()
        n_strats = len(strategy_labels)

        # Dynamically size figure based on number of strategies
        fig, axes = plt.subplots(1, n_strats, figsize=(4 * n_strats, 7.41), sharey=True)
        if n_strats == 1:
            axes = [axes]

        for ax, strategy in zip(axes, strategy_labels):
            subset_df = combined_dfs[combined_dfs["strategy"] == strategy]

            # Plot individual replicate traces in the background
            for rep in subset_df["replicate"].unique():
                rep_df = subset_df[subset_df["replicate"] == rep]

                # Forward individual replicate
                ax.plot(
                    rep_df["data_fraction"],
                    rep_df["Forward"],
                    color="#007FFF",
                    alpha=0.25,
                    linewidth=1,
                    zorder=1,
                )

                # Backward individual replicate
                ax.plot(
                    rep_df["data_fraction"],
                    rep_df["Backward"],
                    color="#FF4500",
                    alpha=0.25,
                    linewidth=1,
                    zorder=1,
                )

            # Plot aggregated mean lines using the numerical data_fraction on the x-axis
            sns.lineplot(
                x="data_fraction",
                y="Forward",
                data=subset_df,
                ax=ax,
                linewidth=2.5,
                color="#007FFF",
                label="Forward Mean",
                errorbar="sd",
                marker="o",
                markersize=6,
                zorder=2,
            )

            sns.lineplot(
                x="data_fraction",
                y="Backward",
                data=subset_df,
                ax=ax,
                linewidth=2.5,
                color="#FF4500",
                label="Backward Mean",
                errorbar="sd",
                marker="s",
                markersize=6,
                zorder=2,
            )

            # Set standard ticks from 0.0 to 1.0
            ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])

            # Formatting
            ax.set_title(strategy, weight="bold", pad=15)
            ax.set_xlabel("Fraction of Data Analyzed", weight="bold", labelpad=10)

            if ax == axes[0]:
                ax.set_ylabel("Convergence ΔG (kcal/mol)", weight="bold")
            else:
                ax.set_ylabel("")

            sns.despine(ax=ax)
            ax.legend(frameon=False, loc="best")

        plt.suptitle(
            f"ΔG Convergence: {edge_id.capitalize()} ({leg.capitalize()})",
            weight="bold",
            y=1.05,
            fontsize=18,
        )
        plt.tight_layout()

        # Save to the specific analysis folder mapped by the orchestrator
        plot_filename = f"{edge_id}_{leg}_convergence"
        plt.savefig(out_dir / f"{plot_filename}.pdf", dpi=300, bbox_inches="tight")
        plt.savefig(out_dir / f"{plot_filename}.png", dpi=300, bbox_inches="tight")

        # Close plot rather than showing it to prevent script blocking
        plt.close("all")

    def compare_edge(self, edge_id: str, out_dir: Path):
        """3. Compares Bound vs Free for the entire edge. (Optional for convergence mapping)"""
        pass
