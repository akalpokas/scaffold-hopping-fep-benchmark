import numpy as np
import pandas as pd
from pathlib import Path
from .base import AnalysisModule
from rich.console import Console
from rich.table import Table
from BioSimSpace.FreeEnergy import Relative

# Attempt to import plotting libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns

    PLOTTING_ENABLED = True
except ImportError:
    plt = None
    sns = None
    PLOTTING_ENABLED = False


class ddGAnalyzer(AnalysisModule):
    def __init__(self):
        super().__init__()
        self.console = Console()
        self.data = {}

        # Notify the user once upon initialization if plotting is disabled
        if not PLOTTING_ENABLED:
            self.console.print(
                "[yellow]Warning: matplotlib or seaborn not found. "
                "PMF plotting functionality will be disabled.[/yellow]"
            )

    @property
    def name(self) -> str:
        return "ddG"

    def _mask_optimization_files(
        self, run_dir: Path, target_folder: str = "OptimizeLambdaProbabilities"
    ):
        """
        Finds all .parquet files in the specified folder and appends .bak
        to their filenames to hide them from globbing API.
        """
        target_dir = Path(run_dir) / target_folder

        if not target_dir.is_dir():
            return

        for file_path in target_dir.glob("*.parquet"):
            file_path.rename(file_path.with_suffix(".bak"))

    def analyze_replicate(
        self,
        edge_id: str,
        run_dir: Path,
        out_dir: Path,
        leg: str,
        rep: int,
        T: float = 300.0,
        **kwargs,
    ):
        """
        Parses and records the dG value and PMF from a single replicate simulation.
        """
        # Dynamically initialize dictionary structures
        if edge_id not in self.data:
            self.data[edge_id] = {}
        if leg not in self.data[edge_id]:
            self.data[edge_id][leg] = {}

        self._mask_optimization_files(run_dir)

        try:
            pmf, overlap = Relative.analyse(str(run_dir))
            dg = Relative.difference(pmf)[0].value()
        except Exception as e:
            print(
                f"Error analyzing edge '{edge_id}', replicate {rep} in leg '{leg}': {e}"
            )
            dg = np.nan
            pmf = None

        # Store the single replicate value and PMF as a dictionary
        self.data[edge_id][leg][rep] = {"dg": dg, "pmf": pmf}

    def aggregate_leg(self, edge_id: str, out_dir: Path, leg: str):
        """
        Calculates the mean and standard deviation for all parsed replicates of a given leg,
        and plots the PMF values for all replicates if plotting is enabled.
        """
        if edge_id not in self.data or leg not in self.data[edge_id]:
            return

        # Filter out NaN values from failed runs before calculating stats
        valid_reps = {
            rep: data
            for rep, data in self.data[edge_id][leg].items()
            if not np.isnan(data["dg"])
        }

        rep_values = [data["dg"] for data in valid_reps.values()]

        if rep_values:
            mean_dg = np.mean(rep_values)
            std_dg = np.std(rep_values, ddof=1) if len(rep_values) > 1 else 0.0
        else:
            mean_dg = np.nan
            std_dg = np.nan

        # Store aggregated metrics
        self.data[edge_id][f"{leg}_mean"] = mean_dg
        self.data[edge_id][f"{leg}_std"] = std_dg

        # Plot PMFs if valid replicates exist AND plotting libraries are installed
        if valid_reps and PLOTTING_ENABLED:
            fig, ax = plt.subplots(figsize=(10, 6))

            for rep, data in valid_reps.items():
                pmf = data["pmf"]
                if pmf is not None:
                    lam_vals = [pmf[i][0] for i in range(len(pmf))]
                    pmf_vals = [pmf[i][1].value() for i in range(len(pmf))]

                    sns.lineplot(x=lam_vals, y=pmf_vals, ax=ax, label=f"Replica {rep}")

                    # Save this replicate's PMF data
                    pmf_df = pd.DataFrame({"lambda": lam_vals, "pmf": pmf_vals})
                    csv_path = out_dir / f"{edge_id}_{leg}_rep_{rep}_pmf.csv"
                    pmf_df.to_csv(csv_path, index=False)

            ax.set_xlabel("Lambda")
            ax.set_ylabel("Free Energy (kcal/mol)")
            ax.set_title(f"PMF Profile: {edge_id} ({leg.capitalize()})")

            # Save the plot
            plot_path = out_dir / f"{edge_id}_{leg}_pmf.png"
            fig.savefig(plot_path, dpi=300, bbox_inches="tight")
            plt.close(fig)

    def compare_edge(self, edge_id: str, out_dir: Path):
        """
        Determines the system type based on simulated legs, calculates the final ddG
        (or dG for nodes), and outputs a formatted table summarizing the results.
        """
        if edge_id not in self.data:
            return

        edge_info = self.data[edge_id]

        # Identify populated legs dynamically (ignoring mean/std keys)
        legs = [k for k in edge_info.keys() if isinstance(edge_info[k], dict)]

        is_relative = False
        calc_label = ""
        calc_name = "dG"

        # Determine calculation logic based on expected leg pairs
        if "free" in legs and "bound" in legs:
            is_relative = True
            leg_ref, leg_target = "free", "bound"
            calc_label = "Bound - Free"
            calc_name = "ddG"
            leg_order = ["free", "bound"]
            title_type = "Binding"
        elif "vacuum" in legs and "solvent" in legs:
            is_relative = True
            leg_ref, leg_target = "vacuum", "solvent"
            calc_label = "Solvent - Vacuum"
            calc_name = "ddG"
            leg_order = ["vacuum", "solvent"]
            title_type = "Hydration"
        else:
            # Single-leg executions (e.g., AHFE or incomplete runs)
            leg_order = sorted(legs)
            title_type = "Node / Single Leg"

        # Execute subtraction and error propagation if relative
        if is_relative:
            ref_mean = edge_info.get(f"{leg_ref}_mean", np.nan)
            ref_std = edge_info.get(f"{leg_ref}_std", np.nan)
            target_mean = edge_info.get(f"{leg_target}_mean", np.nan)
            target_std = edge_info.get(f"{leg_target}_std", np.nan)

            if not np.isnan(ref_mean) and not np.isnan(target_mean):
                final_val = target_mean - ref_mean
                final_std = np.sqrt(ref_std**2 + target_std**2)
            else:
                final_val, final_std = np.nan, np.nan

            self.data[edge_id][calc_name] = final_val
            self.data[edge_id][f"{calc_name}_std"] = final_std

        # ==========================================
        # Render Table
        # ==========================================
        table = Table(
            title=f"{title_type} Free Energy Profile: [bold cyan]{edge_id}[/bold cyan]",
            show_lines=False,
            title_justify="center",
        )

        table.add_column("Leg", style="cyan", justify="left")
        table.add_column("Replicate / Summary", style="magenta", justify="center")
        table.add_column("dG / ddG (kcal/mol)", style="green", justify="right")

        # Dynamically build table sections based on `leg_order`
        for leg in leg_order:
            first = True
            leg_mean = edge_info.get(f"{leg}_mean", np.nan)
            leg_std = edge_info.get(f"{leg}_std", np.nan)

            for rep, data in sorted(edge_info.get(leg, {}).items()):
                val = data["dg"]
                leg_label = f"[bold]{leg.capitalize()}[/bold]" if first else ""
                first = False
                val_str = f"{val:.2f}" if not np.isnan(val) else "NaN"
                table.add_row(leg_label, f"Rep {rep}", val_str)

            mean_str = (
                f"[bold]{leg_mean:.2f} ± {leg_std:.2f}[/bold]"
                if not np.isnan(leg_mean)
                else "[bold red]NaN[/bold red]"
            )
            table.add_row("", "[bold]Mean ± Std[/bold]", mean_str)
            table.add_section()

        # Footer calculation formatting
        if is_relative:
            val_str = (
                f"[bold yellow]{final_val:.2f} ± {final_std:.2f}[/bold yellow]"
                if not np.isnan(final_val)
                else "[bold red]NaN[/bold red]"
            )
            table.add_row(
                "[bold yellow]Overall Edge[/bold yellow]",
                f"[bold yellow]{calc_label}[/bold yellow]",
                val_str,
            )
        elif len(legs) == 1:
            # For AHFE nodes or single-leg runs, the node result is just the isolated leg's average
            single_leg = legs[0]
            single_val = edge_info.get(f"{single_leg}_mean", np.nan)
            single_std = edge_info.get(f"{single_leg}_std", np.nan)
            val_str = (
                f"[bold yellow]{single_val:.2f} ± {single_std:.2f}[/bold yellow]"
                if not np.isnan(single_val)
                else "[bold red]NaN[/bold red]"
            )
            table.add_row(
                "[bold yellow]Overall Node[/bold yellow]",
                f"[bold yellow]{single_leg.capitalize()} Only[/bold yellow]",
                val_str,
            )

        self.console.print(table)
