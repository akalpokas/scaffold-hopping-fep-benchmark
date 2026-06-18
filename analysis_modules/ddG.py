import numpy as np
import pandas as pd
from pathlib import Path
from .base import AnalysisModule
from rich.console import Console
from rich.table import Table
from BioSimSpace.FreeEnergy import Relative


class ddGAnalyzer(AnalysisModule):
    def __init__(self):
        super().__init__()
        self.console = Console()
        # Dictionary to store edge -> leg -> replicate data
        # Example shape: { "edge_1": { "free": { 1: 1.2, 2: 1.1 }, "bound": { ... } } }
        self.data = {}

    @property
    def name(self) -> str:
        return "ddG"

    def _mask_optimization_files(self, run_dir: Path, target_folder: str = "OptimizeLambdaProbabilities"):
        """
        Finds all .parquet files in the specified folder and appends .bak
        to their filenames to hide them from globbing API.
        """
        target_dir = Path(run_dir) / target_folder

        if not target_dir.is_dir():
            return

        parquet_files = list(target_dir.glob("*.parquet"))

        if not parquet_files:
            return

        for file_path in parquet_files:
            new_file_path = file_path.with_suffix(".bak")
            file_path.rename(new_file_path)

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
        Parses and records the dG value from a single replicate simulation.
        """
        # Initialize dictionary structure for this edge if it doesn't exist
        if edge_id not in self.data:
            self.data[edge_id] = {"free": {}, "bound": {}}

        self._mask_optimization_files(run_dir)
        
        try:
            pmf, overlap = Relative.analyse(str(run_dir))
            dg = Relative.difference(pmf)[0].value()
        except Exception as e:
            print(f"Error analyzing edge '{edge_id}', replicate {rep} in leg '{leg}': {e}")
            dg = np.nan

        # Store the single replicate value
        self.data[edge_id][leg][rep] = dg

    def aggregate_leg(self, edge_id: str, out_dir: Path, leg: str):
        """
        Calculates the mean and standard deviation for all parsed replicates of a given leg.
        """
        if edge_id not in self.data or leg not in self.data[edge_id]:
            return

        # Filter out NaN values from failed runs before calculating stats
        rep_values = [val for val in self.data[edge_id][leg].values() if not np.isnan(val)]

        if rep_values:
            mean_dg = np.mean(rep_values)
            # Use ddof=1 for sample standard deviation if >1 replicates, otherwise 0
            std_dg = np.std(rep_values, ddof=1) if len(rep_values) > 1 else 0.0
        else:
            mean_dg = np.nan
            std_dg = np.nan

        # Store aggregated metrics
        self.data[edge_id][f"{leg}_mean"] = mean_dg
        self.data[edge_id][f"{leg}_std"] = std_dg

    def compare_edge(self, edge_id: str, out_dir: Path):
        """
        Calculates the final ddG and outputs a formatted table summarizing the replicates,
        the leg aggregates, and the overall edge calculation.
        """
        if edge_id not in self.data:
            return

        edge_info = self.data[edge_id]

        free_mean = edge_info.get("free_mean", np.nan)
        free_std = edge_info.get("free_std", np.nan)
        bound_mean = edge_info.get("bound_mean", np.nan)
        bound_std = edge_info.get("bound_std", np.nan)

        # Calculate ddG (Bound - Free) and propagate the standard error
        if not np.isnan(free_mean) and not np.isnan(bound_mean):
            ddg = bound_mean - free_mean
            ddg_std = np.sqrt(free_std**2 + bound_std**2)
        else:
            ddg = np.nan
            ddg_std = np.nan

        # Save the result to the dictionary if downstream scripts/files need it
        self.data[edge_id]["ddG"] = ddg
        self.data[edge_id]["ddG_std"] = ddg_std

        # ==========================================
        # Render Table
        # ==========================================
        table = Table(
            title=f"Binding Free Energy Profile: [bold cyan]{edge_id}[/bold cyan]",
            show_lines=False,
            title_justify="center"
        )
        
        table.add_column("Leg", style="cyan", justify="left")
        table.add_column("Replicate / Summary", style="magenta", justify="center")
        table.add_column("dG / ddG (kcal/mol)", style="green", justify="right")

        # 1. Free Leg Sub-section
        first_free = True
        for rep, val in sorted(edge_info.get("free", {}).items()):
            leg_label = "[bold]Free[/bold]" if first_free else ""
            first_free = False
            val_str = f"{val:.2f}" if not np.isnan(val) else "NaN"
            table.add_row(leg_label, f"Rep {rep}", val_str)

        table.add_row("", "[bold]Mean ± Std[/bold]", f"[bold]{free_mean:.2f} ± {free_std:.2f}[/bold]")
        table.add_section() # Adds a horizontal rule

        # 2. Bound Leg Sub-section
        first_bound = True
        for rep, val in sorted(edge_info.get("bound", {}).items()):
            leg_label = "[bold]Bound[/bold]" if first_bound else ""
            first_bound = False
            val_str = f"{val:.2f}" if not np.isnan(val) else "NaN"
            table.add_row(leg_label, f"Rep {rep}", val_str)

        table.add_row("", "[bold]Mean ± Std[/bold]", f"[bold]{bound_mean:.2f} ± {bound_std:.2f}[/bold]")
        table.add_section() # Adds a horizontal rule

        # 3. Overall Edge ddG Calculation
        ddg_str = f"[bold yellow]{ddg:.2f} ± {ddg_std:.2f}[/bold yellow]" if not np.isnan(ddg) else "[bold red]NaN[/bold red]"
        table.add_row("[bold yellow]Overall Edge[/bold yellow]", "[bold yellow]Bound - Free[/bold yellow]", ddg_str)

        self.console.print(table)