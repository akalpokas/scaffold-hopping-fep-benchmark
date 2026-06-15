import pandas as pd
from pathlib import Path
from .base import AnalysisModule
from rich.console import Console
from rich.table import Table
from BioSimSpace.FreeEnergy import Relative


class ddGAnalyzer(AnalysisModule):

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

        # 1. Check if the directory actually exists first
        if not target_dir.is_dir():
            print(f"Notice: Directory '{target_dir}' not found. Skipping masking.")
            return

        # 2. Find all .parquet files
        parquet_files = list(target_dir.glob("*.parquet"))

        if not parquet_files:
            print(f"Notice: No .parquet files found in '{target_dir}'.")
            return

        # 3. Rename them by appending .bak
        print(f"Found {len(parquet_files)} .parquet file(s). Masking them now...")
        for file_path in parquet_files:
            # Create the new path by appending .bak to the current name
            new_file_path = file_path.with_suffix(".bak")

            # Perform the rename
            file_path.rename(new_file_path)

    def analyze_replicate(
        self,
        run_dir: Path,
        out_dir: Path,
        leg: str,
        rep: int,
        T: float = 300.0,
        **kwargs,
    ):
        """
        Reports the dG values from the individual edge simulations.

        Parameters:
        - run_dir (Path): The directory containing the parquet files.
        - T (float): The temperature at which the simulations were performed. Default is 300.0.

        Returns:
        - lengths (list): A list of lengths of dataframes.
        """

        self.console = Console()
        self.table = Table(title="dG Values", show_lines=True)
        self.table.add_column("Path", style="cyan", no_wrap=True)
        self.table.add_column("dG Value", style="magenta")

        self._mask_optimization_files(run_dir)
        try:
            pmf, overlap = Relative.analyse(str(run_dir))
            dg = Relative.difference(pmf)[0].value()
        except Exception as e:
            print(f"Error analyzing replicate {rep} in leg '{leg}': {e}")
            dg = float('nan')

        self.table.add_row(str(run_dir), str(f"{dg:.2f} kcal/mol"))
        self.console.print(self.table)

    def aggregate_leg(self, out_dir: Path, leg: str):
        pass

    def compare_edge(self, out_dir: Path):
        pass
