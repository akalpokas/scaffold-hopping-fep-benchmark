import pandas as pd
from pathlib import Path
from .base import AnalysisModule
from rich.console import Console
from rich.table import Table

class EnergyTrajAnalyzer(AnalysisModule):

    @property
    def name(self) -> str:
        return "energy_traj"

    def analyze_replicate(self, run_dir: Path, out_dir: Path, leg: str, rep: int, T: float = 300.0, **kwargs):
        """
        Reports the energy trajectories from the simulations.

        Parameters:
        - run_dir (Path): The directory containing the parquet files.
        - T (float): The temperature at which the simulations were performed. Default is 300.0.

        Returns:
        - lengths (list): A list of lengths of dataframes.
        """
        from glob import glob
        from pathlib import Path
        from BioSimSpace.FreeEnergy import Relative

        self.console = Console()
        self.table = Table(title="Energy Trajectories", show_lines=True)
        self.table.add_column("Path", style="cyan", no_wrap=True)
        self.table.add_column("Energy Trajectory Length", style="magenta")

        files = glob(f"{run_dir}/*.parquet")

        glob_path = Path(run_dir)
        

        analysed_df_list = []


        files = sorted(glob_path.glob("**/*.parquet"))
        for f in files:
            path = Path(f)
            analysed_df = Relative._somd2_extract(path, T=T, **kwargs)
            analysed_df_list.append(analysed_df)

        # Populate the table with the results
        # check if all dataframes contain the same number of entries
        lengths = [len(df) for df in analysed_df_list]
        if len(set(lengths)) != 1:
            self.console.print(f"Not all dataframes contain the same number of entries: {lengths}")

        self.table.add_row(str(run_dir), str(lengths))
        self.console.print(self.table)

    def aggregate_leg(self, out_dir: Path, leg: str):
        pass

    def compare_edge(self, out_dir: Path):
        pass

