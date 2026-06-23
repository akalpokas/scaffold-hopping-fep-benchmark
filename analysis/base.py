from pathlib import Path


class AnalysisModule:
    """Base class for all analysis plugins."""

    @property
    def name(self) -> str:
        raise NotImplementedError

    def analyze_replicate(self, run_dir: Path, out_dir: Path, leg: str, rep: int):
        """1. Analyzes a single trajectory (e.g., calculates one RMSF array)."""
        pass

    def aggregate_leg(self, out_dir: Path, leg: str):
        """2. Aggregates the 3 replicates for a single leg (e.g., mean/std)."""
        pass

    def compare_edge(self, out_dir: Path):
        """3. Compares Bound vs Free for the entire edge (e.g., delta RMSF)."""
        pass
