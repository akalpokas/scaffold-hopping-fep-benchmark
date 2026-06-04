import json
import argparse
from pathlib import Path

# Import analysis modules
from analysis_modules.ddG import ddGAnalyzer
from analysis_modules.energy_traj import EnergyTrajAnalyzer

# Registry of available analyses
AVAILABLE_MODULES = {
    "energy_traj": EnergyTrajAnalyzer(),
    "ddG": ddGAnalyzer(),
}

def run_analysis(network_path: Path, active_modules: list, protocol: str, k: str, de: str):
    with open(network_path, 'r') as f:
        edges_data = json.load(f)

    print(f"Running analyses: {[m.name for m in active_modules]}\n")

    # Determine which replicates to run
    replicates_to_run = [args.replicate] if args.replicate is not None else [1, 2, 3]

    # Determine which legs to run
    legs_to_run = [args.leg_name.lower()] if args.leg_name.lower() in ["free", "bound"] else ["free", "bound"]

    for edge in edges_data:
        edge_id = edge["edge_id"]
        edge_dir = Path(edge["output_dir"])
        
        # Centralized analysis output folder
        out_dir = edge_dir / "analysis"
        out_dir.mkdir(exist_ok=True)
        
        print(f"=== Analyzing Edge: {edge_id} ===")

        # 1 & 2: Replicate and Leg Level
        for leg in legs_to_run:
            for rep in replicates_to_run:
                run_folder = f"{leg}_k_{k}_de_{de}_{protocol}_protocol_repl_{rep}"
                run_dir = edge_dir / run_folder
                
                if run_dir.exists():
                    # Pass the directory to every active module
                    for module in active_modules:
                        module.analyze_replicate(run_dir, out_dir, leg, rep)
                else:
                    print(f"    [!] Missing run directory: {run_folder}")
            
            # After all replicates for this leg are done, trigger aggregation
            for module in active_modules:
                module.aggregate_leg(out_dir, leg)

        # 3: Edge Level (Cross-leg comparison)
        for module in active_modules:
            module.compare_edge(out_dir)
            
        print("") # Formatting newline

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extensible Analysis Router")
    parser.add_argument("--network", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=str)
    
    # Let the user choose which analyses to run
    parser.add_argument("--modules", nargs="+", choices=AVAILABLE_MODULES.keys(), default=["energy_traj"],
                        help="List of analysis modules to run (e.g., --modules energy_traj convergence)")
    
    parser.add_argument("--k", default="125", type=str)
    parser.add_argument("--de", default="150", type=str)
    parser.add_argument("--leg_name", required=False, default="both", type=str, choices=["free", "bound", "both"],
                        help="Leg to run: 'free', 'bound', or 'both' (default)")
    
    parser.add_argument("--replicate", required=False, type=int, help="Optional: Specific replicate number to run (default: runs 1, 2, and 3)")
    
    args = parser.parse_args()
    
    # Instantiate the requested modules
    active_modules = [AVAILABLE_MODULES[m] for m in args.modules]
    
    run_analysis(args.network, active_modules, args.protocol, args.k, args.de)