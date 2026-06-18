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

def run_analysis(network_path: Path, active_modules: list, protocol: str, k: str, de: str, modifiers: str, specific_edge: str = None):
    with open(network_path, 'r') as f:
        edges_data = json.load(f)

    # Filter for a specific edge if provided
    if specific_edge:
        edges_data = [e for e in edges_data if e.get("edge_id") == specific_edge]
        if not edges_data:
            print(f"Error: Edge '{specific_edge}' not found in network.")
            return

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
                run_folder = f"{leg}_k_{k}_{modifiers}_de_{de}_{protocol}_protocol_repl_{rep}"
                run_dir = edge_dir / run_folder
                
                if run_dir.exists():
                    # Pass the directory to every active module
                    for module in active_modules:
                        module.analyze_replicate(edge_id, run_dir, out_dir, leg, rep)
                else:
                    print(f"    [!] Missing run directory: {run_folder}")
            
            # After all replicates for this leg are done, trigger aggregation
            for module in active_modules:
                module.aggregate_leg(edge_id, out_dir, leg)

        # 3: Edge Level (Cross-leg comparison)
        for module in active_modules:
            module.compare_edge(edge_id, out_dir)
        print("") # Formatting newline

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Network-wide analysis script for RBFE benchmarks. Run folders should be organized as: <edge_output_dir>/<leg>_k_<k>_<modifiers>_de_<de>_<protocol>_protocol_repl_<replicate>/")
    
    # Required arguments
    parser.add_argument("--network", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=str)
    
    # Let the user choose which analyses to run
    parser.add_argument("--modules", nargs="+", choices=AVAILABLE_MODULES.keys(), default=["energy_traj"],
                        help="List of analysis modules to run (e.g., --modules energy_traj convergence)")
    
    # Optional arguments with smart defaults
    parser.add_argument("--edge-id", default=None, type=str, help="Optional: Analyze only a single edge")
    
    parser.add_argument("--k", default="125", type=str, help="Soft Morse bond strength parameter for the run. Default is 125 kcal/mol/A^2.")
    parser.add_argument("--de", default="50", type=str, help="Soft Morse dissociation energy parameter for the run. Default is 50 kcal/mol.")
    
    parser.add_argument("--leg_name", required=False, default="both", type=str, choices=["free", "bound", "both"],
                        help="Leg to run: 'free', 'bound', or 'both' (default)")
    
    parser.add_argument("--replicate", required=False, type=int, help="Optional: Specific replicate number to run (default: runs 1, 2, and 3)")

    parser.add_argument(
        "--modifiers",
        type=str,
        required=False,
        default="",
        help="Optional string of modifiers to include in the run folder name (e.g., 'ghost' if --ghost_modifications is set)"
    )
    
    args = parser.parse_args()
    
    # Instantiate the requested modules
    active_modules = [AVAILABLE_MODULES[m] for m in args.modules]
    
    run_analysis(
        network_path=args.network, 
        active_modules=active_modules, 
        protocol=args.protocol, 
        k=args.k, 
        de=args.de, 
        modifiers=args.modifiers,
        specific_edge=args.edge_id
    )