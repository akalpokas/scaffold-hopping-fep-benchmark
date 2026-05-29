import json
import argparse
import subprocess
import time
from pathlib import Path

def submit_jobs(network_path: Path, specific_edge: str = None, protocol: str = "prod", leg_name: str = "free", replicate: int = None, de_strength: int = 150):
    Path("logs").mkdir(exist_ok=True)
    with open(network_path, 'r') as f:
        edges_data = json.load(f)
        
    if specific_edge:
        edges_data = [e for e in edges_data if e.get("edge_id") == specific_edge]
        if not edges_data:
            print(f"Error: Edge '{specific_edge}' not found in network.")
            return

    # Determine which replicates to run
    replicates_to_run = [replicate] if replicate is not None else [1, 2, 3]

    print(f"Preparing to submit {len(edges_data)} edges for replicates {replicates_to_run} to SLURM...")

    for edge in edges_data:
        edge_id = edge["edge_id"]
        
        # Inner loop handles the replicates for each edge
        for rep in replicates_to_run:
            
            # Build the sbatch command
            sbatch_cmd = [
                "sbatch",
                f"--job-name=rbfe_{edge_id}_rep{rep}",
                f"--output=logs/{edge_id}_rep{rep}.slurm.out",
                "slurm_run_singularity_prod.sh",    
                str(network_path),
                str(edge_id),             
                str(protocol),            
                str(leg_name),            
                str(de_strength),         
                str(rep)           
            ]
            
            try:
                # Run the command and capture the output
                result = subprocess.run(sbatch_cmd, check=True, capture_output=True, text=True)
                print(f"Submitted {edge_id} (Replicate {rep}): {result.stdout.strip()}")
                
                # Brief pause to avoid overwhelming the SLURM scheduler daemon
                time.sleep(0.2)
                
            except subprocess.CalledProcessError as e:
                print(f"Failed to submit {edge_id} (Replicate {rep}). Error: {e.stderr}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit RBFE network to SLURM")
    parser.add_argument("--network", required=True, type=Path, help="Path to network.json")
    parser.add_argument("--edge", default=None, type=str, help="Optional: Submit only a single edge")
    parser.add_argument("--protocol", required=True, type=str, help="Protocol to run (e.g. 'testing', 'prod', 'long')")
    parser.add_argument("--leg_name", required=True, type=str, help="Leg to run (e.g. 'free' or 'bound')")
    parser.add_argument("--replicate", required=False, type=int, help="Replicate number to run (e.g. 1)")
    parser.add_argument("--de_strength", default=150, type=int, help="DE strength parameter for the run")
    
    args = parser.parse_args()
    
    submit_jobs(
        network_path=args.network, 
        specific_edge=args.edge, 
        protocol=args.protocol, 
        leg_name=args.leg_name, 
        replicate=args.replicate,
        de_strength=args.de_strength
    )