import json
import argparse
import subprocess
import time
from pathlib import Path


def submit_jobs(
    network_path: Path,
    specific_edge: str = None,
    protocol: str = "prod",
    leg_name: str = "both",
    replicate: int = None,
    bond_strength: int = 125,
    de_strength: int = 150,
    ghost_modifications: bool = False,
):

    Path("slurm_logs").mkdir(exist_ok=True)
    with open(network_path, "r") as f:
        edges_data = json.load(f)

    if specific_edge:
        edges_data = [e for e in edges_data if e.get("edge_id") == specific_edge]
        if not edges_data:
            print(f"Error: Edge '{specific_edge}' not found in network.")
            return

    # Determine which replicates to run
    replicates_to_run = [replicate] if replicate is not None else [1, 2, 3]

    # Determine which legs to run
    legs_to_run = (
        [leg_name.lower()]
        if leg_name.lower() in ["free", "bound"]
        else ["free", "bound"]
    )

    total_jobs = len(edges_data) * len(legs_to_run) * len(replicates_to_run)
    print(f"Preparing to submit {total_jobs} total jobs to SLURM...")
    print(
        f"Edges: {len(edges_data)} | Legs: {legs_to_run} | Replicates: {replicates_to_run}\n"
    )

    for edge in edges_data:
        edge_id = edge["edge_id"]

        # Loop through the assigned legs
        for leg in legs_to_run:

            # Loop through the assigned replicates
            for rep in replicates_to_run:

                # Build the sbatch command
                sbatch_cmd = [
                    "sbatch",
                    f"--job-name=rbfe_{edge_id}_{leg}_rep{rep}",
                    f"--output=slurm_logs/{edge_id}_{leg}_rep{rep}.slurm.out",
                    "slurm_run_apptainer_prod.sh",
                    str(network_path),
                    str(edge_id),
                    str(protocol),
                    str(leg),
                    str(bond_strength),
                    str(de_strength),
                    str(rep),
                    str(ghost_modifications),
                ]

                try:
                    # Run the command and capture the output
                    result = subprocess.run(
                        sbatch_cmd, check=True, capture_output=True, text=True
                    )
                    print(
                        f"Submitted {edge_id} | {leg.upper()} | Rep {rep} -> {result.stdout.strip()}"
                    )

                    # Brief pause to avoid overwhelming the SLURM scheduler daemon
                    time.sleep(0.2)

                except subprocess.CalledProcessError as e:
                    print(
                        f"Failed {edge_id} | {leg.upper()} | Rep {rep} -> Error: {e.stderr}"
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit RBFE network to SLURM")

    # Required arguments
    parser.add_argument(
        "--network", required=True, type=Path, help="Path to network.json"
    )
    parser.add_argument(
        "--protocol",
        required=True,
        type=str,
        help="Protocol to run (e.g. 'testing', 'prod', 'prod_2fs',  'long')",
    )

    # Optional arguments with smart defaults
    parser.add_argument(
        "--edge-id", default=None, type=str, help="Optional: Submit only a single edge"
    )

    parser.add_argument(
        "--leg_name",
        default="both",
        type=str,
        choices=["free", "bound", "both"],
        help="Leg to run: 'free', 'bound', or 'both' (default)",
    )

    parser.add_argument(
        "--replicate",
        required=False,
        type=int,
        help="Optional: Specific replicate number to run (default: runs 1, 2, and 3)",
    )
    parser.add_argument(
        "--bond_strength",
        type=int,
        default=125,
        help="Soft Morse bond strength parameter for the run. Default is 125 kcal/mol/A^2.",
    )
    parser.add_argument(
        "--de_strength",
        type=int,
        default=50,
        help="Soft Morse dissociation energy parameter for the run. Default is 50 kcal/mol.",
    )

    parser.add_argument(
        "--ghost_modifications",
        action="store_true",
        help="Whether to apply ghost modifications to the system. If not provided, ghost modifications will not be applied.",
    )

    args = parser.parse_args()

    submit_jobs(
        network_path=args.network,
        specific_edge=args.edge_id,
        protocol=args.protocol,
        leg_name=args.leg_name,
        replicate=args.replicate,
        bond_strength=args.bond_strength,
        de_strength=args.de_strength,
        ghost_modifications=args.ghost_modifications,
    )
