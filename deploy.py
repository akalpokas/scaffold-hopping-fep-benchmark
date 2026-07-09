import json
import argparse
import subprocess
import time
from pathlib import Path
from pipeline._utils import validate_protocol


def determine_legs(edge_dict: dict, requested_leg: str) -> list[str]:
    """
    Determines which thermodynamic legs to run based on the JSON schema
    and the user's requested leg.
    """
    # If protein paths exist, it's RBFE. Otherwise, it's RHFE/AHFE.
    if "protein_paths" in edge_dict:
        valid_legs = ["free", "bound"]
    else:
        # If there is no ligand_b_paths, it's an AHFE node. Otherwise, it's an RHFE edge.
        if "ligand_b_paths" in edge_dict:
            valid_legs = ["solvent", "vacuum"]
        else:
            valid_legs = ["solvent"]

    if requested_leg.lower() == "all":
        return valid_legs
    elif requested_leg.lower() in valid_legs:
        return [requested_leg.lower()]
    else:
        # If the user requested 'bound' but this is an RHFE edge, return empty
        return []


def submit_jobs(
    network_path: Path,
    specific_edge: str = None,
    protocol: str = "prod",
    leg_name: str = "all",
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
            print(f"Error: ID '{specific_edge}' not found in network.")
            return

    replicates_to_run = [replicate] if replicate is not None else [1, 2, 3]

    jobs_submitted = 0
    print(f"Preparing to submit jobs to SLURM...")

    for edge in edges_data:
        edge_id = edge.get("edge_id", edge.get("node_id", "unknown"))

        # Dynamically determine legs for this specific edge/node
        legs_to_run = determine_legs(edge, leg_name)

        if not legs_to_run:
            print(
                f"Skipping {edge_id}: Requested leg '{leg_name}' is not valid for this system type."
            )
            continue

        for leg in legs_to_run:
            for rep in replicates_to_run:
                # Build the sbatch command
                sbatch_cmd = [
                    "sbatch",
                    f"--job-name=fep_{edge_id}_{leg}_rep{rep}",
                    f"--output=slurm_logs/{edge_id}_{leg}_rep{rep}.slurm.out",
                    "slurm_run_apptainer_prod.sh",
                    str(network_path),
                    str(edge_id),
                    str(protocol),
                    str(leg),
                    str(bond_strength),
                    str(de_strength),
                    str(rep),
                    str(ghost_modifications),  # This is argument $8
                ]

                try:
                    result = subprocess.run(
                        sbatch_cmd, check=True, capture_output=True, text=True
                    )
                    print(
                        f"Submitted {edge_id} | {leg.upper()} | Rep {rep} -> {result.stdout.strip()}"
                    )
                    jobs_submitted += 1
                    time.sleep(0.1)  # Brief pause for the SLURM daemon

                except subprocess.CalledProcessError as e:
                    print(
                        f"Failed {edge_id} | {leg.upper()} | Rep {rep} -> Error: {e.stderr}"
                    )

    print(f"\nTotal jobs successfully submitted: {jobs_submitted}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit FEP network to SLURM")

    parser.add_argument(
        "--network", required=True, type=Path, help="Path to network.json"
    )
    parser.add_argument(
        "--protocol",
        required=True,
        type=validate_protocol,
        help="Protocol to run (e.g. 'prod', 'tucker_long', 'prod_rest2_2fs')",
    )
    parser.add_argument(
        "--edge-id",
        default=None,
        type=str,
        help="Optional: Submit only a single edge/node",
    )
    parser.add_argument(
        "--leg_name",
        default="all",
        type=str,
        choices=["free", "bound", "solvent", "vacuum", "all"],
        help="Leg to run: 'free', 'bound', 'solvent', 'vacuum', or 'all' (default)",
    )
    parser.add_argument(
        "--replicate",
        required=False,
        type=int,
        help="Specific replicate number to run (default: 1, 2, 3)",
    )
    parser.add_argument(
        "--bond_strength",
        type=int,
        default=125,
        help="Soft Morse bond strength. Default 125.",
    )
    parser.add_argument(
        "--de_strength",
        type=int,
        default=50,
        help="Soft Morse dissociation energy. Default 50.",
    )
    parser.add_argument(
        "--ghost_modifications", action="store_true", help="Apply ghost modifications."
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
