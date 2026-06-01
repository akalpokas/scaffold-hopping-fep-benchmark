import json
import argparse
import shutil
import sys
from pathlib import Path


def clean_runs(
    network_path: Path, protocol: str, specific_edge: str = None, force: bool = False
):
    with open(network_path, "r") as f:
        edges_data = json.load(f)

    if specific_edge:
        edges_data = [e for e in edges_data if e.get("edge_id") == specific_edge]

    directories_to_delete = []

    print(f"Scanning for '{protocol}' runs to clean...")

    # 1: Collect Targets
    for edge in edges_data:
        edge_dir = Path(edge["output_dir"])

        if not edge_dir.exists():
            continue

        # Look for any DIRECTORY inside the edge folder that contains the protocol name
        for item in edge_dir.iterdir():
            if item.is_dir() and protocol in item.name:
                directories_to_delete.append(item)

    if not directories_to_delete:
        print("No directories found matching the criteria. Nothing to clean.")
        return

    # 2: Display and Confirm
    print(f"\nFound {len(directories_to_delete)} directories slated for removal:")
    for d in directories_to_delete:
        print(f"  - {d}")

    if not force:
        # Prompt the user if --force was NOT provided
        response = (
            input(
                f"\nAre you sure you want to permanently delete these {len(directories_to_delete)} directories? [y/N]: "
            )
            .strip()
            .lower()
        )
        if response not in ["y", "yes"]:
            print("Cleanup aborted by user. No files were modified.")
            sys.exit(0)
        print("Proceeding with cleanup...")
    else:
        print("\n[!] --force flag detected. Bypassing user confirmation.")

    # 3: Execute Deletion
    deleted_count = 0
    for item in directories_to_delete:
        try:
            shutil.rmtree(item)
            deleted_count += 1
        except Exception as e:
            print(f"Failed to delete {item}. Error: {e}")

    print(
        f"\nCleanup complete. Successfully removed {deleted_count}/{len(directories_to_delete)} directories."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Clean specific MD runs from an RBFE network"
    )

    # Required arguments
    parser.add_argument(
        "--network", required=True, type=Path, help="Path to network.json"
    )
    parser.add_argument(
        "--protocol",
        required=True,
        type=str,
        help="Protocol name to target for deletion (e.g., 'testing')",
    )

    # Optional arguments
    parser.add_argument(
        "--edge", default=None, type=str, help="Optional: Clean only a single edge"
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Skip user confirmation prompt and delete immediately",
    )

    args = parser.parse_args()

    clean_runs(
        network_path=args.network,
        protocol=args.protocol,
        specific_edge=args.edge,
        force=args.force,
    )
