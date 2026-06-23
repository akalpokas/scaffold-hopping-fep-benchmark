"""
Routes and orchestrates Hydration Free Energy pipelines, generating the solvent
and vacuum thermodynamic legs for both relative (RHFE) and absolute (AHFE) inputs.
"""

import json
import argparse
from pathlib import Path
from typing import Union
from pydantic import ValidationError

import BioSimSpace as BSS
import sire as sr

from pipeline.schemas import RHFEEdge, AHFENode
from pipeline import core

# ==========================================
# Setup Functions
# ==========================================


def build_relative_system(config: RHFEEdge):
    """
    Sets up a Relative Hydration Free Energy (RHFE) calculation.
    Generates a merged alchemical ligand in both vacuum and solvent.
    """
    print(f"[{config.id}] Detected RHFE Edge.")

    # 1. Load and Parameterize
    lig_a_files = [str(p) for p in config.ligand_a_paths]
    lig_b_files = [str(p) for p in config.ligand_b_paths]

    lig_a = BSS.IO.readMolecules(lig_a_files)[0]
    lig_b = BSS.IO.readMolecules(lig_b_files)[0]

    lig_a = core.parameterize_ligand(lig_a, config.ligand_ff)
    lig_b = core.parameterize_ligand(lig_b, config.ligand_ff)

    # 2. Align and Merge
    int_mapping = {int(k): int(v) for k, v in config.mapping.items()}
    lig_a_aligned = BSS.Align.rmsdAlign(lig_a, lig_b, int_mapping)

    transformation_type = str(config.metadata.get("notes", "unknown")).lower()
    allow_ring_breaking = "bond" in transformation_type

    merged = BSS.Align.merge(
        lig_a_aligned,
        lig_b,
        int_mapping,
        allow_ring_breaking=allow_ring_breaking,
        allow_ring_size_change=allow_ring_breaking,
        force=allow_ring_breaking,
    )

    # 3. Setup and Save Vacuum Leg
    vacuum_path = str(config.output_dir / "vacuum")
    BSS.Stream.save(merged.toSystem(), vacuum_path)

    # 4. Setup and Save Solvent Leg
    solvated_system = core.solvate_and_neutralize(
        system_to_solvate=merged,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )
    solvent_path = str(config.output_dir / "solvent")
    BSS.Stream.save(solvated_system, solvent_path)


def build_absolute_system(config: AHFENode):
    """
    Sets up an Absolute Hydration Free Energy (AHFE) calculation.
    Decouples a single ligand into dummy atoms in both vacuum and solvent.
    """
    print(f"[{config.id}] Detected AHFE Node.")

    # 1. Load and Parameterize
    lig_a_files = [str(p) for p in config.ligand_a_paths]
    lig_a = BSS.IO.readMolecules(lig_a_files)[0]
    lig_a = core.parameterize_ligand(lig_a, config.ligand_ff)

    # Helper function to apply the Sire decoupling logic
    def _apply_decoupling(bss_system, output_path: Path):
        sire_sys = sr.convert.biosimspace_to_sire(bss_system)

        # In a pure system/solvated system, find the ligand
        # We assume the ligand is the first molecule, or we find it by size/properties
        mols = sire_sys.molecules()
        ligand_mol = mols[0]  # Usually the ligand is at index 0 before water

        perturbed_mol = sr.morph.decouple(ligand_mol, as_new_molecule=False)
        sire_sys.update(perturbed_mol)

        sr.stream.save(sire_sys, str(output_path))

    solvated_system = core.solvate_and_neutralize(
        system_to_solvate=lig_a,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )
    solvent_path = config.output_dir / "solvent.bss"
    _apply_decoupling(solvated_system, solvent_path)


# ==========================================
# Pipeline Orchestration
# ==========================================


def load_and_route_network(
    json_path: Path, specific_id: str = None
) -> list[Union[RHFEEdge, AHFENode]]:
    """Loads the JSON and categorizes each entry as Relative or Absolute."""
    with open(json_path, "r") as f:
        raw_data = json.load(f)

    if specific_id:
        raw_data = [d for d in raw_data if d.get("edge_id") == specific_id]

    parsed_configs = []
    for entry in raw_data:
        # Try parsing as a Relative Edge first
        try:
            config = RHFEEdge(**entry)
            parsed_configs.append(config)
            continue
        except ValidationError:
            pass

        # If it fails, try parsing as an Absolute Node
        try:
            config = AHFENode(**entry)
            parsed_configs.append(config)
        except ValidationError as e:
            print(f"Skipping {entry.get('edge_id', 'unknown')}: Invalid schema.\n{e}")

    return parsed_configs


def run_mapping_stage(network_path: Path, specific_id: str = None):
    configs = load_and_route_network(network_path, specific_id)

    for config in configs:
        config.create_output_dir()

        # Only RHFE edges have a mapping to visualize
        if isinstance(config, RHFEEdge):
            try:
                lig_a = BSS.IO.readMolecules([str(p) for p in config.ligand_a_paths])[0]
                lig_b = BSS.IO.readMolecules([str(p) for p in config.ligand_b_paths])[0]

                # Convert string keys to ints for BSS
                int_mapping = {int(k): int(v) for k, v in config.mapping.items()}

                core.visualize_mcs(
                    lig_a, lig_b, int_mapping, output_dir=config.output_dir
                )
                print(f"Mapped: {config.id}")
            except Exception as e:
                print(f"Map Failed for {config.id}: {str(e)}")
        else:
            print(
                f"Skipping Mapping for {config.id} (AHFE nodes do not require mapping)."
            )


def run_setup_stage(network_path: Path, specific_id: str = None):
    configs = load_and_route_network(network_path, specific_id)

    for config in configs:
        config.create_output_dir()
        try:
            if isinstance(config, RHFEEdge):
                build_relative_system(config)
            elif isinstance(config, AHFENode):
                build_absolute_system(config)

            print(f"Successfully Built: {config.id}")
        except Exception as e:
            print(f"Setup Failed for {config.id}: {str(e)}")


# ==========================================
# Command Line Interface
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HFE Pipeline Router (RHFE & AHFE)")
    parser.add_argument("stage", choices=["map", "setup"], help="Pipeline stage to run")
    parser.add_argument(
        "--network", required=True, type=Path, help="Path to JSON network"
    )
    parser.add_argument("--edge-id", default=None, help="Run a specific ID")

    args = parser.parse_args()

    if args.stage == "map":
        run_mapping_stage(args.network, args.edge_id)
    elif args.stage == "setup":
        run_setup_stage(args.network, args.edge_id)
