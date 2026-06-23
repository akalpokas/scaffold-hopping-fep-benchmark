"""
Orchestrates the Relative Binding Free Energy (RBFE) pipeline, generating
the alchemically merged bound (complex) and free (solvent) thermodynamic legs.
"""

import json
import argparse
from pathlib import Path
from typing import List
from pydantic import ValidationError

import BioSimSpace as BSS
import sire as sr

from pipeline.schemas import RBFEEdge
from pipeline import core

# ==========================================
# Setup Functions
# ==========================================


def build_rbfe_system(config: RBFEEdge):
    """
    Sets up a Relative Binding Free Energy (RBFE) calculation.
    Generates a merged alchemical ligand in both a free (solvent) and bound (complex) state.
    """
    print(f"[{config.id}] Building RBFE system...")

    # 1. Load and Parameterize Ligands
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

    # Extract changed bonds for diagnostic purposes
    sire_system = sr.convert.biosimspace_to_sire(merged)
    perturbable_mol = sire_system.molecules("molecule property is_perturbable")[0]
    pert = perturbable_mol.perturbation(map={"coordinates": "coordinates0"})
    pert_omm = pert.to_openmm()

    changed_bonds_df = pert_omm.changed_bonds()
    changed_bonds_df.to_csv(config.output_dir / "changed_bonds.csv", index=False)

    # 3. Setup and Save Free Leg
    solvated_free = core.solvate_and_neutralize(
        system_to_solvate=merged,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )

    free_path = str(config.output_dir / "free")
    BSS.Stream.save(solvated_free, free_path)

    # 4. Load Protein, Parameterize, and Save Bound Leg
    protein_files = [str(p) for p in config.protein_paths]
    protein_xtal = BSS.IO.readMolecules(protein_files)[0]

    protein, xtal_waters = core.parameterize_protein(protein_xtal, config.protein_ff)

    # Combine into holo system
    if xtal_waters is not None:
        protein_holo = protein + merged + xtal_waters
    else:
        protein_holo = protein + merged

    solvated_bound = core.solvate_and_neutralize(
        system_to_solvate=protein_holo,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )

    bound_path = str(config.output_dir / "bound")
    BSS.Stream.save(solvated_bound, bound_path)


# ==========================================
# Pipeline Orchestration
# ==========================================


def load_network(json_path: Path, specific_id: str = None) -> List[RBFEEdge]:
    """Loads the JSON and parses it strictly as RBFE edges."""
    with open(json_path, "r") as f:
        raw_data = json.load(f)

    if specific_id:
        raw_data = [d for d in raw_data if d.get("edge_id") == specific_id]

    parsed_configs = []
    for entry in raw_data:
        try:
            config = RBFEEdge(**entry)
            parsed_configs.append(config)
        except ValidationError as e:
            print(
                f"Skipping {entry.get('edge_id', 'unknown')}: Invalid RBFE schema.\n{e}"
            )

    return parsed_configs


def run_mapping_stage(network_path: Path, specific_id: str = None):
    configs = load_network(network_path, specific_id)
    print(f"Generating mapping visualizations for {len(configs)} edges...")

    for config in configs:
        config.create_output_dir()
        try:
            lig_a = BSS.IO.readMolecules([str(p) for p in config.ligand_a_paths])[0]
            lig_b = BSS.IO.readMolecules([str(p) for p in config.ligand_b_paths])[0]
            int_mapping = {int(k): int(v) for k, v in config.mapping.items()}

            core.visualize_mcs(lig_a, lig_b, int_mapping, output_dir=config.output_dir)
            print(f"Mapped: {config.id}")
        except Exception as e:
            print(f"Map Failed for {config.id}: {str(e)}")


def run_setup_stage(network_path: Path, specific_id: str = None):
    configs = load_network(network_path, specific_id)
    print(f"Building systems for {len(configs)} edges...")

    for config in configs:
        config.create_output_dir()
        try:
            build_rbfe_system(config)
            print(f"Successfully Built: {config.id}")
        except Exception as e:
            print(f"Setup Failed for {config.id}: {str(e)}")


# ==========================================
# Command Line Interface
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RBFE Preparation Pipeline")
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
