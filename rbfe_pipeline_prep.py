import json
import argparse
import traceback
from unittest.mock import patch
from pathlib import Path
from enum import Enum
from typing import Any, Dict, Optional
from pydantic import BaseModel, FilePath, Field
import BioSimSpace as BSS
import sire as sr

# ==========================================
# Schema Definitions
# ==========================================


class LigandForcefield(str, Enum):
    OPENFF = "openff"
    GAFF2 = "gaff2"


class ProteinForcefield(str, Enum):
    AMBER14SB = "amber14"


class RBFEEdge(BaseModel):
    edge_id: str = Field(
        ...,
        description="Unique identifier for this transformation, e.g., 'lig1_to_lig2'",
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict,
        description="Optional tracking info like experimental ddG, Lomap scores, etc.",
    )

    # FilePath ensures the files actually exist before the pipeline does anything
    ligand_a_path: FilePath
    ligand_b_path: FilePath
    protein_path: FilePath

    mapping: Dict[int, int] = Field(
        ..., description="Atom index mapping from Ligand A to Ligand B"
    )

    ligand_ff: LigandForcefield = LigandForcefield.GAFF2
    protein_ff: ProteinForcefield = ProteinForcefield.AMBER14SB

    solvent_padding_nm: float = Field(
        1.5, ge=0.0, description="Padding for solvent box in nanometers"
    )
    ionic_strength_molar: float = Field(
        0.15, ge=0.0, description="Ionic strength in molar for solvation"
    )

    output_dir: Path = Field(..., description="Directory to save the outputs")

    def create_output_dir(self):
        """Utility to ensure the working directory exists."""
        self.output_dir.mkdir(parents=True, exist_ok=True)


# ==========================================
# Core Processing Functions
# ==========================================


def _solvate_and_neutralize(system_to_solvate, padding_nm: float, ion_conc: float):
    """
    Helper function to calculate the bounding box, apply padding,
    and solvate a system in a truncated octahedron.
    """
    # Initialize base shape to extract the correct angles for a truncated octahedron
    _, angles = BSS.Box.truncatedOctahedron(10 * BSS.Units.Length.nanometer)

    # Calculate bounding box
    box_min, box_max = system_to_solvate.getAxisAlignedBoundingBox()
    box_size = [y - x for x, y in zip(box_min, box_max)]

    # Apply padding (using the schema's parameter converted to BSS nanometers)
    padding = padding_nm * BSS.Units.Length.nanometer
    box_length = max(box_size) + 2 * padding

    # Solvate
    solvated = BSS.Solvent.tip3p(
        molecule=system_to_solvate,
        box=3 * [box_length],
        angles=angles,
        ion_conc=ion_conc,
        is_neutral=True,
    )
    solvated.reduceBoxVectors()

    return solvated


def visualize_mcs(config: "RBFEEdge"):
    """
    Generates and saves the BioSimSpace 2D mapping visualization.
    """
    mol_a = BSS.IO.readMolecules(str(config.ligand_a_path))[0]
    mol_b = BSS.IO.readMolecules(str(config.ligand_b_path))[0]

    mapping = config.mapping

    # Store the original notebook state and force it to True
    original_state = getattr(BSS, "_is_notebook", False)
    BSS._is_notebook = True

    try:
        # Intercept the IPython display call
        # patch() replaces the display function with a 'Mock' object that records
        # what was passed into it instead of actually printing to the screen.
        with patch("IPython.display.display") as mock_display:

            BSS.Align.viewMapping(mol_a, mol_b, mapping=mapping, pixels=1000)

            # Extract the captured image
            if mock_display.called:
                # mock_display.call_args[0][0] gets the first argument of the first call,
                # which is the IPython 'Image' object BSS created.
                img_obj = mock_display.call_args[0][0]

                # The Image object stores raw file bytes in its .data attribute
                out_path = config.output_dir / "mapping_vis.png"
                with open(out_path, "wb") as f:
                    f.write(img_obj.data)
            else:
                raise RuntimeError(
                    f"Failed to capture image for {config.output_dir.name}: "
                    "IPython display() was never called."
                )
    finally:
        # Clean up by restoring the original state
        BSS._is_notebook = original_state


def setup_alchemical_system(config: "RBFEEdge"):
    """
    Parameterizes ligands and protein, merges the alchemical ligands,
    and generates both the free and bound legs of the RBFE transformation.
    """
    # ==========================================
    # Load and Parameterize Ligands
    # ==========================================
    lig_a = BSS.IO.readMolecules([str(config.ligand_a_path)])[0]
    lig_b = BSS.IO.readMolecules([str(config.ligand_b_path)])[0]

    # Map the schema Enum to the specific BSS function calls
    if config.ligand_ff.name == "GAFF2":
        lig_a = BSS.Parameters.gaff2(lig_a).getMolecule()
        lig_b = BSS.Parameters.gaff2(lig_b).getMolecule()
    else:
        lig_a = BSS.Parameters.openff_unconstrained_2_2_1(
            lig_a, use_nagl=False
        ).getMolecule()
        lig_b = BSS.Parameters.openff_unconstrained_2_2_1(
            lig_b, use_nagl=False
        ).getMolecule()

    # ==========================================
    # Align, Merge, and Extract Diagnostics
    # ==========================================
    # We must explicitly cast dict keys to integers in case they were loaded purely from JSON
    int_mapping = {int(k): int(v) for k, v in config.mapping.items()}

    lig_a_aligned = BSS.Align.rmsdAlign(lig_a, lig_b, int_mapping)
    transformation_type = str(config.metadata.get("notes", "unknown")).lower()
    if "bond" in transformation_type:
        merged = BSS.Align.merge(
            lig_a_aligned,
            lig_b,
            int_mapping,
            allow_ring_breaking=True,
            allow_ring_size_change=True,
            force=True
        )
    elif "standard" in transformation_type:
        merged = BSS.Align.merge(
            lig_a_aligned,
            lig_b,
            int_mapping,
        )
    else:
        raise ValueError(
            f"Unrecognized transformation type in metadata notes: '{transformation_type}'. "
            "Please include 'bond' or 'standard' in the notes to indicate the expected transformation."
        )

    # Diagnostic extraction using Sire
    # merged_system = merged.toSystem()
    merged_system = merged
    sire_system = sr.convert.biosimspace_to_sire(merged_system)

    mol = sire_system.molecules("molecule property is_perturbable")
    perturbable_mol = mol[0]
    pert = perturbable_mol.perturbation(map={"coordinates": "coordinates0"})
    pert_omm = pert.to_openmm()

    # Save the changed bonds diagnostic to the output directory immediately
    changed_bonds_df = pert_omm.changed_bonds()
    changed_bonds_df.to_csv(config.output_dir / "changed_bonds.csv", index=False)

    # ==========================================
    # Setup and Save Free Leg
    # ==========================================
    solvated_free = _solvate_and_neutralize(
        system_to_solvate=merged,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )

    # Stream.save automatically appends the necessary extensions (.prm7/.rst7)
    free_path = str(config.output_dir / "free")
    BSS.Stream.save(solvated_free, free_path)

    # ==========================================
    # Load Protein, Parameterize, and Save Bound Leg
    # ==========================================
    protein_xtal = BSS.IO.readMolecules([str(config.protein_path)])[0]

    protein = protein_xtal.extract(
        [atom.index() for atom in protein_xtal.search("not resname WAT").atoms()]
    )
    try:
        xtal_waters = protein_xtal.extract(
            [atom.index() for atom in protein_xtal.search("resname WAT").atoms()]
        )
    except Exception as e:
        print(
            f"Warning: Failed to extract crystal waters from {config.protein_path.name}. "
            "Proceeding without them. Error details: "
        )
        print(e)
        xtal_waters = None
    # Parameterize using the protein FF defined in your schema
    if config.protein_ff.name == "AMBER14SB":
        protein = BSS.Parameters.ff14SB(protein, ensure_compatible=False).getMolecule()
        if xtal_waters is not None:
            xtal_waters = BSS.Parameters.ff14SB(
                xtal_waters, water_model="tip3p", ensure_compatible=False
            ).getMolecule()
    elif config.protein_ff.name == "AMBER99SB_ILDN":
        protein = BSS.Parameters.ff99SBildn(
            protein, ensure_compatible=False
        ).getMolecule()
        if xtal_waters is not None:
            xtal_waters = BSS.Parameters.ff99SBildn(
                xtal_waters, water_model="tip3p", ensure_compatible=False
            ).getMolecule()

    # Combine into holo system
    # protein_holo = protein.toSystem() + merged + xtal_waters
    if xtal_waters is not None:
        protein_holo = protein + merged + xtal_waters
    else:
        protein_holo = protein + merged

    solvated_bound = _solvate_and_neutralize(
        system_to_solvate=protein_holo,
        padding_nm=config.solvent_padding_nm,
        ion_conc=config.ionic_strength_molar,
    )

    bound_path = str(config.output_dir / "bound")
    BSS.Stream.save(solvated_bound, bound_path)


# ==========================================
# Pipeline Orchestration
# ==========================================


def load_network(json_path: Path, specific_edge_id: Optional[str] = None):
    with open(json_path, "r") as f:
        edges_data = json.load(f)

    if specific_edge_id:
        edges_data = [e for e in edges_data if e.get("edge_id") == specific_edge_id]

    return edges_data


def run_mapping_stage(network_path: Path, specific_edge_id: str = None):
    edges_data = load_network(network_path, specific_edge_id)
    print(
        f"Generating mapping visualizations for {len(edges_data)} edges in {network_path.name}..."
    )

    for edge_dict in edges_data:
        edge_name = Path(edge_dict.get("output_dir", "unknown")).name
        try:
            config = RBFEEdge(**edge_dict)
            config.create_output_dir()
            visualize_mcs(config)
            print(f"Mapped: {edge_name}")
        except Exception as e:
            print(f"Map Failed: {edge_name} - {str(e)}")


def run_setup_stage(network_path: Path, specific_edge_id: str = None):
    edges_data = load_network(network_path, specific_edge_id)
    print(f"Building systems for {len(edges_data)} edges in {network_path.name}...")

    for edge_dict in edges_data:
        edge_name = Path(edge_dict.get("output_dir", "unknown")).name
        try:
            config = RBFEEdge(**edge_dict)
            config.create_output_dir()
            setup_alchemical_system(config)
            print(f"System Built: {edge_name}")
        except Exception as e:
            print(f"Setup Failed: {edge_name} - {str(e)}")


# ==========================================
# Command Line Interface
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RBFE Benchmark Preparation Pipeline")

    # Positional argument: map or setup
    parser.add_argument(
        "stage", choices=["map", "setup"], help="Which stage of the pipeline to run"
    )

    # Required network file
    parser.add_argument(
        "--network",
        required=True,
        type=Path,
        help="Path to the JSON network file to process",
    )

    # Optional edge filter
    parser.add_argument(
        "--edge-id",
        default=None,
        help="Optional: Run a specific edge by its unique ID (e.g., 'lig1_to_lig2')",
    )

    args = parser.parse_args()

    if args.stage == "map":
        run_mapping_stage(args.network, args.edge_id)
    elif args.stage == "setup":
        run_setup_stage(args.network, args.edge_id)
