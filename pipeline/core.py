"""
Houses the shared BioSimSpace and Sire operations, providing
reusable functions for ligand parameterization, system solvation, and mapping.
"""

from __future__ import annotations
import traceback
from pathlib import Path
from typing import Optional, Tuple, Dict
from unittest.mock import patch

import BioSimSpace as BSS

# Import the Enums from previously defined schemas file
from .schemas import LigandForcefield, ProteinForcefield

# ==========================================
# Parameterization Utilities
# ==========================================


def parameterize_ligand(
    ligand: BSS.Types.Molecule, ff_type: LigandForcefield
) -> BSS.Types.Molecule:
    """
    Applies the specified forcefield to a ligand.
    """
    if ff_type == LigandForcefield.GAFF2:
        return BSS.Parameters.gaff2(ligand).getMolecule()

    elif ff_type == LigandForcefield.OPENFF:
        # use_nagl=False forces standard AM1-BCC charge calculation
        return BSS.Parameters.openff_unconstrained_2_2_1(
            ligand, use_nagl=False
        ).getMolecule()

    elif ff_type == LigandForcefield.PRE_PARAMETRIZED:
        return ligand

    else:
        raise ValueError(f"Unsupported ligand forcefield: '{ff_type.value}'")


def parameterize_protein(
    protein_xtal: BSS.Types.Molecule, ff_type: ProteinForcefield
) -> Tuple[BSS.Types.Molecule, Optional[BSS.Types.Molecule]]:
    """
    Extracts the protein and crystal waters from a complex, parameterizing them.
    Returns a tuple of (parameterized_protein, parameterized_waters).
    """
    # Extract protein (everything that is not water)
    protein = protein_xtal.extract(
        [atom.index() for atom in protein_xtal.search("not resname WAT").atoms()]
    )

    # Extract crystal waters
    try:
        xtal_waters = protein_xtal.extract(
            [atom.index() for atom in protein_xtal.search("resname WAT").atoms()]
        )
    except Exception as e:
        print(
            f"Warning: Failed to extract crystal waters. Proceeding without them. Error: {e}"
        )
        xtal_waters = None

    # Parameterize
    if ff_type == ProteinForcefield.AMBER14SB:
        protein = BSS.Parameters.ff14SB(protein, ensure_compatible=False).getMolecule()
        if xtal_waters is not None:
            xtal_waters = BSS.Parameters.ff14SB(
                xtal_waters, water_model="tip3p", ensure_compatible=False
            ).getMolecule()

    elif ff_type == ProteinForcefield.PRE_PARAMETRIZED:
        pass

    else:
        raise ValueError(f"Unsupported protein forcefield: '{ff_type.value}'")

    return protein, xtal_waters


# ==========================================
# Solvation & Box Setup
# ==========================================


def solvate_and_neutralize(
    system_to_solvate: BSS.Types.System, padding_nm: float, ion_conc: float
) -> BSS.Types.System:
    """
    Calculates the bounding box, applies padding, and solvates a system
    in a truncated octahedron, adding neutralizing ions.
    """
    # Initialize base shape to extract the correct angles for a truncated octahedron
    _, angles = BSS.Box.truncatedOctahedron(10 * BSS.Units.Length.nanometer)

    # Calculate bounding box
    box_min, box_max = system_to_solvate.getAxisAlignedBoundingBox()
    box_size = [y - x for x, y in zip(box_min, box_max)]

    # Apply padding
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


# ==========================================
# Visualization Utilities
# ==========================================


def visualize_mcs(
    mol_a: BSS.Types.Molecule,
    mol_b: BSS.Types.Molecule,
    mapping: Dict[int, int],
    output_dir: Path,
    filename: str = "mapping_vis.png",
):
    """
    Generates and saves the BioSimSpace 2D mapping visualization invisibly.
    """
    # Store the original notebook state and force it to True
    original_state = getattr(BSS, "_is_notebook", False)
    BSS._is_notebook = True

    try:
        # Intercept the IPython display call
        with patch("IPython.display.display") as mock_display:
            BSS.Align.viewMapping(mol_a, mol_b, mapping=mapping, pixels=1000)

            # Extract the captured image
            if mock_display.called:
                img_obj = mock_display.call_args[0][0]
                out_path = output_dir / filename
                with open(out_path, "wb") as f:
                    f.write(img_obj.data)
            else:
                raise RuntimeError("IPython display() was never called.")
    finally:
        # Clean up by restoring the original state
        BSS._is_notebook = original_state
