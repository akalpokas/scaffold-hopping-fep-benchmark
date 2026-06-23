"""
Defines the data structures, input validation rules, and JSON configuration
schemas for all thermodynamic cycles (RBFE, RHFE, AHFE) using Pydantic.
"""

from __future__ import annotations
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List
from pydantic import BaseModel, Field, FilePath, field_validator

# ==========================================
# Enums
# ==========================================


class LigandForcefield(str, Enum):
    OPENFF = "openff"
    GAFF2 = "gaff2"
    PRE_PARAMETRIZED = "pre_parametrized"


class ProteinForcefield(str, Enum):
    AMBER14SB = "amber14"
    PRE_PARAMETRIZED = "pre_parametrized"


# ==========================================
# Base Schema (Shared Configuration)
# ==========================================


class BaseConfig(BaseModel):
    """
    Core configuration shared across all thermodynamic cycles (RBFE, RHFE, AHFE).
    """

    id: str = Field(
        ...,
        alias="edge_id",
        description="Unique identifier (e.g., 'lig1_to_lig2' or 'lig1')",
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict,
        description="Optional tracking info like experimental ddG, Lomap scores, etc.",
    )

    ligand_a_paths: List[FilePath] = Field(
        ..., description="Paths to the topology and coordinate files for Ligand A"
    )
    ligand_ff: LigandForcefield = LigandForcefield.GAFF2

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
# Relative Calculations (RBFE & RHFE)
# ==========================================


class RelativeEdge(BaseConfig):
    """
    Shared schema for any relative alchemical transformation (requires two ligands).
    """

    ligand_b_paths: List[FilePath] = Field(
        ..., description="Paths to the topology and coordinate files for Ligand B"
    )
    mapping: Dict[int, int] = Field(
        ..., description="Atom index mapping from Ligand A to Ligand B"
    )

    @field_validator("metadata")
    @classmethod
    def validate_metadata_notes(cls, v: Dict[str, Any]) -> Dict[str, Any]:
        """Validates that the metadata dictates the morphing strategy."""
        if "notes" not in v:
            raise ValueError("The 'metadata' dictionary must contain a 'notes' key.")

        notes = str(v["notes"]).lower()
        valid_phrases = ["standard morph", "bond annihilation", "bond creation"]

        if not any(phrase in notes for phrase in valid_phrases):
            raise ValueError(
                f"Metadata 'notes' must contain one of: {valid_phrases}. "
                f"Received: '{v['notes']}'"
            )
        return v


class RBFEEdge(RelativeEdge):
    """
    Specific schema for Relative Binding Free Energy.
    Adds protein-specific requirements.
    """

    protein_paths: List[FilePath] = Field(
        ..., description="Paths to the protein topology and coordinate files"
    )
    protein_ff: ProteinForcefield = ProteinForcefield.AMBER14SB


class RHFEEdge(RelativeEdge):
    """
    Specific schema for Relative Hydration Free Energy.
    Inherits from RelativeEdge, but requires no protein.
    """

    pass


# ==========================================
# Absolute Calculations (AHFE)
# ==========================================


class AHFENode(BaseConfig):
    """
    Specific schema for Absolute Hydration Free Energy.
    Deals with a single ligand (node), so no ligand_b or mapping is required.
    """

    @field_validator("metadata")
    @classmethod
    def validate_metadata_notes(cls, v: Dict[str, Any]) -> Dict[str, Any]:
        """Validates that metadata explicitly mentions decoupling/annihilation."""
        notes = str(v.get("notes", "")).lower()
        valid_phrases = ["annihilation", "decoupling"]

        if not any(phrase in notes for phrase in valid_phrases):
            raise ValueError(
                f"Metadata 'notes' for AHFE must contain one of: {valid_phrases} "
                "to ensure the correct lambda schedule is applied."
            )
        return v
