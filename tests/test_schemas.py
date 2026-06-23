import pytest
from pathlib import Path
from pydantic import ValidationError

from pipeline.schemas import (
    RBFEEdge,
    RHFEEdge,
    AHFENode,
    LigandForcefield,
    ProteinForcefield,
)

# ==========================================
# Resolve Paths Dynamically
# ==========================================
# This gets the absolute path to the 'tests/inputs' directory
TEST_DIR = Path(__file__).resolve().parent
INPUTS_DIR = TEST_DIR / "inputs"

LIG_A = str(INPUTS_DIR / "dummy_A.sdf")
LIG_B = str(INPUTS_DIR / "dummy_B.sdf")
PROT = str(INPUTS_DIR / "dummy_protein.pdb")

# ==========================================
# Reusable Mock Data
# ==========================================


@pytest.fixture
def base_valid_relative_data():
    return {
        "edge_id": "ligA_to_ligB",
        "metadata": {"notes": "standard morph"},
        "ligand_ff": "gaff2",
        "output_dir": "outputs/test_edge",
        "mapping": {"0": 1, "1": 2},
        "ligand_a_paths": [LIG_A],
        "ligand_b_paths": [LIG_B],
    }


@pytest.fixture
def base_valid_ahfe_data():
    return {
        "edge_id": "ligA_node",
        "metadata": {"notes": "standard annihilation"},
        "ligand_ff": "openff",
        "output_dir": "outputs/test_node",
        "ligand_a_paths": [LIG_A],
    }


# ==========================================
# Schema Tests
# ==========================================


def test_valid_rhfe_edge(base_valid_relative_data):
    edge = RHFEEdge(**base_valid_relative_data)
    assert edge.id == "ligA_to_ligB"
    assert edge.ligand_ff == LigandForcefield.GAFF2
    assert len(edge.mapping) == 2


def test_invalid_relative_metadata(base_valid_relative_data):
    bad_data = base_valid_relative_data.copy()
    bad_data["metadata"] = {"notes": "I am just morphing stuff"}

    with pytest.raises(ValidationError) as exc_info:
        RHFEEdge(**bad_data)
    assert "must contain one of" in str(exc_info.value)


def test_missing_mapping_in_relative(base_valid_relative_data):
    bad_data = base_valid_relative_data.copy()
    del bad_data["mapping"]

    with pytest.raises(ValidationError) as exc_info:
        RHFEEdge(**bad_data)
    assert "mapping" in str(exc_info.value)


def test_valid_ahfe_node(base_valid_ahfe_data):
    node = AHFENode(**base_valid_ahfe_data)
    assert node.id == "ligA_node"
    assert not hasattr(node, "mapping")


def test_invalid_ahfe_metadata(base_valid_ahfe_data):
    bad_data = base_valid_ahfe_data.copy()
    bad_data["metadata"] = {"notes": "standard morph"}

    with pytest.raises(ValidationError) as exc_info:
        AHFENode(**bad_data)
    assert "annihilation" in str(exc_info.value)


def test_valid_rbfe_edge(base_valid_relative_data):
    rbfe_data = base_valid_relative_data.copy()
    rbfe_data["protein_paths"] = [PROT]
    rbfe_data["protein_ff"] = "amber14"

    edge = RBFEEdge(**rbfe_data)
    assert edge.protein_ff == ProteinForcefield.AMBER14SB
