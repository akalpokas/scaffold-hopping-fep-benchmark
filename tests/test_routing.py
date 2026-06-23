import json
import pytest
from pathlib import Path

from setup_hfe import load_and_route_network
from pipeline.schemas import RHFEEdge, AHFENode

TEST_DIR = Path(__file__).resolve().parent
INPUTS_DIR = TEST_DIR / "inputs"

LIG_A = str(INPUTS_DIR / "dummy_A.sdf")
LIG_B = str(INPUTS_DIR / "dummy_B.sdf")


def test_hfe_router(tmp_path):
    """
    Test that the JSON router correctly identifies RHFE vs AHFE
    and drops invalid schemas without crashing the pipeline.
    """
    mixed_network = [
        # 1. Valid RHFE Edge
        {
            "edge_id": "morph_1",
            "metadata": {"notes": "bond annihilation"},
            "ligand_ff": "gaff2",
            "output_dir": "out/rhfe_1",
            "mapping": {"0": 1},
            "ligand_a_paths": [LIG_A],
            "ligand_b_paths": [LIG_B],
        },
        # 2. Valid AHFE Node
        {
            "edge_id": "node_1",
            "metadata": {"notes": "decoupling"},
            "ligand_ff": "gaff2",
            "output_dir": "out/ahfe_1",
            "ligand_a_paths": [LIG_A],
        },
        # 3. Invalid Entry (No paths)
        {
            "edge_id": "broken_edge",
            "metadata": {"notes": "something random"},
            "ligand_ff": "gaff2",
            "output_dir": "out/broken",
        },
    ]

    # Write the network dictionary to a temp JSON file
    json_path = tmp_path / "test_network.json"
    with open(json_path, "w") as f:
        json.dump(mixed_network, f)

    parsed_configs = load_and_route_network(json_path)

    assert (
        len(parsed_configs) == 2
    ), "Router should have kept exactly 2 valid configs and dropped 1"
    assert isinstance(parsed_configs[0], RHFEEdge)
    assert parsed_configs[0].id == "morph_1"
    assert isinstance(parsed_configs[1], AHFENode)
    assert parsed_configs[1].id == "node_1"
