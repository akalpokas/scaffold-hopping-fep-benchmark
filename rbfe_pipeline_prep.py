import json
import argparse
import traceback
from unittest.mock import patch
from pathlib import Path
from enum import Enum
from typing import Any, Dict, Optional
from pydantic import BaseModel, FilePath, Field
import BioSimSpace as BSS

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
        description="Unique identifier for this transformation, e.g., 'lig1_to_lig2'"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict, 
        description="Optional tracking info like experimental ddG, Lomap scores, etc."
    )

    # FilePath ensures the files actually exist before the pipeline does anything
    ligand_a_path: FilePath
    ligand_b_path: FilePath
    protein_path: FilePath
    
    mapping: Dict[int, int] = Field(
        ...,
        description="Atom index mapping from Ligand A to Ligand B"
    )
    
    ligand_ff: LigandForcefield = LigandForcefield.OPENFF
    protein_ff: ProteinForcefield = ProteinForcefield.AMBER14SB
    
    solvent_padding_ang: float = Field(15, ge=0.0, description="Padding for solvent box in Angstroms")
    ionic_strength_molar: float = Field(0.15, ge=0.0 , description="Ionic strength in molar for solvation")
    
    output_dir: Path = Field(..., description="Directory to save the outputs")

    def create_output_dir(self):
        """Utility to ensure the working directory exists."""
        self.output_dir.mkdir(parents=True, exist_ok=True)


# ==========================================
# Core Processing Functions
# ==========================================

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
        # 5. Clean up by restoring the original state
        BSS._is_notebook = original_state

def setup_alchemical_system(config: RBFEEdge):
    # Example placeholder logic
    # out_path = config.output_dir / "system.xml"
    pass


# ==========================================
# Pipeline Orchestration
# ==========================================

def load_network(json_path: Path, specific_edge_id: Optional[str] = None):
    with open(json_path, 'r') as f:
        edges_data = json.load(f)
    
    if specific_edge_id:
        edges_data = [e for e in edges_data if e.get('edge_id') == specific_edge_id]
        
    return edges_data

def run_mapping_stage(network_path: Path, specific_edge_id: str = None):
    edges_data = load_network(network_path, specific_edge_id)
    print(f"Generating mapping visualizations for {len(edges_data)} edges in {network_path.name}...")

    for edge_dict in edges_data:
        edge_name = Path(edge_dict.get('output_dir', 'unknown')).name
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
        edge_name = Path(edge_dict.get('output_dir', 'unknown')).name
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
    parser.add_argument("stage", choices=["map", "setup"], 
                        help="Which stage of the pipeline to run")
    
    # Required network file
    parser.add_argument("--network", required=True, type=Path,
                        help="Path to the JSON network file to process")
    
    # Optional edge filter
    parser.add_argument("--edge-id", default=None,
                        help="Optional: Run a specific edge by its unique ID (e.g., 'lig1_to_lig2')")
    
    args = parser.parse_args()
    
    if args.stage == "map":
        run_mapping_stage(args.network, args.edge_id)
    elif args.stage == "setup":
        run_setup_stage(args.network, args.edge_id)