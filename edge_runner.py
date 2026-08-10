import argparse
import json
from pathlib import Path
import re
from pipeline._utils import validate_protocol

from somd2.config import Config
from alchemate.manager import WorkflowManager
from alchemate.context import SimulationContext
from alchemate.steps.preprocessing import OptimizeLambdaProbabilities
from alchemate.steps.base import RunBasicCalculation
from alchemate.logger import setup_logging

import sire as sr
from pipeline.schemas import RBFEEdge, RHFEEdge, AHFENode


def get_user_input():
    parser = argparse.ArgumentParser(
        description="Unified Simulation Runner (RBFE, RHFE, AHFE)"
    )
    parser.add_argument(
        "--network", type=str, required=True, help="Path to the JSON network file"
    )
    parser.add_argument(
        "--edge-id", type=str, required=True, help="Edge or Node ID to run"
    )
    parser.add_argument(
        "--protocol",
        type=validate_protocol,
        required=True,
        help="Protocol string consisting of a base and optional modifiers separated by underscores.",
    )
    parser.add_argument(
        "--leg_name",
        type=str,
        required=True,
        choices=["bound", "free", "solvent", "vacuum"],
    )
    parser.add_argument("--replicate", type=int, required=True, help="Replicate number")

    # Optional parameters (used for ring-breaking RBFE/RHFE)
    parser.add_argument(
        "--de_strength", type=str, default="50", help="DE strength for Morse potential"
    )
    parser.add_argument(
        "--bond_strength",
        type=str,
        default="125",
        help="Bond strength for soft morse potential",
    )
    parser.add_argument(
        "--ghost_modifications", action="store_true", help="Apply ghost modifications"
    )

    return parser.parse_args()


def parse_config(edge_dict: dict):
    """Dynamically parses the dictionary into the correct Pydantic schema."""
    if "protein_paths" in edge_dict:
        return RBFEEdge(**edge_dict)
    elif "ligand_b_paths" in edge_dict:
        return RHFEEdge(**edge_dict)
    else:
        return AHFENode(**edge_dict)


def main():
    args = get_user_input()

    # Load the edge information
    with open(args.network, "r") as f:
        network_data = json.load(f)

    edge_dict = next(
        (item for item in network_data if item["edge_id"] == args.edge_id), None
    )
    if not edge_dict:
        raise ValueError(f"ID '{args.edge_id}' not found in {args.network}")

    # Parse into the correct schema (RBFE, RHFE, or AHFE)
    config = parse_config(edge_dict)

    # Base Configuration
    somd2_config = Config()
    somd2_config.lambda_values = [
        0.0,
        0.1,
        0.2,
        0.3,
        0.35,
        0.4,
        0.425,
        0.45,
        0.5,
        0.55,
        0.575,
        0.6,
        0.625,
        0.637,
        0.65,
        0.675,
        0.7,
        0.8,
        0.9,
        1.0,
    ]

    # Load System
    system_path = config.output_dir / f"{args.leg_name}.bss"
    if not system_path.exists():
        raise FileNotFoundError(
            f"System file not found: {system_path}. Did you run setup for this leg?"
        )

    sire_system = sr.stream.load(str(system_path))

    # ==========================================
    # Lambda Schedule & Restraints Logic
    # ==========================================
    metadata_notes = str(config.metadata.get("notes", "")).lower()
    bond_alchemy = False

    if isinstance(config, (RBFEEdge, RHFEEdge)):
        if "bond annihilation" in metadata_notes:
            somd2_config.lambda_schedule = "ring_break_morph"
            bond_alchemy = True
        elif "bond creation" in metadata_notes:
            somd2_config.lambda_schedule = "reverse_ring_break_morph"
            bond_alchemy = True
        elif "standard morph" in metadata_notes:
            somd2_config.lambda_schedule = "standard_morph"
        else:
            raise ValueError(f"Unrecognized Relative morph schedule: {metadata_notes}")

    elif isinstance(config, AHFENode):
        if "decouple" in metadata_notes or "decoupling" in metadata_notes:
            somd2_config.lambda_schedule = "decouple"
        else:
            raise ValueError(f"Unrecognized Absolute schedule: {metadata_notes}")

    # Apply Morse Potential only if ring breaking is detected
    if bond_alchemy:
        # Validate that the user-specified transformation matches the detected bond changes
        # Note that this needs to be done before DMR as it will delete the detectable bond
        mol = sire_system.molecules("molecule property is_perturbable")
        ref_mol = sr.morph.link_to_reference(mol)
        perturbable_mol = ref_mol[0]
        pert = perturbable_mol.perturbation()
        pert_omm = pert.to_openmm()
        changed_bonds_df = pert_omm.changed_bonds(to_pandas=True)
        n_bonds_created = (changed_bonds_df["k0"] == 0).sum()
        n_bonds_annihilated = (changed_bonds_df["k1"] == 0).sum()

        if n_bonds_created == 1:
            assert "bond creation" in metadata_notes
        elif n_bonds_annihilated == 1:
            assert "bond annihilation" in metadata_notes
        else:
            raise ValueError(
                "The user specified transformation does not match the detected bond changes."
            )
        hard_restraints, sire_system = sr.restraints.morse_potential(
            sire_system,
            de="150 kcal mol-1",
            auto_parametrise=True,
            direct_morse_replacement=True,
            name="morse_hard",
        )
        soft_restraints, _ = sr.restraints.morse_potential(
            sire_system,
            atoms0=hard_restraints[0].atom0(),
            atoms1=hard_restraints[0].atom1(),
            r0=hard_restraints[0].r0(),
            k=f"{args.bond_strength} kcal mol-1 A-2",
            auto_parametrise=False,
            de=f"{args.de_strength} kcal mol-1",
            name="morse_soft",
        )
        somd2_config.restraints = [hard_restraints, soft_restraints]

        # Now figure out which angles are being annihilated and use the atom numbers for tempering
        changed_angles_df = pert_omm.changed_angles(to_pandas=True)

        # Filter the DataFrame for rows where k0 == 0 OR k1 == 0 (for reverse and forward transformations)
        filtered_df = changed_angles_df[
            (changed_angles_df["k0"] == 0) | (changed_angles_df["k1"] == 0)
        ]

        # Extract the atom numbers
        angle_atom_numbers = (
            filtered_df["angle"]
            .apply(lambda x: [int(n) for n in re.findall(r":(\d+)", str(x))])
            .tolist()
        )

        # Flatten the list and get unique values
        unique_annihilated_atom_num = sorted(
            set(atom for sublist in angle_atom_numbers for atom in sublist)
        )

    # ==========================================
    # MD Protocol Settings
    # ==========================================
    DEFAULT_PARAMS = {
        "equilibration_timestep": "2fs",
        "energy_frequency": "1ps",
        "cutoff": "10A",
        "cutoff_type": "PME",
        "equilibration_constraints": True,
        "num_energy_neighbours": 5,
        "h_mass_factor": 3,
        "rest2_scale": 1,
        "replica_exchange": True,
        "log_level": "debug",
        "save_xml": True,
        "constraint": "bonds",
        "timeout": "30s",
        "shift_delta": "1.5A",
        "shift_coulomb": "1A",
        "save_energy_components": False,
    }

    # Define base setups
    BASE_PROTOCOLS = {
        "testing": {
            "equilibration_time": "100ps",
            "runtime": "1000ps",
            "frame_frequency": "250ps",
            "checkpoint_frequency": "500ps",
            "save_crash_report": True,
            "save_energy_components": True,
        },
        "prod": {
            "equilibration_time": "500ps",
            "runtime": "10000ps",
            "frame_frequency": "250ps",
            "checkpoint_frequency": "1000ps",
        },
        "tucker": {
            "equilibration_time": "250ps",
            "runtime": "10000ps",
            "frame_frequency": "100ps",
            "checkpoint_frequency": "500ps",
            "timestep": "2fs",
            "equilibration_timestep": "1fs",
            "energy_frequency": "2ps",
            "cutoff": "14A",
            "shift_delta": "2.25A",
        },
    }

    # Define modular overrides
    MODIFIERS = {
        "long": {
            "equilibration_time": "1000ps",
            "runtime": "25000ps",
        },
        "short": {
            "equilibration_time": "1000ps",
            "runtime": "5000ps",
        },
        "rest2": {
            "rest2_scale": 2,
        },
        "targetAngleRest2": {
            "rest2_scale": 2,
        },
        "2fs": {
            "timestep": "2fs",
        },
    }

    # Parse the user's requested protocol (e.g., "tucker_long_rest2")
    # We assume the first word is the base, and anything after an underscore is a modifier.
    protocol_parts = args.protocol.split("_")
    base_name = protocol_parts[0]
    requested_modifiers = protocol_parts[1:]

    if base_name not in BASE_PROTOCOLS:
        raise ValueError(
            f"Unknown base protocol: '{base_name}'. Available bases: {list(BASE_PROTOCOLS.keys())}"
        )

    # Build the final configuration sequentially
    # Start with defaults, overwrite with base protocol, then overwrite with modifiers
    final_config = DEFAULT_PARAMS.copy()
    final_config.update(BASE_PROTOCOLS[base_name])

    for mod in requested_modifiers:
        if mod not in MODIFIERS:
            raise ValueError(
                f"Unknown modifier: '{mod}'. Available modifiers: {list(MODIFIERS.keys())}"
            )
        final_config.update(MODIFIERS[mod])

    # Apply all parameters to the somd2_config object
    for param_name, param_value in final_config.items():
        setattr(somd2_config, param_name, param_value)

    # For rest2, add rest2 angle atoms if they exist
    if "targetAngleRest2" in requested_modifiers:
        if unique_annihilated_atom_num:
            somd2_config.rest2_selection = (
                "property is_perturbable and atomnum "
                + ", ".join(map(str, unique_annihilated_atom_num))
            )
        else:
            raise ValueError(
                "No unique annihilated atom numbers found for rest2 angle selection."
            )

    # ==========================================
    # Directory Naming & Ghost Mods
    # ==========================================
    mods_prefix = "ghost_mods_" if args.ghost_modifications else ""
    somd2_config.ghost_modifications = args.ghost_modifications

    # Clean up output directory name depending on whether bond alchemy was used
    if bond_alchemy:
        out_name = f"{args.leg_name}_k_{int(args.bond_strength)}_{mods_prefix}de_{int(args.de_strength)}_{args.protocol}_repl_{args.replicate}"
    else:
        out_name = f"{args.leg_name}_{mods_prefix}{args.protocol}_repl_{args.replicate}"

    somd2_config.output_directory = str(config.output_dir / out_name)

    # create the output directory if it doesn't exist
    Path(somd2_config.output_directory).mkdir(parents=True, exist_ok=True)

    # ==========================================
    # Execution
    # ==========================================

    context = SimulationContext(system=str(system_path), somd2_config=somd2_config)

    setup_logging(log_path=f"{context.somd2_config.output_directory}/alchemate.log")

    simulation_workflow = [
        OptimizeLambdaProbabilities(
            optimization_attempts=10,
            optimization_target="overlap_matrix",
            optimization_threshold=0.1,
            optimization_runtime="500ps",
            vacuum_optimization=False,
        ),
        RunBasicCalculation(calculation_runtime=f"{somd2_config.runtime}"),
    ]

    manager = WorkflowManager(context=context, workflow_steps=simulation_workflow)
    manager.execute()


if __name__ == "__main__":
    main()
