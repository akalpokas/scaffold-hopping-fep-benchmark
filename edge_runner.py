import argparse
import json

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
        type=str,
        required=True,
        choices=["testing", "prod", "prod_2fs", "long"],
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
    somd2_config.lambda_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

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

    # ==========================================
    # MD Protocol Settings
    # ==========================================
    somd2_config.timestep = "4fs"

    if args.protocol == "testing":
        equib_time, prod_time, frame_freq, checkpoint_freq = 100, 1000, 250, 500
        somd2_config.save_crash_report = True
        somd2_config.save_energy_components = True
    elif args.protocol == "prod":
        equib_time, prod_time, frame_freq, checkpoint_freq = 500, 10000, 250, 1000
    elif args.protocol == "prod_2fs":
        equib_time, prod_time, frame_freq, checkpoint_freq = 500, 10000, 250, 1000
        somd2_config.timestep = "2fs"
    elif args.protocol == "long":
        equib_time, prod_time, frame_freq, checkpoint_freq = 1000, 25000, 250, 1000

    somd2_config.equilibration_time = f"{equib_time}ps"
    somd2_config.runtime = f"{prod_time}ps"
    somd2_config.frame_frequency = f"{frame_freq}ps"
    somd2_config.checkpoint_frequency = f"{checkpoint_freq}ps"

    # Base parameters
    somd2_config.equilibration_timestep = "2fs"
    somd2_config.energy_frequency = "1ps"
    somd2_config.cutoff = "10A"
    somd2_config.cutoff_type = "PME"
    somd2_config.equilibration_constraints = True
    somd2_config.num_energy_neighbours = 5
    somd2_config.h_mass_factor = 3
    somd2_config.rest2_scale = 1
    somd2_config.replica_exchange = True
    somd2_config.log_level = "debug"
    somd2_config.save_xml = True
    somd2_config.constraint = "bonds"
    somd2_config.timeout = "30s"
    somd2_config.shift_delta = "1.5A"
    somd2_config.shift_coulomb = "1A"

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
        RunBasicCalculation(calculation_runtime=f"{prod_time}ps"),
    ]

    manager = WorkflowManager(context=context, workflow_steps=simulation_workflow)
    manager.execute()


if __name__ == "__main__":
    main()
