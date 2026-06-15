import argparse
import json
from shutil import copy2

from somd2.config import Config
from alchemate.manager import WorkflowManager
from alchemate.context import SimulationContext
from alchemate.steps.preprocessing import OptimizeLambdaProbabilities
from alchemate.steps.base import RunBasicCalculation
from alchemate.logger import setup_logging

# Import schema from the setup rbfe pipeline
from rbfe_pipeline_prep import RBFEEdge

import sire as sr


def get_user_input():
    parser = argparse.ArgumentParser(description="Simulation parameters")
    parser.add_argument(
        "--network",
        type=str,
        required=True,
        help="Path to the network file containing the edges to run",
    )
    parser.add_argument(
        "--edge-id",
        type=str,
        required=True,
        help="Edge ID to run. Should match an edge_id in the network file",
    )
    parser.add_argument(
        "--protocol",
        type=str,
        required=True,
        help="Type of runtime protocol to use. Options: 'testing', 'prod', 'long'",
    )
    parser.add_argument("--de_strength", type=str, required=True, help="DE strength")
    parser.add_argument(
        "--leg_name",
        type=str,
        required=True,
        help="Leg to run. Options: 'free' or 'bound'",
    )
    parser.add_argument(
        "--replicate",
        type=int,
        required=True,
        help="Replicate to run. If not provided, all replicates will be run",
    )
    parser.add_argument(
        "--bond_strength",
        type=str,
        required=False,
        help="Bond strength for soft morse potential. If not provided, a default value of 125 kcal/mol/A^2 will be used",
        default=125,
    )
    parser.add_argument(
        "--ghost_modifications",
        action="store_true",
        help="Whether to apply ghost modifications to the system. If not provided, ghost modifications will not be applied",
    )

    args = parser.parse_args()
    return (
        args.network,
        args.edge_id,
        args.protocol,
        args.de_strength,
        args.leg_name,
        args.replicate,
        args.bond_strength,
        args.ghost_modifications,
    )


def main():
    (
        network,
        edge_id,
        protocol,
        de_strength,
        leg_name,
        replicate,
        bond_strength,
        ghost_modifications,
    ) = get_user_input()

    # Load the edge information from the network file
    with open(network, "r") as f:
        network_data = json.load(f)
    edge_dict = next((item for item in network_data if item["edge_id"] == edge_id), None)

    if not edge_dict:
        raise ValueError(f"Edge {edge_id} not found in {network}")

    # validate schema
    edge_config = RBFEEdge(**edge_dict)

    somd2_config = Config()
    somd2_config.lambda_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    sire_system = sr.stream.load(f"{edge_config.output_dir}/{leg_name}.bss")

    # attempt to determine the lambda schedule based on edge metadata
    metadata = edge_config.metadata

    if (
        metadata
        and "notes" in metadata
        and "bond annihilation" in metadata["notes"].lower()
    ):
        somd2_config.lambda_schedule = "ring_break_morph"
        bond_alchemy = True
    elif metadata and "notes" in metadata and "bond creation" in metadata["notes"].lower():
        somd2_config.lambda_schedule = "ring_make_morph_reverse"
        bond_alchemy = True
    elif metadata and "notes" in metadata and "standard morph" in metadata["notes"].lower():
        somd2_config.lambda_schedule = "standard_morph"
        bond_alchemy = False
    else:
        raise ValueError(f"""Unable to determine lambda schedule from edge metadata.
                        Please ensure that the 'notes' field in the edge metadata contains one of the following:
                        'bond annihilation', 'bond creation', or 'standard morph'.
                        Edge metadata: {metadata}""")

    if bond_alchemy:
        hard_restraints, sire_system = sr.restraints.morse_potential(
            sire_system,
            de="150 kcal mol-1",
            auto_parametrise=True,
            direct_morse_replacement=True,
            name="morse_hard",
        )
        print(hard_restraints)

        soft_restraints, _ = sr.restraints.morse_potential(
            sire_system,
            atoms0=hard_restraints[0].atom0(),
            atoms1=hard_restraints[0].atom1(),
            r0=hard_restraints[0].r0(),
            k=f"{bond_strength} kcal mol-1 A-2",
            auto_parametrise=False,
            de=f"{de_strength} kcal mol-1",
            name="morse_soft",
        )
        print(soft_restraints)
        somd2_config.restraints = [hard_restraints, soft_restraints]
        
    print(f"Sire system type:{type(sire_system)}")

    somd2_config.timestep = "4fs"
    if protocol == "testing":
        equib_time = 100
        prod_time = 1000
        frame_freq = 250
        checkpoint_freq = 500

        somd2_config.save_crash_report = True
        somd2_config.save_energy_components = True
    elif protocol == "prod":
        equib_time = 500
        prod_time = 10000
        frame_freq = 250
        checkpoint_freq = 1000
    elif protocol == "prod_2fs":
        equib_time = 500
        prod_time = 10000
        frame_freq = 250
        checkpoint_freq = 1000
        somd2_config.timestep = "2fs"
    elif protocol == "long":
        equib_time = 1000
        prod_time = 25000
        frame_freq = 250
        checkpoint_freq = 1000
        
    else:
        raise ValueError(
            f"Invalid protocol: {protocol}. Options are: 'testing', 'prod', 'prod_2fs', 'long'"
        )

    if ghost_modifications:
        somd2_config.ghost_modifications = True
        mods_prefix = "ghost_mods"
    else:
        somd2_config.ghost_modifications = False
        mods_prefix = ""

    somd2_config.output_directory = f"{edge_config.output_dir}/{leg_name}_k_{int(bond_strength)}_{mods_prefix}_de_{int(de_strength)}_{protocol}_protocol_repl_{replicate}"
    somd2_config.equilibration_time = f"{equib_time}ps"
    somd2_config.runtime = f"{prod_time}ps"
    somd2_config.frame_frequency = f"{frame_freq}ps"
    somd2_config.checkpoint_frequency = f"{checkpoint_freq}ps"

    
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

    context = SimulationContext(
        system=f"{edge_config.output_dir}/{leg_name}.bss", somd2_config=somd2_config
    )

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
    final_context = manager.execute()

if __name__ == "__main__":
    main()