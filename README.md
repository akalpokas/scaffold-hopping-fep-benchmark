# scaffold-hopping-fep-benchmark

> [!CAUTION]
> The code in this repository is under active development and is not yet ready for general use.

# Instructions

Tools - Use them to process a network from start to finish.  
Templates - Modify them once to align to specific HPC or simulation processing needs.

| Script    | Type | Purpose |
| -------- | -------- | ------- |
| `rbfe_pipeline_prep.py` | Tool | Constructs the alchemical inputs & production ready SOMD2 run files from raw inputs |
| `deploy.py` | Tool | Deploys the network and built alchemical inputs for processing |
| `analyse.py` |  Tool | Runs various analysis workflows on the processed network |
| `slurm_run_apptainer_prod.sh` | Template | Sets up the slurm template for individual HPC needs |
| `edge_runner.py` | Template | Runs the alchemical transformation using alchemate workflows and SOMD2 |
| `clean_runs.py` | Tool | Convenience script for removing previously ran simulations with a specific protocol |

## Alchemical input generation
To process the mappings for a given network:
```python
python rbfe_pipeline_prep.py map --network network.json
```
To process the mappings for a given network and a specific edge:
```python
python rbfe_pipeline_prep.py map --network network.json --edge ligA_to_ligB

# For example:
python rbfe_pipeline_prep.py map --network networks/zou_network.json --edge-id chk1_c20_to_c17
```
To run the full parametrisation stage and alchemical setup:
```python
python rbfe_pipeline_prep.py setup --network network.json
```

## Network Deployment
To deploy a single replicate testing protocol run for free and bound legs:
```python
python deploy.py --network networks/zou_network.json --protocol testing --leg both --replicate 1
```

## Analysis

To run a basic analysis workflow on all 1st free leg replicates in a given network:
```python
python analyse.py --network networks/zou_network.json --modules energy_traj --protocol testing --leg_name free --k 125 --de 150 --replicate 1
```

# Data Flow
```mermaid
flowchart TD
    classDef file fill:#f9f9f9,stroke:#333,stroke-width:2px;
    classDef script fill:#e1f5fe,stroke:#0288d1,stroke-width:2px;
    classDef check fill:#fff3e0,stroke:#f57c00,stroke-width:2px;

    subgraph Inputs ["1. Raw Inputs"]
        JSON(network.json):::file
        SDF(Ligand SDFs/MOLs):::file
        PDB(Protein PDBs):::file
    end

    subgraph Prep ["2. Preparation Pipeline"]
        Schema{{Pydantic Schema}}:::check
        PipeMap[rbfe_pipeline.py stage: map]:::script
        PipeSetup[rbfe_pipeline.py stage: setup]:::script
        
        Schema --> PipeMap
        Schema --> PipeSetup
    end

    subgraph LocalOut ["3. Edge Outputs"]
        Vis(mapping_vis.png):::file
        CSV(changed_bonds.csv):::file
        Free(Free Leg: bss):::file
        Bound(Bound Leg: bss):::file
    end

    subgraph Cluster ["4. HPC Deployment"]
        Deploy[deploy.py]:::script
        Slurm((SLURM))
        Bash[slurm_run_apptainer_prod.sh]:::script
        Run[edge_runner.py]:::script
        MD[(Raw MD Data)]
    end

    subgraph Analysis ["5. Analysis"]
        Analyse[analyse.py]:::script
        FinalData[(Processed MD Data)]
    end

    %% Step 1 to 2
    JSON --> Schema
    SDF --> Schema
    PDB --> Schema

    %% Step 2 to 3
    PipeMap -->|Generates| Vis
    PipeSetup -->|Generates| CSV
    PipeSetup -->|Generates| Free
    PipeSetup -->|Generates| Bound

    %% Step 3 to 4
    JSON -.->|Reads Edge IDs| Deploy
    Deploy -->|Loops edges, legs, reps| Slurm
    Slurm -->|Uses the template slurm script| Bash
    Bash --> Run
    
    Free -.->|Loaded by| Run
    Bound -.->|Loaded by| Run
    
    Run -->|Executes SOMD2| MD

    %% Step 4 to 5
    JSON -->|Analysed by| Analyse
    MD -->|Analysed by| Analyse
    Analyse -->| Produces | FinalData
```

# Software versioning

## Setup software
`biosimspace - 2026.1.0.dev0` (conda install)
`sire - 2026.1.0.dev0`        (conda install)

## Run images
- Default: [2026.05.20](https://hub.docker.com/repository/docker/akalpokas/alchemate_rb/tags/2026.05.20)
- REST2 Angle Tempering: 