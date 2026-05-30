# scaffold-hopping-fep-benchmark

# Instructions

| Script    | Purpose |
| -------- | ------- |
| `rbfe_pipeline_prep.py` | Constructs the alchemical inputs & ready SOMD2 run files from raw inputs |
| `deploy.py` | Deploys the network and built alchemical inputs to HPC processing |
| `slurm_run_singularity_prod.sh` | Sets up the slurm template for individual HPC needs |
| `alchemate_run_edge` | Runs the alchemical transformation using alchemate workflows and SOMD2 |

To process the mappings for a given network
```python
python rbfe_pipeline_prep.py map --network network.json
```
To process the mappings for a given network and a specific edge
```python
python rbfe_pipeline_prep.py map --network network.json --edge ligA_to_ligB

# For example:
python rbfe_pipeline_prep.py map --network networks/zou_network.json --edge-id chk1_c20_to_c17
```
To run the full parametrisation stage and alchemical setup
```python
python rbfe_pipeline_prep.py setup --network network.json
```


# Data Flow
```mermaid
flowchart TD
    classDef file fill:#f9f9f9,stroke:#333,stroke-width:1px;
    classDef script fill:#e1f5fe,stroke:#0288d1,stroke-width:2px;
    classDef check fill:#fff3e0,stroke:#f57c00,stroke-width:2px;

    subgraph Inputs ["1. Raw Inputs"]
        JSON(network.json):::file
        SDF(Ligand SDFs):::file
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
        Bash[slurm_run_singularity_prod.sh]:::script
        Run[run_edge.py]:::script
        MD[(Outputs)]
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
    
    Run -->|Executes SOM2| MD
```

# Software versioning

## Setup software
`biosimspace - 2026.1.0.dev0` (conda install)
`sire - 2026.1.0.dev0`        (conda install)

## Run images
- Default: [2026.05.20](https://hub.docker.com/repository/docker/akalpokas/alchemate_rb/tags/2026.05.20)
- REST2 Angle Tempering: 