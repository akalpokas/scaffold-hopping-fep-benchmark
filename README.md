# scaffold-hopping-fep-benchmark

# Instructions

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


# Run scripts

# Software versioning

## Setup software
`biosimspace - 2026.1.0.dev0` (conda install)
`sire - 2026.1.0.dev0`        (conda install)

## Run images
- Default: [2026.05.20](https://hub.docker.com/repository/docker/akalpokas/alchemate_rb/tags/2026.05.20)
- REST2 Angle Tempering: 