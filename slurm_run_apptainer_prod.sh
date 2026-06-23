#!/bin/bash
#SBATCH -o isa_somd2.%j.slurm.out
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# Exit immediately if a command exits with a non-zero status
set -e

NETWORK=$1
EDGE_ID=$2
PROTOCOL=$3
LEG=$4
BOND_STRENGTH=$5
DE_STRENGTH=$6
REPLICATE=$7
GHOST_MODS=$8

# Conditionally construct the ghost modifications flag
GHOST_FLAG=""
if [ "$GHOST_MODS" = "True" ]; then
    GHOST_FLAG="--ghost_modifications"
fi

echo "Running on node: ${SLURM_NODELIST}"
echo "Using network file: ${NETWORK}"
echo "Starting job for edge: ${EDGE_ID}"
echo "Using run script: ${PROTOCOL}"
echo "Leg: ${LEG}"
echo "Soft Morse Bond strength: ${BOND_STRENGTH}"
echo "Soft Morse DE strength: ${DE_STRENGTH}"
echo "Replicate: ${REPLICATE}"
echo "Ghost Modifications: ${GHOST_MODS}"

export NUMEXPR_MAX_THREADS=32

apptainer run --nv --bind /path/to/benchmark/ \
 /path/to/apptainer/image.sif \
 python edge_runner.py \
 --network "$NETWORK" \
 --edge-id "$EDGE_ID" \
 --protocol "$PROTOCOL" \
 --leg_name "$LEG" \
 --bond_strength "$BOND_STRENGTH" \
 --de_strength "$DE_STRENGTH" \
 --replicate "$REPLICATE" \
 $GHOST_FLAG
