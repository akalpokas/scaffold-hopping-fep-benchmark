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

NETWORK_FILE=$1
EDGE_ID=$2
PROTOCOL=$3
LEG=$4
DE_STRENGTH=$5
REPLICATE=$6


echo "Running on node: ${SLURM_NODELIST}"
echo "Using network file: ${NETWORK_FILE}"
echo "Starting job for edge: ${EDGE_ID}"
echo "Using run script: ${PROTOCOL}"
echo "Leg: ${LEG}"
echo "DE strength: ${DE_STRENGTH}"
echo "Replicate: ${REPLICATE}"

export NUMEXPR_MAX_THREADS=32


singularity run --nv --bind /scratch/path_to_work_dir/ \
 /scratch/path_to_run_sif_image/alchemate_feat_dmr.sif \
 python alchemate_run_edge.py --network-file "$NETWORK_FILE" --edge-id "$EDGE_ID" --protocol "$PROTOCOL" --leg_name "$LEG"  --de_strength "$DE_STRENGTH" --replicate "$REPLICATE"
