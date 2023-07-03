#!/bin/bash
#SBATCH --array=1-100%5
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --partition=parallel
#SBATCH --time=15:00
#SBATCH --job-name=cluster_search
#SBATCH -o cluster_search_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

# pHMM directories (Uncomment depending on which model set you'd like to apply)
phmms="curated_phmms"
#phmms="swissprot_phmms"

time python3 cluster_search_v3.py \
-g ${SLURM_ARRAY_TASK_ID} \
-p $phmms \
-o ${SLURM_ARRAY_TASK_ID}/output