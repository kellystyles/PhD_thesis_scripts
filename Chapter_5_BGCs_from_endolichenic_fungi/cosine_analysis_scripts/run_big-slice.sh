#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem=128G
#SBATCH --partition=bigmem
#SBATCH --time=3-00:00:00
#SBATCH -o bigslice_%A.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=username@vuw.ac.nz

"""
This script will run BiG-SLiCE on an `input_folder` containing:
 - a `datasets.tsv` metadata file
 - a directory (as named in the metadata) containing:
    - a LIST file of the BGC GBK files
    - the BGC GBK files
"""

bigslice --query input_folder --n_ranks 1 -t ${SLURM_CPUS_PER_TASK} --run_id 1 --query_name ELF_BGC output_folder