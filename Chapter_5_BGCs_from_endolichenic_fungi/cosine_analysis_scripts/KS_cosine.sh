#!/bin/bash
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH --partition=bigmem
#SBATCH --time=2-00:00:00
#SBATCH -o min_cosine_conversion_%A.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=peng.hou@vuw.ac.nz

# Writes all split CSV filenames into a list
ls file*.csv >> ELF_ref.list

# uses GNU parallel to pass each CSV to `min-avg_cosine.py` for processing.
cat ELF_ref.list | parallel python3 mini-avg_cosine.py {} ELF_identifier.csv
