#!/bin/bash
#SBATCH --cpus-per-task=30
#SBATCH --mem=128G
#SBATCH --partition=bigmem
#SBATCH --time=1-00:00:00
#SBATCH -o parse_cosine_%A.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=username@vuw.ac.nz

python parse_cosine_1.py ELF.db
python parse_cosine_2.py


