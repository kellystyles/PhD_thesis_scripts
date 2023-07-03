#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=parallel
#SBATCH --time=4:00:00
#SBATCH -o FGD_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kellystyles@vuw.ac.nz

path=$(pwd)
n=10000

# downloads latest genome assembly list
# Change this address to the taxa level you would like to download genomes from
# browse taxa at https://ftp.ncbi.nlm.nih.gov/genomes/genbank
if [ -f assembly_summary.txt ]; then
	(date && echo "Genomes list already exists!")
else
	(date && echo "Getting latest genomes list!")
	wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
	# counts the number of genomes for each species and sorts from highest to lowest
	awk -F '\t' '{print $8}' assembly_summary.txt | uniq -c | sort -g -r > info.txt
	# removes non-representative genomes
	# remore this command if you want all genomes for all strains, not just representative genomes
	(date && echo "Culling the fluff...")
	grep "representative genome" assembly_summary.txt > assembly_summary_trimmed.txt  
fi

# download random selection of 'n' genomes
# these will be BLAST checked for IDT proteins
shuf -n ${n} assembly_summary_trimmed.txt > random_lines.txt

# cuts the 20th column (with the ftp addresses for download) and removes the leading "ftp://"
cut -d$'\t' -f20 random_lines.txt | cut -c 9- >  ftp.txt

# for loop to download genome data, uncompress, then delete the compressed file
(date && echo "Preparing to download genomes...")
mkdir -p ${path}/genomes/
cd ${path}/genomes/

NUM_CORES=4

# Define a function to process each file in parallel
process_file() {
    file=$1
	echo $file
    name=$(echo $file | cut -d"." -f2 | cut -d"_" -f1 | cut -d"/" -f2)
	echo $name
    if [ ! -f "$name".bam ]; then
        singularity exec /nfs/scratch/styleske/singularity_files/fungiflow_v3.sif \
            bwa index $file
        singularity exec /nfs/scratch/styleske/singularity_files/fungiflow_v3.sif \
            bwa mem $file ../../raw/"$name"_1.fq.gz ../../raw/"$name"_2.fq.gz > "$name".sam
        singularity exec /nfs/scratch/styleske/singularity_files/fungiflow_v3.sif \
            samtools view -b "$name".sam > "$name".bam
    else
        echo "Skipping $name as BAM exists"
    fi
}

# Export the function to make it available to GNU Parallel
export -f process_file

# Use GNU Parallel to process the files in parallel
find . -maxdepth 1 -name "*.fasta" | parallel -j $NUM_CORES process_file {}