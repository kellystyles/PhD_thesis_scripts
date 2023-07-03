#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --partition=parallel
#SBATCH --time=12:00:00
#SBATCH --job-name=genome_checks
#SBATCH -o genome_checks_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

# This script utilises the Singularity container `fungiflow_v3.sif`, which contains
# an installation of BLAST. A local installation of BLAST can be used by omitting 
# all lines `singularity exec ${singularity} \`.

# If using the Singularity container, edit the path variable below

(echo "Script begun " && date)
module purge
module load singularity/3.7.3

singularity="/path/to/fungiflow_v3.sif"
threshold="-20"
# prepared IDT blastdb consisting of B,C,M protein seqs only
singularity exec ${singularity} \
makeblastdb -in idtB_curated.fasta -dbtype prot -input_type fasta -title idtB -out dbs/idtB
singularity exec ${singularity} \
makeblastdb -in idtC_curated.fasta -dbtype prot -input_type fasta -title idtC -out dbs/idtC
singularity exec ${singularity} \
makeblastdb -in idtM_curated.fasta -dbtype prot -input_type fasta -title idtM -out dbs/idtM

# BLASTp each *.faa genome that has an associated GBK file (to match expected behaviour of pipeline)
# If the *.faa genome file has a B, C, and M hit then consider it a 'non-negative genome' and move to folder, otherwise move the 'gbff' file to the 'negative' directory
mkdir -p non-negative && mkdir -p negative && mkdir -p blast_out

NUM_CORES=4

# Define a function to process each file in parallel
process_file() {
    file=$1
    echo $file
    name=$(echo $file | rev | cut -d"/" -f2 | rev)
    echo $name
    if [ ! -f blast_out/"$name"_B.txt ]; then
        singularity exec ${singularity}  \
            blastp -db dbs/idtB -query "$file" -outfmt 6 -evalue 1e-20 > blast_out/"$name"_B.txt
    else
        echo "Skipping B BLASTp for $name as output exists"
    fi
    if [ ! -f blast_out/"$name"_C.txt ]; then
        singularity exec ${singularity}  \
            blastp -db dbs/idtB -query "$file" -outfmt 6 -evalue 1e-20 > blast_out/"$name"_C.txt
    else
        echo "Skipping C BLASTp for $name as output exists"
    fi

    if [ ! -f blast_out/"$name"_M.txt ]; then
        singularity exec ${singularity}  \
            blastp -db dbs/idtB -query "$file" -outfmt 6 -evalue 1e-20 > blast_out/"$name"_M.txt
    else
        echo "Skipping $name as BLASTp output exists"
    fi

    if [ ! -f blast_out/"$name"_M.txt ]; then
        singularity exec ${singularity}  \
            blastp -db dbs/idtG -query "$file" -outfmt 6 -evalue 1e-20 > blast_out/"$name"_G.txt
    else
        echo "Skipping $name as BLASTp output exists"
    fi

	B=$(head -n1 blast_out/"$name"_B.txt | cut -f11 | cut -d"e" -f2 ) 		
	C=$( head -n1 blast_out/"$name"_C.txt | cut -f11 | cut -d"e" -f2 ) 	
	M=$( head -n1 blast_out/"$name"_M.txt | cut -f11 | cut -d"e" -f2 ) 	
	G=$( head -n1 blast_out/"$name"_G.txt | cut -f11 | cut -d"e" -f2 ) 

    ## Here we assume that an IDT BGC contains hits to IdtBCM, whilst for my thesis 
    ## I only used IdtCM as I observed poor annotation of IdtB in several of my 
    ## assemblies. This would require ommision of the first 
    ## `(( "$B" < ${working_threshold} ))` block in the code below.
	if (( "$B" < ${working_threshold} )) && (( "$C" < ${working_threshold} )) && (( "$M" < ${working_threshold} )); then
	if (( "$B" < ${threshold} )) && (( "$C" < ${threshold} )); then
    	echo "$name contains BLASTp hits to B, C, and M proteins. Moving to 'positive' directory"
		mv genomes/"$name" positive/
        designation="positive"
	else
		echo "$name contains NO significant BLASTp hits to B, C, and M proteins. Moving to 'negative' directory"
		mv genomes/"$name" negative/
        designation="negative"
	fi
    echo "${name}\t${B}\t${C}\t${designation}" >> results.txt
}

# Export the function to make it available to GNU Parallel
export -f process_file

# Use GNU Parallel to process the files in parallel
find . -maxdepth 3 -wholename "*/*/*.faa" | parallel -j $NUM_CORES process_file {}
(echo "Script finished on " && date)