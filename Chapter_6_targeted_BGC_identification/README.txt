# Chapter 6 - Targeted identification of IDT BGCs in low-coverage ELF genomes

Dependencies:
- a Python 3.8+ environment with:
	- Singularity v3.8+ (optional)
	- pandas
	- seaborn
	- upsetplot
	- pyhmmer
- the above environment could be prepared using the commands:
```
conda create -n chapter_6 python3
conda install -n chapter_6 -c bioconda pyhmmer pandas seaborn singularity upsetplot
source activate chapter_6
```

This chapter discusses the various ways to identify BGCs and the challenges of doing so with highly fragmented genome assemblies.
It discusses the use of profile Hidden Markov Models (pHMMs) as one method to quickly identify a similar sequence from a large volume of sequence data.
The chapter then discusse the development of the `cluster_search` script which utilises pHMMs to search for input pHMMs and extract a GenBank file, from annotated GenBank assemblies.

A key experiment was the determination of which genome assemblies possessed an indole diterpene BGC. BLASTp was used as the ground truth for this and the efficacy of various pHMM sets (prepared from curated proteins or iteratively) were compared to this using `cluster_search`.

The overall workflow of the results for this chapter are as follows:
1. Run Jupyter Notebook `iterative_pyhmmer_build.ipynb` to obtain iterative SwissProt models and plots.
2. Download all fungal genomes with annotated assemblies from GenBank using the script `genome_getter.sh`.
3. BLASTp each genome against a BLAST DB of curated IDT protein sequences with `genome_checker.sh` script.
 - multiFASTAs to generate each IDT BLAST DB are the same as the ones used to generate the curated pHMMs.
4. Prepare test databases using the `prepare_test_dbs.py` script.
 - `python3 prepare_synthetic_dbs.py 100`
5. Each set of pHMMs was applied using the `cluster_search.sh` script.
 - run this script once for each pHMM set
6. Plots were generated using the Jupyter Notebook `identify_IDT_genomes.ipynb`.

