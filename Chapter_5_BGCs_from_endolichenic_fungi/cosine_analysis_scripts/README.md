# Cosine analysis

**All the scripts and Notebooks discussed here were written by Dr Peng Hou and were modified for clarity by myself.**

One way to assess BGC novelty is to compare newly identified BGCs to known BGCs, both characterized and uncharacterized. 

The ELF BGCs were input into BiG-SLiCE and compared to the 1.2 million BGCs in the BiG-FAM BGC database. The BiG-SLiCE Euclidean distances were converted to cosine vectors as performed in the Paoli et a., 2022[^1] that examined the novelty of BGCs from a global ocean microbiome. Various plots were prepared using a Jupter Notebook (`cosine_plots.ipynb` in the parent directory).

The steps to obtain to obtain the cosine novelty scores, the scripts were run in the following order:
1. Run BiG-SLiCE using the `run_big-slice.sh` SLURM script. This took ~2.5 days to run due to the large number of comparisons required. This will output a large matrix file that needs to be split into smaller and more manageable sections.
2. This file was split into 123 files of (~320 Mb each) containing XX rows of data. This could be performed in a number of ways, including the BASH one-liner below.
```
awk -v lines=24244 'NR==1 || (NR-1)%lines==0 { file = "file_" ++count ".csv" } { print > file }' large_output_file.csv
```
3. Run the `cosine.sh` SLURM script to apply the `min-avg_cosine.py` script to each split file in parallel.
4. The `mincosine.ipynb` Jupyter Notebook can then pool the output files into a single file and provide the number of BGCs that were novel (min d>0.2).

[^1]: Paoli, L. et al. Biosynthetic potential of the global ocean microbiome. Nature 607, 111-118, doi:10.1038/s41586-022-04862-3 (2022).