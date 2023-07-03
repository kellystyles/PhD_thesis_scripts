# Cosine analysis

**All the scripts and Notebooks discussed here were written by Dr Peng Hou and were modified for clarity by myself.**

One way to assess BGC novelty is to compare newly identified BGCs to known BGCs, both characterized and uncharacterized. 

The ELF BGCs were input into BiG-SLiCE and compared to the 1.2 million BGCs in the BiG-FAM BGC database. The BiG-SLiCE Euclidean distances were converted to cosine vectors as performed in the Paoli et a., 2022[^1] that examined the novelty of BGCs from a global ocean microbiome. Various plots were prepared using the Jupter Notebook, `cosine_plots.ipynb`.

The steps to obtain to obtain the cosine novelty scores, the scripts were run in the following order:
1. Run BiG-SLiCE using the `run_big-slice.sh` SLURM script. This took ~2.5 days to run due to the large number of comparisons required. This will output a large matrix file that needs to be split into smaller and more manageable sections.
2. This file was split into 123 files of (~320 Mb each) containing 24,244 rows of data. This could be performed in a number of ways, including the BASH one-liner below.
```
awk -v lines=24244 'NR==1 || (NR-1)%lines==0 { file = "file_" ++count ".csv" } { print > file }' large_output_file.csv
```
3. Run the `cosine.sh` SLURM script to apply the `min-avg_cosine.py` script to each split file in parallel.
4. The `mincosine.ipynb` Jupyter Notebook can then pool the output files into a single file and provide the number of BGCs that were novel (mininum distance score >0.2).
5. The script `parse_cosine_1.py` prepares a list of only ELF BGCs and the associated PFAM hits. The script `parse_cosine_2.py` will then perform agglomerative clustering using scikit-learn, to group BGCs into GCFs and GCCs, using a threshold of 0.2 and 0.8, respectively. Both scripts can be run using the SLURM script `parse_cosine.sh`.
6. The `metatable.py` script will parse the SQLite 3 database to obtain more nuanced details of each BGC. Run as so:
```
python3 metatable.py ELF.db ELF
```
7. The `metatable.ipynb` Jupyter Notebook then merges the output of `metatable.py` and the mincosine scores as prepared from the concatenated output CSV of `mincosine.ipynb`. 
8. The Jupyter Notebook `cosine_plots.ipynb` then concatenates all the relevant data and generates plots.


[^1]: Paoli, L. et al. Biosynthetic potential of the global ocean microbiome. Nature 607, 111-118, doi:10.1038/s41586-022-04862-3 (2022).
