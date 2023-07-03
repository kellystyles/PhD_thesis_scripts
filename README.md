# PhD Jupyter Notebooks

This GitHub is the repository for all the figures generated for my PhD thesis discussing the genetics and secondary metabolism of endolichenic fungi (ELF) from New Zealand.
Each chapter in my thesis has its own sub-directory, within which is contained one to two Jupyter Notebooks that will perform the data analyses and generate plots from the data.

## Dependencies

The Notebooks in this repository use Python 3.9. Aside from standard Python libraries, the required dependencies used are Python packages listed below:
| Package | Version | Reference |
|---------:|:---------:|-----------:|
| pandas | 1.5.3 | McKinney, 2010[^McKinney,2010] |
| seaborn | 0.11.2 | Waskom, 2021[^Waskom,2021] |
| plotly | 5.13.1 | Plotly Technologies Inc., 2015[^PlotlyTechnologiesInc.,2015] |
| prince | 0.8.3 | Halford, 2022[^Halford,2022] |
| scipy | 1.10.1 | Virtanen et al., 2020[^Virtanenetal.,2020] |
| scikit-learn | 1.2.1 | Pedregosa et al., 201[^Pedregosaetal.,2011] |
| upsetplot | 0.8.0 | Nothman, 2018[^Nothman,2018] |
| statannotations | 0.5 | Charlier et al., 2022[^Charlieretal.,2022] |
| pyhmmer | 0.7.1 | Larralde & Zeller, 2023[^Larralde&Zeller,2023] |

To open these Notebooks, I would recommend using Anaconda Navigator to install the dependencies and Jupyter Notebook, from which these Notebooks can be run.
[^McKinney,2010]: McKinney, W. (2010). Data structures for statistical computing in python. Proceedings of the 9th Python in Science Conference
[^Waskom,2021]: Waskom, M. L. (2021). Seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021. 
[^PlotlyTechnologiesInc.,2015]: Plotly Technologies Inc. (2015). Collaborative data science. In. Montréal, QC: Plotly Technologies Inc.
[^Halford,2022]: Halford, M. Prince. In. https://github.com/MaxHalford/prince 
[^Virtanenetal.,2020]: Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., & Bright, J. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261-272. 
[^Pedregosaetal.,2011]: Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., ... & Duchesnay, É. (2011). Scikit-learn: Machine learning in Python. the Journal of machine Learning research, 12, 2825-2830.
[^Nothman,2018]: Nothman, J. (2018). UpSetPlot. In. https://github.com/jnothman/UpSetPlot.
[^Charlieretal.,2022]: Charlier, F., Weber, M., Izak, D., Harkin, E., Magnus, M., Lalli, J., Fresnais, L., Chan, M., Markov, N., & Amsalem, O. (2022). trevismd/statannotations: v0. 5. Zenodo10, 5281. 
[^Larralde&Zeller,2023]: Larralde, M., & Zeller, G. (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad214

## Overview of each Chapter

### Chapter 3 - ELF isolation, sequencing, and taxonomy

This chapter discusses the breadth of the lichen collection surveyed, as well as the ELF isolated from some of these lichens. Several statistical tests (Pearson's Chi-Squared Mean) are performed in an attempt to link the observation of ELF taxa to lichen taxa, form, and macroclimatic factors.

How to run:
1. `lichen_plots.ipynb` will prepare descriptive plots of the lichen and ELF dataset using both seaborn[^Waskom,2021] and plotly[^PlotlyTechnologiesInc., 2015], by parsing data from the `ELF_master_results.csv` file.
2. `statistical_tests.ipynb` will perform the statistical analyses found in this chapter. The lichen and ELF dataset, from `ELF_master_results.csv`, are parsed and data manipulated to perform Pearson's Chi-Squared tests using the scipy[^Virtanenetal.,2020] package. An abundance table is also prepared from the taxonomy data, from which rarefaction curves can be prepared using [iNEXT Online](https://chao.shinyapps.io/iNEXTOnline/).

### Chapter 4 - ELF assembly metrics

In this chapter, aspects of the [Fungiflow](https://github.com/kellystyles/fungiflow) pipeline are examined. Test results of the pipeline against several datasets are inspected and the ELF assembly quality are examined. Aspects of ELF assembly metrics and quality are plotted and discussed. Finally, comparative phylogenomics of the ELF is performed and plots prepared.

How to run:
1. `pipeline_tests.ipynb` will prepare plots of the results of the Fungiflow pipeline tests against taxonomically diverse synthetic genomes. Data is parsed from `chapter_4_fungiflow_workflow/master_results.csv`.
2. Assembly metrics from `ELF_master_results.csv` are plotted and compared using `ELF_data.ipynb`. It then continues and plots the estimated and actual coverages of the two short read Illumina sequencing runs. It finishes with a Prinicipal Components Analysis of various assembly metrics to identify correlated metrics, using two separate PCA methods (prince[^Halford,2022] and sklearn[^Pedregosaetal.,2011]).
3. Comparative phylogenomics are performed from the output CSV files of `funannotate compare`[^Palmeretal.,2021] using the `comparative_phylogenomics/comparative_phylogenomics.ipynb` Jupyter Notebook. This parses each of the relevant summary CSV files from `funannotate compare` and generates a number of bar plots and clustered heatmaps.

[^Palmeretal.,2021]: Palmer, J., & Stajich, J. (2021). Funannotate v1. 8.3: eukaryotic genome annotation (Version 1.8. 3). Zenodo. doi, 10. 

### Chapter 5 - ELF BGCs

Chapter 5 discusses the diversity and novelty of ELF bisoynthetic gene clusters (BGCs). Most of the plots here describe the type and diversity of BGC data captured by the Fungiflow pipeline.
*Note that the cosine analysis of BiG-SLiCE data is only represented by the scripts used, as the data is too large (38.9 GB) to be uploaded here and takes several days to process.*

How to run:
1. Run `BGC_dataplots.ipynb` to prepare descriptive plots of BGCs from the ELF, parsed from `ELF_all_bgcs.csv`.
2. Run `cosine_plots.ipynb` to generate plots that assess the novelty of ELF BGCs from the BiG-SLiCE data prepared by the cosine analysis scripts (*does not perform the actual cosine analysis*).
3. `bioassay_plots.ipynb` generates plots of preliminary bioassay data from the crude extracts of several ELF. It uses the statannotations[^Charlieretal.,2022] and scipy[^Virtanenetal.,2020] packages to determine the statistical significance of various ELF to the negative control.

### Chapter 6 - targeted BGC identification

Chapter 6 contains scripts that outline the search for BGCs using profile Hidden Markov Models (pHMMs). It describes two pHMM sets, a curated set found in the `curated_pHMMs` directory and an iteratively generated pHMM set found in the `SwissProt_pHMMs` directory. These searches focus on the fungal indole diterpene (IDT) pathway.

*There are 3 versions of the `cluster_search.py` script, with each later version being an improved version of the previous.*

How to run:
1. Run Jupyter Notebook `iterative_pyhmmer_build.ipynb` to obtain iterative SwissProt models using pyhmmer[^Larralde&Zeller,2023] and descriptive plots of this process.
2. Download all fungal genomes with annotated assemblies from GenBank using the script `genome_getter.sh`.
3. BLASTp each genome against a BLAST DB of curated IDT protein sequences with `genome_checker.sh` script.
  *multiFASTAs to generate each IDT BLAST DB are the same as the ones used to generate the curated pHMMs.*
4. Prepare test databases using the `prepare_test_dbs.py` script.
    ```
    python3 prepare_synthetic_dbs.py 100
    ```
5. Each set of pHMMs was applied using the `cluster_search_v3.sh` script.
  *run this script once for each pHMM set*
   ```
   python3 cluster_search_v3.py -g 'gbk_dir' -p 'phmm_dir'
   ```
7. Plots were generated using the Jupyter Notebook `identify_IDT_genomes.ipynb`.
