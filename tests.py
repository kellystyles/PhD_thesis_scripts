#!/usr/bin/env python
import subprocess
import tempfile


def exec_notebook(path):
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        print(fout.name)
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=1000",
                "--output", fout.name, path]
        subprocess.check_call(args)


def test():
# these notebook must be run in order
    exec_notebook('./Chapter_3_endolichenic_fungi_isolation/lichen_plots.ipynb')
    exec_notebook('./Chapter_3_endolichenic_fungi_isolation/statistical_tests.ipynb')
    exec_notebook('./Chapter_4_fungiflow_workflow/pipeline_tests.ipynb')
    exec_notebook('./Chapter_4_fungiflow_workflow/ELF_data.ipynb')
    exec_notebook('./Chapter_4_fungiflow_workflow/comparative_phylogenomics.ipynb')
    exec_notebook('./Chapter_5_BGCs_from_endolichenic_fungi/BGC_data_plots.ipynb')
    exec_notebook('./Chapter_5_BGCs_from_endolichenic_fungi/cosine_analysis_scripts/cosine_plots.ipynb')
    exec_notebook('./Chapter_5_BGCs_from_endolichenic_fungi/bioassay_plots.ipynb')
    exec_notebook('./Chapter_6_targeted_BGC_identification/cosine_analysis_scripts/iterative_pyhmmer_build.ipynb')
    exec_notebook('./Chapter_6_targeted_BGC_identification/identify_IDT_genomes.ipynb')

if __name__ == '__main__':
    test()
