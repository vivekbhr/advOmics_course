# Single-cell sortChIC analysis to identify differential chromatin modifications in mouse bone marrow

## System setup

 - **Machine:** At least 16 gb of RAM.
 - **OS:** mac/linux. A linux subsystem for windows might work, but is untested.
 - **Language:** Python, Bash, R (optional, if you use DESeq2 instead of pyDESeq2)
 - **Slack channel name:** project_sortchic

### Install conda

Install the minimal version of **conda**, either as [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main), or from [miniforge](https://github.com/conda-forge/miniforge)

You can read all about conda on its [User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).

In a nutshell:
   - Add your favorite channels to the list of channels where conda searches for packages conda config --add channels bioconda
   - Create an environemnt specific for your project (e.g. you want to use an older version of python) conda create -n my_project python=3.7
   - Activate it: conda activate my_project
   - Install tools you want to use (notice the change in prompt, prefixed with (my_project) to show we are in the activated environment (my_project)$ conda install -c bioconda bwa samtools
   - Deactivate the environment (my_project)$ conda deactivate

### install required packages via conda

Create a new conda environment: advomics, and install:

 - sincei (also installs scanpy)
 - jupyter-notebook
 - loompy
 - pyDESEQ2

**Launch a new notebook within the advomics env and load the above packages to ensure that everything works.**


### Tips

 - Loompy is not well maintained, and might be incompatible with the latest version of scanpy/anndata. It might need a separate environment. I suggest installing sincei and jupyter first, then trying to load loompy and the dataset. If you get an error, try a fresh install of loompy in a new environment.

## Objectives

The aim of this project is to partially replicate the findings of [Zeller, Yueng et al (2023)](https://www.nature.com/articles/s41588-022-01260-3), using a combination of tools: sincei, scanpy, and pyDESeq2.

In particular we will try to replicate the results of cell clustering (shown in **Fig. 2 a,b**) and recover the top cell-type specific genes regulated by the activate/repressive chromatin states (**Fig. 2 c-f**).  

Each project group would be assigned a chromatin mark and would therefore focus on replicating the results shown for their own chromatin mark.

 - Group 1: H3K4me1
 - Group 2: H3K27me3
 - Group 3: H3K4me3
 - Group 4: H3K9me3

## Instructions

 - The data is available from figshare, using this link: https://figshare.com/s/31c3d4676830a614d886. Have a look at the data description on figshare for more details about the files.

 - Download the .loom file named with your selected hPTM (in case of H3K4me1, the file is `scCounts_genebodies_k4me1.loom`)

 - Follow the [sortChIC data analysis tutorial from sincei](https://sincei.readthedocs.io/en/latest/content/tutorials/sincei_tutorial_sortChIC.html) (step 4 onward) to filter and cluster the cells using the gene-body signal of the selected hPTM.

 - After clustering, load the clustered .loom file in python, and apply the T-test in scanpy to find the top genes per cluster. Annotate clusters based on the top genes (compared to those shown in the sortChIC manuscript)


### Expected analysis for the project

 - Compare your clustering/annotation results with the metadata provided with the files. Do we see the same cell types as the authors, or are there differences? How can we explain the differences (is that related to the biology, or an effect of preprocessing)?

 - Identify genes that are driving the cluster-specific differences. First we will try a simpler method: T-test or Wilcox test. Next, we will try the bulk differential enrichment (DE) analysis using DESeq2. Next, we will compare the 2 methods and also compare the results with the original study. Which genes are (re)discovered? Could you identify novel genes? Which analysis works better?

### Tips

 - The logic for pyDESeq2 analysis in our data is same as that for differential expression analysis (introduced in the **Transcriptomics** lectures). And if you prefer, you can use the original version of DESeq2 implemented in R, instead of pyDESeq2, and follow the workflow of the transcriptomics practical.

 - As pyDESEq2 is aimed at (pseudo)bulk differential expression analysis, this would require creating arbitrary pseudo-bulk replicates from single-cell counts for each cluster/celltype by sampling cells (without replacement) from each group.
