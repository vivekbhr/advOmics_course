# Single-cell sortChIC analysis to identify differential chromatin modifications in mouse bone marrow

## System setup

 - **Machine:** At least 16 gb of RAM.
 - **OS:** mac/linux. A linux subsystem for windows might work, but is untested.
 - **Language:** Python, Bash, R (optional, if you use DESeq2 instead of pyDESeq2)
 - **Slack channel name:** project_sortchic

## Introduction

Organism development is a complicated process in which cells follow defined differentiation paths to form each cell type. Epigenetic mechanisms such as histone modifications have been shown to have a crucial role in controlling differentiation decisions, but studying the changes in histone mark distributions during the process of cell differentiation has so far been very challenging.

To tackle this challenge, [Zeller, Yueng et al (2023)](https://www.nature.com/articles/s41588-022-01260-3) developed a method that utilizes detailed knowledge about cell-type enrichment and identification by surface markers, together with highly sensitive histone mark profiling in single cells, called ‘sort-assisted single-cell chromatin immunocleavage’ (sortChIC). For each cell, they recorded its histone modification by sequencing antibody-targeted micrococcal nuclease cuts, as well as surface marker abundance by fluorescence-activated cell sorting (FACS). They applied this method to mouse bone marrow, for which antibody combinations have been established for the identification and purification of many mature and progenitor cell types, down to very rare hematopoietic stem cells (HSCs) that make up less than 0.1% of bone marrow cells (Fig. 1a). This enabled them to systematically describe the dynamic changes of two active histone modifications (H3K4me1 and H3K4me3) and two repressive histone modifications (H3K9me3 and H3K27me3) along the full hematopoietic differentiation trajectory. They used histone-modification abundance on known marker genes, as well as FACS-recorded surface-marker information, to unambiguously identify cell types for all four histone modifications (Fig. 1b).

![sortChIC background](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41588-022-01259-w/MediaObjects/41588_2022_1259_Fig1_HTML.png?as=webp)

## Objectives

The aim of this project is to partially replicate the findings of [Zeller, Yueng et al (2023)](https://www.nature.com/articles/s41588-022-01260-3), using a combination of tools: sincei, scanpy, and pyDESeq2.

In particular we will try to replicate the results of cell clustering (shown in **Fig. 2 a,b**) and recover the top cell-type specific genes regulated by the activate/repressive chromatin states (**Fig. 2 c-f**).  

Each project group would be assigned a chromatin mark and would therefore focus on replicating the results shown for their own chromatin mark.

 - Group 1: H3K4me1
 - Group 2: H3K27me3
 - Group 3: H3K4me3
 - Group 4: H3K9me3

## Installation instructions

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

Create a  new conda environment: advomics, and install:

 - jupyter
 - loompy
 - sincei (also installs scanpy)
 - pyDESEQ2 (from pip)

**Launch a new notebook within the advomics env and load the above packages to ensure that everything works.**


### Tips
 - installing packages from conda requires a bit of back and forth. Not all packages will install right away, as the fetched package binaries depend on your OS. For example, the following order worked on an M3 macbook (architecture: OSX-ARM64)

 ```
 conda create -n advomics_sincei -c bioconda -c conda-forge jupyter
 conda activate advomics
 pip install sincei
 pip install pydeseq2
 ```

### Getting around loompy

 - Loompy is not well maintained, and might be incompatible with the latest version of scanpy/anndata. It might need a separate environment. I suggest installing sincei and jupyter first, then trying to load loompy and the dataset. If you get an error, try a fresh install of loompy in a new environment.

We created a version of sincei that doesn't rely on the loompy and loom files and instead uses the .h5ad files for the data (available on figshare). You can install it via github:

```
conda create -n advomics2 -c conda-forge scanpy gensim
conda activate advomics2
(sincei): pip install git+https://github.com/bhardwaj-lab/sincei.git@vivek#egg=sincei
# downgrade scipy/numpy (otherwise gensim might complain)
pip install scipy==1.12
# test
scClusterCells --help
```

### Installing multiqc

To install the version of multiqc that works with sincei (and summarises the QC results), you can try:

```
pip install --upgrade --force-reinstall git+https://github.com/vivekbhr/MultiQC.git
```


## Analysis instructions

 - The data is available from figshare, using this link: https://figshare.com/s/31c3d4676830a614d886. Have a look at the data description on figshare for more details about the files.

 - Download the .loom/.h5ad file named with your selected hPTM (in case of H3K4me1, the file is `scCounts_genebodies_k4me1.loom`)

 - Follow the [sortChIC data analysis tutorial from sincei](https://sincei.readthedocs.io/en/latest/content/tutorials/sincei_tutorial_sortChIC.html) (step 4 onward) to filter and cluster the cells using the gene-body signal of the selected hPTM.

 - After clustering, load the clustered .loom/.h5ad file in python, and apply the T-test in scanpy to find the top genes per cluster. Look at the description of [anndata file format](https://anndata.readthedocs.io/en/stable/tutorials/notebooks/getting-started.html) to understand how to work with and extract information from this object. Since you start from raw counts, you need to first normalize and log-transform your data. Have a look at the [normalization](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#normalization) and [marker gene set]((https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#differentially-expressed-genes-as-markers) steps from the scanpy manual for this analysis. Annotate clusters based on the top genes (compared to those shown in the sortChIC manuscript).

 - You can use the file `allgenes_symbols_biotypes.tsv` to find the gene symbols corresponding to the ensembl gene IDs indicated in the top regions of your clusters. Hopefully the top ~20 gene list have some of the same genes are the authors show in Fig 2.


### Expected analysis for the project

 - Compare your clustering/annotation results with the metadata provided with the files. Do we see the same cell types as the authors, or are there differences? How can we explain the differences (is that related to the biology, or an effect of preprocessing)?

 - Identify genes that are driving the cluster-specific differences. You can start with a simpler method: T-test or Wilcox test. Next, you can try the bulk differential enrichment (DE) analysis using DESeq2. The idea is to compare the 2 methods and also compare the results with the original study. Which genes are (re)discovered? Could you identify novel genes? Which analysis works better?

### Tips

 - If you get errors with the memory usage, you can subset your data for a single chromosome using the feature names in the .loom objects (which are in the format `chromosome_start_end::geneID`).

 - The logic for pyDESeq2 analysis in our data is same as that for differential expression analysis (introduced in the **Transcriptomics** lectures). And if you prefer, you can use the original version of DESeq2 implemented in R, instead of pyDESeq2, and follow the workflow of the transcriptomics practical.

 - As pyDESEq2 is aimed at (pseudo)bulk differential expression analysis, this would require creating arbitrary (3-5) pseudo-bulk counts from single-cell counts for each cluster/celltype by sampling cells (without replacement) from each of them. Recall from the single-cell genomics lecture that the pseudo-bulk DE methods should work very well for such analysis.

 - H3K4me1 and H3K4me3 are **active** histone modifications, therefore they correlate with active transcription, while H3K27me3 and H3K9me3 are **silencing** histone modifications, enriched on silenced genes. Therefore, we expect *opposite* results for group 1/3 vs group 2/4.
