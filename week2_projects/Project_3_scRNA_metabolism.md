# Inferring metabolism from gene expression

## System setup

- **Machine:** At least 8 gb of RAM.
- **OS:** Windows/Mac/Linux.
- **Language:** R. Tested on version 4.4.2.
- **Slack channel name:** project_metabolism_from_gene_expression

## Introduction and dataset

Metabolomics has important limitations, such as lower throughput compared to transcriptomics, and the fact that single-cell analysis is only available for non-polar metabolites. As such, important processes such as glycolysis cannot be directly measured at the single-cell level using metabolomics.

Somewhat surprisingly, it is possible to infer metabolic information of transcriptomic data. This is because cellular metabolism is regulated by enzymes, of which the transcripts can be detected with RNAseq. For example, one can use a gene set of glycolytic enzymes to infer glycolytic activity. A major limitation though is that there's only a weak correlation between transcript counts and enzymatic activity.

## Objective

In this exercise, we will support the main conclusions of one paper using the data of another, completely unrelated paper. This should give you an idea of how powerful gene set analysis can be, and should prepare you to apply this technique yourself in your own research projects.

The dataset paper that we will look at is the paper of Haber et al., 2017. The paper can be viewed at https://doi.org/10.1038/nature24489 . It is a widely used single-cell RNAseq dataset. According to the editorial summary, it is "a comprehensive analysis of the epithelial cell composition of mouse small intestines when healthy and after infection". In this assignment, we will only use data of the healthy mouse. We will start by downloading, loading and pre-processing the data.

The results paper that we will look at is the paper of Rodriguez-Colman et al., 2017, "Interplay between metabolic identities in the intestinal crypt supports stem cell function" ( https://www.nature.com/articles/nature21673 ). This paper, written by our lab, describes the metabolic differences in cell types in the intestinal crypt. The abstract contains a few conclusions, which we will try to replicate using the data of the Haber et al. paper.

## Obtaining the data

First, go the Haber et al. 2017 paper linked above, and locate the file "GSE92332_atlas_UMIcounts.txt.gz" within the provided source data. Warning: the download is just 15.2 MB, but if you extract it, you'll end up with a text file of 226 MB. For this assignment it doesn't matter whether you leave it compressed or not. To save disk space, you can leave it compressed, so that the file name ends with ".txt.gz" instead of ".txt".

Use `data_table <- read.table(data_file_path, header = TRUE, sep = "\t", row.names = 1)` to load the data, which is a large table of raw transcript counts.The data should have 15 971 rows (genes) and 7216 columns (cells). We're going to do this assignment in Seurat, so make sure it's installed, and load it using `library(Seurat)`. Put your data in a Seurat object with `data <- CreateSeuratObject(data_table)`.

If you're not familiar with Seurat, the PMBC3K introductory tutorial at https://satijalab.org/seurat/articles/pbmc3k_tutorial is really helpful to follow. The documentation of Seurat (same website) is very comprehensive. Because Seurat is so widely-used ,AI tools should also provide useful guidance in teaching you how to use the package.

## Filtering for epithelial cells
A straightforward way to filter for epithelial cells is to look at the expression of epithelial markers. Use `data_subset <- subset(data, Krt8 > 0 | Epcam > 0, slot = "counts")` to filter for cells that have expression of either Krt8 or Epcam. What percentage of cells have expression of either Krt8 or Epcam? Do you think filtering for epithelial cells is necessary for this dataset?

## Quality checks

In the PMBC3K Seurat tutorial at https://satijalab.org/seurat/articles/pbmc3k_tutorial , the authors recommend to calculate the percentage of mitochondrial genes using `data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")`. Why will this approach not work for this dataset?

Now read the relevant part of the methods section (we are just assuming you hadn't done that yet), and describe how they did the filtering.

## Normalization, scaling, variable features, cell type assignment

Normalize and log-transform the data like Haber *et al*. Find the 2000 most variable genes. Which are the top 10 most variable genes? The 2000 features will be used for dimensionality reduction in the next steps. Scale and center the data with the default settings.

We're also going to assign the cell types. We could in principle do this ourselves like in the tutorial, so by doing automated clustering and looking at the marker genes of each cluster. Fortunately, the authors already did this for us. If you do `head(colnames(data))`, you'll see that cells have names like `B1_AAACCGTGCCAAGT_Tuft`. The first part is the batch, the second part the barcode of the cell, and the third part the cell type, as assigned by the authors. Seurat automatically extracted the first part, assuming it would be cell types, into `Idents(data)`. You can see all cell types using `levels(Idents(data)))`, which should currently still give something like "B1", "B10", "B2", etc.

The code for assigning the correct cell types into `Idents(data)` is a bit tricky, involving functions like `str_split` or similar. Still, it is necessary to perform this step before procedding. Extract the cell types, and verify that `levels(Idents(data))` now shows the cell types instead of the batch numbers.

## Dimensionality reduction, visualization and gene set scores

Now that all the pre-processing is over, we can analyze the data. Run and visualize a PCA, then a UMAP.

On Wikipedia, there's a list of intestinal stem cell marker genes: https://en.wikipedia.org/wiki/List_of_intestinal_stem_cell_marker_genes . Plot a few of them per cell type, keeping in mind that mouse gene names might differ slightly from the human gene names listed on Wikipedia, and that Seurat is case-sensitive.

Now we're going to leave the PMBC3K Seurat tutorial. For the rest of the assignment, we will make extensive use of gene sets and the `AddModuleScore` function. Using the Seurat documentation, describe in your own words what it does.

Example usage of this function is as follows:

```
epithelial_gene_set <- c("Epcam", "Krt8")
data <- AddModuleScore(data, features = epithelial_gene_set, name = "epithelial_gene_set")
FeaturePlot(data, features = "epithelial_gene_set1", reduction = "umap")
```
The one gotcha here is that you need to append "1" to the name in `FeaturePlot`. This is because `AddModuleScore` can in principle also take multiple gene sets at once, and it will append a number to the name for each gene set.

Now use `AddModuleScore` with the stem cell genes you found above.

The benefit of using gene sets over single genes, is that it is less noisy. If sequencing failed to detect any `Lgr5` transcripts in a cell, either because of technical issues or because the cell was not expressing Lgr5 at that moment, the cell can still be detected as a stem cell if it expresses enough other stem cell markers.

At https://www.gsea-msigdb.org/gsea/msigdb/index.jsp , there are large collections of gene sets, for both humans and mice. On that page, click on "M2" ("curated gene sets", part of "Mouse Collections"). This brings you to a table (which also includes other mouse collections). The curated gene sets are subdivided into various collections. Find the "CP: Canonical pathways" collection, and click on the "browse 1730 gene sets" (number might be updated in the future). This brings you a page with a long list of gene sets.

Use Ctrl + F to search for the gene set "REACTOME_SIGNALING_BY_NOTCH". Click through to the page of this gene set, and download the gene set as "TSV metadata". Open the file in a text editor or Excel, and copy the gene names to a vector in R. (Even better would be to write a function that reads in the file and extracts the gene names, but this is a bit more advanced.)

Then use `AddModuleScore` with this gene set, and make a `FeaturePlot` like before. The color scaling can be a bit unfortunate, since color scaling is optimized for genes. By adjusting the scaling a bit, you can obtain nice results, like in this example:

```
FeaturePlot(data, features = "reactome_signaling_by_notch1", reduction = "umap", min.cutoff = 0.1, max.cutoff = 0.4)
```

There's also the function `VlnPlot` to visualize gene sets, which uses the cell types we've assigned earlier. Example usage of this function, if you just want to compare Paneth and stem cells, is as follows:

```
VlnPlot(data, features = c("reactome_signaling_by_notch1", "reactome_signaling_by_wnt1"), assay = "RNA", idents = c("Paneth", "Stem"))
```

Here, we also included a gene set for Wnt signaling, which you can find on the same page as the Notch signaling gene set.

As stated in the introduction of this assignment, we were interested in inferring metabolism, so we're interested in gene sets like `REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE` or `WP_ELECTRON_TRANSPORT_CHAIN`.


## A little research project

To show the power of using gene sets, we (or rather: you) are going to replicate some of the conclusions of another paper, which was published in the same year as the Haber et al. paper. This paper, about research by the UMC Utrecht, is called "Interplay between metabolic identities in the intestinal crypt supports stem cell function" ( https://www.nature.com/articles/nature21673 ). Read the abstract, list all the conclusions, and see which ones you can investigate with the Haber et al. dataset, and do so for a few of them.

To pass the assignment, you just need to find supporting (or contrasting) evidence (one gene set each) for at least two of the conclusions of the paper. You can also try to find evidence for more conclusions, but this is not necessary. You can use gene sets from any source, but in principle using only gene sets on the page you went to before should be enough. You can also compare multiple gene sets, to check whether they give similar results.

Here, we used `AddModuleScore` to calculate gene set scores. An alternative method is `UCell`. What is the difference between both methods? Could switching to `UCell` affect your conclusions from the previous paragraph?

We used gene set scores here, which is the most basic method of metabolic interference from transcriptomics. There are also more advanced tools such as Compass from Wagner *et al*. We considered that tool too computationally expensive for this assignment. Based on their paper/documentation, what is the advantage of using Compass over using gene scores? Could using Compass affect your conclusions?

And one final question: you'll probably notice that the differences between cell types are present, but quite small. Why could this be?
