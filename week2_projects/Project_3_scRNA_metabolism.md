# Inferring metabolism from gene expression

## System setup

- **Machine:** At least 8 gb of RAM.
- **OS:** Windows/Mac/Linux.
- **Language:** R. Tested on version 4.4.2.
- **Slack channel name:** project_metabolism_from_gene_expression

## Introduction and dataset

It seems you should use transcriptomics if you want to look at transcripts, and metabolomics to look at metabolic processes. While true, you can also use (single cell) transcriptomics to look at metabolic processes. While metabolites do not have RNA transcripts, cellular metabolism is regulated by enzymes, which do have transcripts that you can analyze. The benefit of using single cell transcriptomics to analyze metabolic processes, is that you also have information about other aspects of the cell, such as the cell type. This can give you a lot more complete picture of the tissue you are looking at.

## Objective

In this exercise, we will support the main conclusions of one paper using the data of another, completely unrelated paper. This should give you an idea of how powerful gene set analysis can be, and should prepare you to apply this technique yourself in your own research projects.

The dataset paper that we will look at is the paper of Haber et al., 2017. The paper can be viewed at https://doi.org/10.1038/nature24489 . It is a widely used single-cell RNAseq dataset. According to the editorial summary, it is "a comprehensive analysis of the epithelial cell composition of mouse small intestines when healthy and after infection". In this assignment, we will only use data of the healthy mouse. We will start by downloading, loading and pre-processing the data.

The results paper that we will look at is the paper of Rodriguez-Colman et al., 2017, "Interplay between metabolic identities in the intestinal crypt supports stem cell function" ( https://www.nature.com/articles/nature21673 ). This paper, written by our lab, describes the metabolic differences in cell types in the intestinal crypt. The abstract contains a few conclusions, which we will try to replicate using the data of the Haber et al. paper.

## Downloading the data

First, go the Haber et al. 2017 paper linked above, find the "Accession codes" section and click on the link to the Gene Expression Omnibus. (The text of that link reads "GSE92332", and it should be the only link in that section.)

If you scroll down on that webpage, you should find a table with download links. Download the "GSE92332_atlas_UMIcounts.txt.gz". Warning: while the download is just 15.2 MB, but if you extract it, you'll end up with a text file of 226 MB. For this assignment it doesn't matter whether you leave it compressed or not. To save disk space, you can leave it compressed, so that the file name ends with ".txt.gz" instead of ".txt".

For this assignment, you'll need to hand in a single PDF file containing *all* your code, the requested figures and the text. It should read like a report. The best way to do this is to write a RMarkdown file, which is really easy to do from RStudio: https://www.geeksforgeeks.org/generate-pdf-from-rmarkdown-file-with-r/ (section: "Using the RStudio IDE").

To verify whether everything works, add the following contents to your RMarkdown file, and run it:


```
data_file_path <- "C:/example/path/to/GSE92332_atlas_UMIcounts.txt.gz"
# Change this path to your actual path. If you place the file in the same folder as the RMarkdown file, you can simply use "GSE92332_atlas_UMIcounts.txt.gz"

if (!file.exists(data_file_path)) {
  stop(paste("File not found: ", normalizePath(data_file_path)))
}
print("File found!")

```

Make sure you have no errors, and that the message "File found!" is printed. If you have errors, you'll need to fix them before you can continue with the assignment.

## Loading in the data
Use `data_table <- read.table(data_file_path, header = TRUE, sep = "\t", row.names = 1)` to load the data, which is a large table of raw transcript counts. Use `dim(data)` to check the dimensions of the data. The data should have 15 971 rows (genes) and 7216 columns (cells). Verify with `rownames(data_table)` and `colnames(data_table)` that the row names are indeed gene names and the column names are cell names.

We're going to do this assignment in Seurat. If you haven't installed Seurat yet, run `install.packages('Seurat')` in the R console. Then, back in the RMarkdown document, load the Seurat library with `library(Seurat)`, and put your data table in a Seurat object with `data <- CreateSeuratObject(data_table)`. This step is necessary to use the Seurat functions on your data.

## Filtering for epithelial cells
A straightforward way to filter for epithelial cells is to look at the expression of epithelial markers. Use `data_subset <- subset(data, Krt8 > 0 | Epcam > 0, slot = "counts")` to filter for cells that have expression of either Krt8 or Epcam. What percentage of cells have expression of either Krt8 or Epcam? Do you think filtering for epithelial cells is necessary for this dataset?

## Quality checks

In the PMBC3K Seurat tutorial at https://satijalab.org/seurat/articles/pbmc3k_tutorial , the authors recommend to calculate the percentage of mitochondrial genes using `data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")`. This percentage can then be visualized across all cells using `VlnPlot(data, features = c("percent.mt"))`. What goes wrong when using this approach in our dataset?

Tip: with `mito.genes <- grep(pattern = "^MT-", x = rownames(data), value = TRUE)`, you can check which genes match the given pattern, "^MT-" in this case. Because of the "^", this pattern matches all gene names that start with "MT-".

Check the result of `VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)`. How can you see that the data has already been filtered?

Now read the relevant part of the methods section (because of course you haven't done that yet), and describe how they did the filtering.

## Normalization, scaling, variable features, cell type assignment

Normalize and compute the log-transformation of the data as described in the Methods section. You'll only need one statement containing the `NormalizeData` function for this. https://satijalab.org/seurat/reference/normalizedata

Next, find the most variable genes using `data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)`. Which are the top 10 most variable genes? You can find how to do this in the PMBC3K Seurat tutorial linked above. These features will be used for dimensionality reduction in the next steps.

Scale the data using `data <- ScaleData(data)`. This step is necessary for the PCA and UMAP steps.

We're also going to assign the cell types. We could in principle do this ourselves, by doing automated clustering and looking at the marker genes of each cluster, but the authors already did this for us. If you do `head(colnames(data))`, you'll see that cells have names like `B1_AAACCGTGCCAAGT_Tuft`. The first part is the batch, the second part the barcode of the cell, and the third part the cell type, as assigned by the authors.

Seurat stores cell types (or cell identities) in `Idents(data)`. You can see all cell types using `levels(Idents(data))`. If you execute that expression now, you'll see that the "cell types" are `"B1"  "B10" "B2"  "B3"  "B4"  "B5"  "B6"  "B7"  "B8"  "B9"`. These are actually the batch numbers, which Seurat automatically extracted  from the cell names. We want to store the cell types in `Idents(data)` instead. (For simplicity, we'll skip batch correction in this assignment.) This code can be a bit tricky, so we'll give you the code for this step:

```
# This splits the cell names by "_", and takes the third part of the splitted name
cell_types_by_cell <- sapply(strsplit(colnames(data), "_"), function(x) x[3])
Idents(data) <- factor(cell_types_by_cell)  # Store these as factors in Idents(data)
```

Verify that `levels(Idents(data))` now shows the cell types instead of the batch numbers.

## Dimensionality reduction, visualization and gene set scores

Now that all the pre-processing is over, we can analyze the data. We're going to do a PCA analysis. For this section, you'll need to find and follow the appropriate steps of the PMBC3K Seurat tutorial linked above. You'll need to at least use the `RunPCA` and `DimPlot` functions.

How many principal components do you roughly need to capture the majority of the true signal in the data? You can use the `ElbowPlot` function to visualize this.

Now use this number in the `RunUMAP` function, and plot a UMAP using `DimPlot`. (We can skip the clustering part of the tutorial - the part with `FindNeighbors` and `FindClusters` - because we already assigned cell types ourselves.)

On Wikipedia, there's a list of intestinal stem cell marker genes: https://en.wikipedia.org/wiki/List_of_intestinal_stem_cell_marker_genes . Plot a few of them using `FeaturePlot`. Keep in mind that the dataset might use slightly different gene names than the ones on Wikipedia, for example we use `Ascl2` instead of `ASCL2`. Verify that the stem cell markers are indeed expressed in the stem cells, as assigned by Haber et al.

Now we're going to leave the PMBC3K Seurat tutorial. For the rest of the assignment, we will make extensive use of gene sets and the `AddModuleScore` function. Google this function to find the documentation, and describe in your own words what it does.

Example usage of this function is as follows:

```
epithelial_gene_set <- c("Epcam", "Krt8")
data <- AddModuleScore(data, features = epithelial_gene_set, name = "epithelial_gene_set")
FeaturePlot(data, features = "epithelial_gene_set1", reduction = "umap")
```

Using the function is in principle straightforward: you pass in your Seurat object `data`, your gene set into the `features` parameter, and give it a name for storing it in the Seurat object. In `FeaturePlot` you then reference this name, as if it were a normal gene. One gotcha is that you need to append "1" to the name in `FeaturePlot`. This is because `AddModuleScore` can also take multiple gene sets at once, and it will append a number to the name for each gene set. (In the example, we only passed one gene set, but we could also have used `features = list(epithelial_gene_set, another_gene_set)`.)

Now use `AddModuleScore` with the stem cell genes you found above, and include this in your report.

The benefit of using gene sets over single genes, is that it is less noisy. If sequencing failed to detect any `Lgr5` transcripts in a cell, either because of technical issues or because the cell was not expressing Lgr5 at that moment, the cell can still be detected as a stem cell if it expresses enough other stem cell markers.

At https://www.gsea-msigdb.org/gsea/msigdb/index.jsp , there are large collections of gene sets, for both humans and mice. On that page, click on "M2" ("curated gene sets", part of "Mouse Collections"). This brings you to a table (which also includes other mouse collections). The curated gene sets are subdivided into various collections. Find the "CP: Canonical pathways" collection, and click on the "browse 1730 gene sets" (number might be updated in the future). This brings you a page with a long list of gene sets.

Use Ctrl + F to search for the gene set "REACTOME_SIGNALING_BY_NOTCH". Click through to the page of this gene set, and download the gene set as "TSV metadata". Open the file in a text editor or Excel, and copy the gene names to a vector in R. You can do this as follows:

```
library(stringr)  # This library is necessary for the str_split function


reactome_signaling_by_notch <- str_split("Kat2b,Akt1,( ... left out ... ),Uba52rt,Uba52", ",")
```

(Even better would be to write a function that reads in the file and extracts the gene names, but this is a bit more advanced, and not necessary for this assignment.)

Then use `AddModuleScore` with this gene set, and make a `FeaturePlot` like before. The color scaling can be a bit unfortunate, since color scaling is optimized for genes. By adjusting the scaling a bit, you can obtain nice results, like in this example:

```
FeaturePlot(data, features = "reactome_signaling_by_notch1", reduction = "umap", min.cutoff = 0.1, max.cutoff = 0.4)
```

Include this plot in your report. There's also the function `VlnPlot` to visualize gene sets, which uses the cell types we've assigned earlier. Example usage of this function, if you just want to compare Paneth and stem cells, is as follows:

```
VlnPlot(data, features = c("reactome_signaling_by_notch1", "reactome_signaling_by_wnt1"), assay = "RNA", idents = c("Paneth", "Stem"))
```

Here, we also included a gene set for Wnt signaling, which you can find on the same page as the Notch signaling gene set. You don't need to include this plot in your report.

As you can see, gene sets can be a powerful tool to analyze single-cell RNA data. They can be used to detect cell types, but also to detect biological processes, such as Notch signaling. You can also search for other processes. In the rest of the assignment, we will focus on metabolic processes, so depending on your choices, you will use gene sets like `REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE` or `WP_ELECTRON_TRANSPORT_CHAIN`. Sometimes, you need to try a few different search keywords to find the gene set you're looking for.


## A little research project

To show the power of using gene sets, we (or rather: you) are going to replicate some of the conclusions of another paper, which was published in the same year as the Haber et al. paper. This paper, about research by the UMC Utrecht, is called "Interplay between metabolic identities in the intestinal crypt supports stem cell function" ( https://www.nature.com/articles/nature21673 ). Read the abstract, list all the conclusions, and see which ones you can investigate with the Haber et al. dataset, and do so for a few of them.

To pass the assignment, you just need to find supporting evidence (one gene set each) for at least two of the conclusions of the paper. You can also try to find evidence for more conclusions, but this is not necessary. You can use gene sets from any source, but in principle using only gene sets on the page you went to before should be enough. You can also compare multiple gene sets, to check whether they give similar results.

And one final question: you'll probably notice that the differences between cell types are not *that* large. Why could this be?
