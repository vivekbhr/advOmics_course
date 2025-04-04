# Identifying mutagenic processes in bone marrow transplant samples

## System set up

-**Machine** No big RAM requirements, any normal laptop should do.

-**OS** MAC/Windows/Linux. Tested on macOS Sequioa 15.3.2

-**Language** R (mainly for MutationalPatterns). Tested with R version 4.4.2

-**Slack channel name** project_mutationalPatterns 

## Installing packages in R

```
install.package('dplyr')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MutationalPatterns")

```

## Introduction

Different mutagenic processes leave distinct mutational patterns on the DNA. Identifying these patterns can teach us a lot about the history of the cell, and for example help determine the etiology of a cancer.
They are also helpful to see if specific drugs are mutagenic, provided that these drugs leave a specific patterns. In this project you are going to take a look at one of these drugs.

## Paper & data
You will be working on data from the paper "Antiviral treatment causes a unique mutational signature in cancers of transplantation recipients" by De Kanter et al.
This paper investigates specific mutational patterns found in the HSPCs of patients who received a hematopoetic stem cell transplantation as cancer treatment. They perform whole genome sequencing both on the donor hematopoietic stem cell and the cell after transplantation.
They provide all mutations identified over 51 HSPCs in Supplementary table 3, which can be downloaded here: https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(21)00337-4#mmc1

## Objectives
The goal is to reproduce the main finding of this paper, which is the unique mutational signature they identify in this paper. You will use the R package MutationalPatterns for this.
Mutational patterns expects vcfs as input data, and the mutational calls are provided in an excel file by the paper. So first, you will have to produce the vcfs files.
A standardized format like vcf is really handy, but it also means it expects very particular input. You can learn more about the structure of vcfs here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format
The easiest way is to define a vcf header and a body and then write it to a file together. To help you along, here is some code that makes a header that works:

```
library(dplyr)
## sample = sample identifier
  vcf_header <- c(
      "##fileformat=VCFv4.2",
      "##source=CustomScript",
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
      paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample)
    )



```

Some of these columns a vcf expects to have is not present in the excel file, but you can assume that the calls provided here have passed quality control and the filter.
Just set these to pass and/or 100. You can add extra INFO if you want if you want to include the VAF for example.

```
 vcf_body <- sample_df %>%
      mutate(ID = ".", QUAL = ".", FILTER = "PASS", INFO = "DP=100", FORMAT = "GT", GENOTYPE = "0/1") %>%
      select(Chromosome, Position, ID, Reference, Alternative, QUAL, FILTER, INFO, FORMAT, GENOTYPE)
```


Then make sure you generate a vcf body that contains all necessary info and write the head and body per sample to a .vcf file and the sample info is ready for use in MutationalPatterns!
As we discussed in the lecture, there are two main ways of identifying signatures: a refit and a de novo extraction. We are interested in a novel mutational signature here, so we will have to
do de novo extraction, but refitting on the COSMIC database can still be really helpful to see which known processes are active. You can also use the COSMIC database to identify signatures from your de novo extraction.
A big challenge in this type of research is really proving that you have found a new process, so focus on what evidence you can give that this signature is truly novel and you are not just misidentifying a known process.
Can you find any connection between contribution of signatures and total mutation load? What could that tell you about the exposure?

## Expected output

- The mutational signature to be identified and evidence of the exposure in the samples, in the form of 96 mutation plot and signature contribution plots
- Evidence of novelty, in the form of (lack of) cosine similarity to known signature, in a cosine similarity plot. You can also show an increase in reconstruction accuracy when the new signature is included vs a refit on known signatures.
- Evidence of (lack of) correlation between mutation load and signatures

## Tips

Be aware that both de novo extraction and refit will always give you a result, even if you are fitting on nonsense. For the de novo extraction the number of signatures you are extracting is very important!
There are several metrics the mutationalPatterns package report to help you choose an appropriate number of signatures to extract, motivate your choice!

Good luck!
