# Single-Cell Transcriptomics Impromptu Project

This repository contains code and resources to replicate a single-cell metatranscriptomics project based on data from the publication "Microbial single-cell RNA sequencing by split-pool barcoding" by [Kuchina et al](https://www.science.org/doi/full/10.1126/science.aba5257).

## Project Overview

This project aims to analyze bacterial cell gene expression using the microbial split-pool ligation transcriptomics (microSPLiT) method. The data used in this project were obtained from the Gene Expression Omnibus (GEO) public repository under accession number [GSE151940](https://geo.metadataplus.biothings.io/geo/query/acc.cgi?acc=GSE151940). This particular repository centers around comparing the gene expresion of _Escherichia coli_ and _Bacillus subtilis_.

## Dependencies

Ensure you have the following R packages installed:

- `edgeR`: for differential expression analysis.
- `ggplot2` and `ggrepel` for visualization of results.

## Files and Directories

- **sc_transc.R**: R script for running the analysis
- **results**: Directory for storing analysis results.
  - `top_degs.csv`: Results of differential expression analysis.
- **figures**: Directory for storing generated figures.
  - `ma_plot.png`: MA plot visualizing differential expression results.
  - `volcano_plot.png`: Volcano plot visualizing differential expression results.

## Reproducing the Analysis

To replicate the analysis, download the processed data from the [supplementary files](https://geo.metadataplus.biothings.io/geo/download/?acc=GSE151940&format=file). There are 3 different experiments, I went for the M11:
  - `GSM4594094_M11_dcm.csv`: Processed read counts data.
  - `GSM4594094_M11_genes.csv.csv`: Gene names.
  - `GSM4594094_M11_barcodes.csv`: Sample annotations (_Escherichia coli_ vs _Bacillus subtilis_).

![volcano_plot](https://github.com/manuelgug/Single-Cell_Transcriptomics/blob/main/figures/volcano_plot.png)
_Figure 1. Volcano plot of top differencially expressed genes. logFC cutoff was set to 0.1 just to see some colors, so don't take this too seriousy..._
