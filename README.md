# Proteomics Data Reanalysis: Z-culisetae Experimental Groups

## Overview

This repository contains the R script used for analyzing Label-Free Quantification (LFQ) mass spectrometry proteomics data. The goal is to identify significantly differentially expressed proteins across four experimental conditions (Control, Dark, Light, Ecdy) using a standard, data-driven workflow based on ANOVA and pairwise comparisons.

The script incorporates essential steps for handling typical proteomics data challenges, including **missing value filtering**, **normalization**, and **imputation**.

---

## Analysis Workflow

The script executes the following steps in sequence:

1.  **Data Loading & Pre-processing**
    * Load the primary data file.
    * Remove the E11 sample due to its status as an outlier with high missingness.
    * Filters proteins that do not have at least 3 valid LFQ values in at least one experimental group.
    * Applies Median Normalization to correct for systematic loading errors.
    * Performs Left-Censored Imputation (down-shifted narrow normal distribution) to replace remaining missing values.

2.  **Quality Control & Visualization**
    * Performs Principal Component Analysis (PCA) on the pre-processed, imputed data.
    * Generates a PCA plot to confirm successful clustering of biological replicates and separation of experimental groups.

3.  **Statistical Testing**
    * **Global Test (ANOVA):** Performs a one-way ANOVA across the four conditions (C, D, L, E) for every protein to find any significant difference.
    * **Multiple Testing Correction:** Applies the Benjamini-Hochberg (BH) procedure to the p-values to control the False Discovery Rate (FDR).
    * **Post-Hoc Analysis:** Calculates Log2 Fold Change (Log2FC) for the comparison of Light (L) vs. Control (C).

4.  **Final Visualization**
    * Generates a Heatmap of all proteins found to be significant by the ANOVA (FDR < 0.05).
    * Generates a Volcano Plot for the Light vs. Control comparison, highlighting proteins that are both statistically significant (FDR < 0.05) and have a strong magnitude of change ($\mid Log2FC \mid > 1$).

---

## How to Run the Script

### Prerequisites

You need the following software and R packages installed:

1.  **R** (latest version recommended)
2.  **RStudio** (recommended IDE)

### Package Installation (Run this once)

Open RStudio and run the following in the Console:

```r
install.packages(c("tidyverse", "pheatmap", "ggrepel"))
```
### AI Acknowledgement
This README was edited by AI for formatting and clarity.