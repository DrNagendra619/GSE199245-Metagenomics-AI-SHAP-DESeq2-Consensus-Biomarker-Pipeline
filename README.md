# GSE199245-Metagenomics-AI-SHAP-DESeq2-Consensus-Biomarker-Pipeline
Integrative Clinical Metagenomic Profiling of the MET4IO Trial GSE199245 A Consensus Framework Combining Explainable AI XGBoost SHAP and Differential Abundance Statistics DESeq2
# Integrative Clinical Metagenomic Profiling of the MET4-IO Trial (GSE199245)

### A Consensus Framework Combining Explainable AI (XGBoost-SHAP) and Differential Abundance Statistics (DESeq2)

![R](https://img.shields.io/badge/R-4.5.0+-blue.svg)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-green.svg)
![Pipeline](https://img.shields.io/badge/Pipeline-Phyloseq%20|%20XGBoost%20|%20DESeq2-orange.svg)
![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)

## 📌 Overview
This repository contains a robust bioinformatics pipeline for analyzing longitudinal microbiome data from the **MET4-IO clinical trial (GSE199245)**. The pipeline integrates **Machine Learning (XGBoost)** with **Explainable AI (SHAP)** and traditional **Biostatistics (DESeq2)** to identify "Gold Standard" biomarkers associated with treatment response.

Unlike standard workflows, this pipeline features **self-correcting metadata inference** and **numeric sanitization**, making it crash-proof against "dirty" real-world clinical data.

## 🔬 Methodology
The pipeline employs a "Consensus Approach" to biomarker discovery:
1.  **Ecological Diversity:** Evaluates Alpha (Shannon) and Beta (Bray-Curtis PCoA) diversity changes over time.
2.  **Predictive Modeling (AI):** Trains an **XGBoost** classifier to distinguish Baseline vs. Post-Treatment samples.
3.  **Explainability (XAI):** Uses **SHAP (SHapley Additive exPlanations)** to calculate the exact contribution of specific bacterial taxa to the model's predictions.
4.  **Statistical Validation:** Runs **DESeq2** (Wald test) to identify statistically significant differentially abundant taxa ($p_{adj} < 0.05$).
5.  **Consensus:** Intersects the AI-identified features with Statistically significant features to define robust biomarkers.

## ⚙️ Pipeline Steps
1.  **Automated Data Retrieval:** Fetches raw count matrices and metadata directly from NCBI GEO.
2.  **Numeric Sanitization:** Scans and removes non-numeric artifacts (hidden text columns) from count matrices.
3.  **Dynamic Metadata Inference:** Automatically infers clinical groups (`Baseline` vs `Post_Tx`) based on sample date codes in headers, bypassing missing/mismatched metadata files.
4.  **Phyloseq Construction:** Builds the master biological object.
5.  **Machine Learning:** Trains XGBoost with optimized hyperparameters for sparse microbiome data.
6.  **Performance Validation:** Generates ROC Curves and AUC scores.
7.  **Visualization:** Produces publication-ready Boxplots, Heatmaps, and Venn Diagrams.

## 📦 Dependencies
This pipeline is optimized for **R 4.5.0** and **Bioconductor 3.22**.

### Required R Packages
```r
# Core Bioinformatics
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "GenomicRanges", "S4Vectors", "IRanges"))

# Machine Learning & Stats
install.packages(c("xgboost", "SHAPforxgboost", "caret", "pROC"))

# Visualization & Data Manipulation
install.packages(c("ggplot2", "plotly", "ggpubr", "ggVennDiagram", "data.table", "dplyr"))
