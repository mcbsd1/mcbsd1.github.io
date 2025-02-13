---
title: "Single Cell RNAseq Data Processing & Analysing Course"
author: "Sumukh Deshpande"
#full-width: true
---

---
## Day 3 - Analysing scRNAseq data in Seurat

---

### Recap from Day 2
---
---

- On 25th Feb, we covered:
  - How to load cellranger data in RStudio.
  - How to process scRNAseq data in Seurat.
    - Loading and creating seurat object
    - Generating QC metrics
    - Filtering data
    - Normalising data
    - Identifying highly variable genes
    - Scaling data
    - Performing Principal Component Analysis
    - Finding clusters (Clustering)
    - Visualising clusters using UMAP
    - Other visualisation functions to create plots for identifying differences in gene expression between clusters.


- Today we will cover:
  - Integrating single-cell datasets.
  - Differential Gene Expression analysis and marker selection.
  - Cell-type annotation using SingleR.


---
---
### scRNAseq Data Integration
---
---

- Integration of single-cell sequencing datasets, for example across experimental batches, donors, or conditions, is often an important step in scRNA-seq workflows.
- Integrative analysis can help to match shared cell types and states across datasets, which can boost statistical power, and most importantly, facilitate accurate comparative analysis across datasets.

---
#### Experiment
---

- Cancer NSCLC processed sequencing data have been taken from GEO database (GSE117570).
- 4 early‐stage non‐small cell lung cancer (NSCLC) patients.
- The Cell Ranger Single Cell Software Suite v.2.0.1 was used to perform sample de‐multiplexing, alignment, filtering, and UMI counting.
- Paper reference: [Qianqian et al. Cancer Medicine. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6558497/)
- <u>Input files:</u>
  - Patient 1 - Normal
    - GSM3304008_P1_Normal_processed_data.txt
  - Patient 1 - Tumor
    - GSM3304007_P1_Tumor_processed_data.txt
  - Patient 2 - Normal
    - GSM3304010_P2_Normal_processed_data.txt
  - Patient 2 - Tumor
    - GSM3304009_P2_Tumor_processed_data.txt
  - Patient 3 - Normal
    - GSM3304012_P3_Normal_processed_data.txt
  - Patient 3 - Tumor
    - GSM3304011_P3_Tumor_processed_data.txt
  - Patient 4 - Normal
    - GSM3304014_P4_Normal_processed_data.txt
  - Patient 4 - Tumor
    - GSM3304013_P4_Tumor_processed_data.txt

---
#### Setup and QC
---


Load datasets:

```

```
