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

- Kidney processed sequencing data have been taken from GEO database (GSE131685).
- Kidney specimens were collected fresh, dissected and digested into single cells from organ donors (two males and one female) aged 57–65 years.
- Paper reference: [Liao et al. Sci Data. 2020](https://www.nature.com/articles/s41597-019-0351-8)
