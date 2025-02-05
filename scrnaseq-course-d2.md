---
title: "Single Cell RNAseq Data Processing & Analysing Course"
author: "Sumukh Deshpande"
#full-width: true
---

---
# Day 2 - Processing scRNAseq reads generated from Nextflow pipeline

---

## Recap from Day 1
---
---

- Yesterday we covered:
  - How to use Unix.
  - How to get onto HAWK.
  - How to fetch the Github repository and download the resources.
  - How to execute the scrnaseq nextflow pipeline on Hawk.

---
## Outputs from the nextflow pipeline
---
---

- The output from the nextflow execution run is stored under `results/cellranger_count_1k_mouse_kidney_CNIK_3pv3`
- Login to your Hawk account and access the folder:

```
cd /scratch/c.123456/scrnaseq-course/scrnaseq/results

cd cellranger_count_1k_mouse_kidney_CNIK_3pv3

```

---
## Downloading data for analysis in Seurat
---
---

- Once the counts are generated using nextflow pipeline, download the data onto your personal computer.
- For this, you will need to use Filezilla to copy the files from remote onto your local machine.
- Follow the steps below:

1. Type in the Host: hawklogin.cf.ac.uk
2. Username: <hawk username>
3. Password: <hawk password>
4. Port: 22

- Click on `Quickconnect`
- Navigate to the folder `scrnaseq-course/scrnaseq/results`
- Create a folder on your local machine, called `scrnaseq-nextflow`.
- Create 4 folders within `scrnaseq-nextflow` folder:
  - bin
  - input
  - output
  - resources
- Now, Drag and Drop the files from the `results` folder in the remote server to your `input` folder within `scrnaseq-nextflow` folder on your local machine.
- Once the files are copied, these should be visible on your local machine under `input` folder.

---
<img src="/assets/img/filezilla-login.png" alt="gui1" width="1200"/>


---
## Processing scrnaseq data (1k_mouse_kidney_CNIK_3pv3) in Seurat 
---
---

- Once the counts have been downloaded into the `input` folder on your local machine, Start **RStudio** on Mac/Windows.

---
<img src="/assets/img/rstudio.png" alt="gui1" width="1200"/>

---
- Click on `New File` --> `R Markdown...`
<img src="/assets/img/rstudio-newfile.png" alt="gui1" width="800"/>

---
- This is will open an `Untiltled` file on the top-left hand side.

---
<img src="/assets/img/rstudio-untilted-rmarkdown.png" alt="gui1" width="800"/>

---
- Delete lines from 7 until 29 except first 5-6 lines.

---
<img src="/assets/img/rmarkdown-2.png" alt="gui1" width="800"/>

<br>

---
#### Load libraries
---

- Copy the line in your markdown document:

```
# load libraries
library(Seurat)
```

- Press the green arrow key on line 7 which says **Run Current Chunk** to load the Seurat library.

<img src="/assets/img/rmarkdown-1.png" alt="gui1" width="800"/>


---
#### Setup the Seurat object
---

- We start by reading in the data. The `Read10X()` function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix.
- The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).

```
# Load the 1k mouse kidney dataset
dat2 <- Read10X(data.dir = "input/filtered_feature_bc_matrix/")
```

- We next use the count matrix to create a Seurat object. 

```
# Initialize the Seurat object with the raw (non-normalized data).
dat2Obj <- CreateSeuratObject(counts = dat2, project = "1KmouseKidney", min.cells = 3, min.features = 200)

dat2Obj
## An object of class Seurat 
## 18321 features across 1375 samples within 1 assay 
## Active assay: RNA (18321 features, 0 variable features)
##  1 layer present: counts
```

---
#### QC and filtering cells for further analysis
---

- Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria.
- These includes:
  - The number of unique genes detected in each cell.
    - Low-quality cells or empty droplets will often have very few genes
    - Cell doublets or multiplets may exhibit an aberrantly high gene count
  - The percentage of reads that map to the mitochondrial genome.
    - Low-quality / dying cells often exhibit extensive mitochondrial contamination
    - We calculate mitochondrial QC metrics with the `PercentageFeatureSet()` function
    - We use the set of all genes starting with `MT-` or `mt-` as a set of mitochondrial genes

```
# The [[ operator can add columns to object metadata.
dat2Obj[["percent.mt"]] <- PercentageFeatureSet(dat2Obj, pattern = "^mt-")
```

- To visualize QC metrics, we use `VlnPlot()` function.

```
# Visualize QC metrics as a violin plot
VlnPlot(dat2Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

<img src="/assets/img/Violin_plot_1k_data.png" alt="gui1" width="1200"/>

- Next, you can plot a FeatureScatter plot for visualising feature-feature relationships.

```
# Plot FeatureScatter plot
plot1 <- FeatureScatter(dat2Obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat2Obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

<img src="/assets/img/FeatureScatter.png" alt="gui1" width="1200"/>

<br>

- Now, we will filter cells that have unique feature counts over 9000 or less than 200
- We filter cells that have >5% mitochondrial counts

```
dat2FilteredObj <- subset(dat2Obj, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 &  percent.mt < 5)

dat2FilteredObj
# An object of class Seurat 
# 18321 features across 1235 samples within 1 assay 
# Active assay: RNA (18321 features, 0 variable features)
#  1 layer present: counts
```

Lets break it down:
- `nFeature_RNA > 200`: Keeps cells that have more than 200 detected genes (features). This removes low-quality cells or empty droplets.
- `Feature_RNA < 9000`: Removes cells with too many detected genes, which might indicate doublets (two cells captured in one droplet).
- `percent.mt < 5`: Filters out cells where mitochondrial gene percentage is greater than 5%. High mitochondrial content often indicates damaged or dying cells.

- With this command, the number of cells reduces to 1235.


---
#### Normalizing the data
---

- After removing unwanted cells from the dataset, the next step is to normalize the data.
- By default, we employ a global-scaling normalization method `LogNormalize` that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
- Normalized values are stored in `dat2NormalisedObj[["RNA"]]$data`

```
dat2NormalisedObj <- NormalizeData(object = dat2FilteredObj, normalization.method = "LogNormalize", scale.factor = 10000)
```

- The use of `SCTransform` replaces the need to run `NormalizeData`, `FindVariableFeatures`, or `ScaleData`.


---
#### Identification of highly variable features (feature selection)
---

- We next calculate a subset of features (genes) that exhibit high cell-to-cell variation in the dataset.
- These are the features those are highly expressed in some cells, and lowly expressed in others.
- It has been found [Brennecke et al. Nat Methods.2013](https://www.nature.com/articles/nmeth.2645) that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

```
dat2Normalised2Obj <- FindVariableFeatures(dat2NormalisedObj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dat2Normalised2Obj), 10)
plot3 <- VariableFeaturePlot(dat2Normalised2Obj)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))
```

<img src="/assets/img/VariableGenesPlot.png" alt="gui1" width="1200"/>

---
#### Scaling the data
---

- The `ScaleData` function in Seurat is used to normalize and scale gene expression values across cells.
- Additionally, it can **regress out unwanted sources of variation**, such as sequencing depth (`nUMI`) or mitochondrial gene percentage (`percent.mito`), helping to remove technical artifacts.
- By default, only variable features are scaled.
- But if you observe, in out script, we have used `rownames(dat2Normalised2Obj)` which applies scaling to all genes which is stored in `all.genes` object.
- This is then applied to `ScaleData` function which uses `all.genes` as features for scaling.

```
all.genes <- rownames(dat2Normalised2Obj)
dat2ScaledObj <- ScaleData(object = dat2Normalised2Obj, features = all.genes)
```

<br>

<details>
<summary><b>How can I remove unwanted sources of variation</b></summary>

- To remove unwanted sources of variation from a single-cell dataset, such as "percent.mt", use "vars.to.regress" to include these:

<pre><span style="color:crimson;">dat2ScaledObj <- ScaleData(object = dat2Normalised2Obj, features = all.genes, vars.to.regress = "percent.mt")</span></pre>

</details>
