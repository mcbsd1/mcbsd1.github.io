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
  - Cell-type annotation using CellMarker.


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

- Create a script called "**NSCLC_dataset.Rmd**" by following the instructions as instructed on **Day-2**.

- Set the working directory to the bin folder where the script named "NSCLC_dataset.Rmd" will be stored.
Setting up working directory:

```
setwd(/path/to/your/bin)
```

Load libraries:

```
library(Seurat)
library(DoubletFinder)
```

Load all the datasets:

```
# Loading P1 Tumor processed data

data <- read.table("../input/GSM3304007_P1_Tumor_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p1tumor <- CreateSeuratObject(counts = data1, project = "GSM3304007")

# Loading P1 Normal processed data

data <- read.table("../input/GSM3304008_P1_Normal_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p1normal <- CreateSeuratObject(counts = data1, project = "GSM3304008")

# Loading P2 Tumor processed data

data <- read.table("../input/GSM3304009_P2_Tumor_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p2tumor <- CreateSeuratObject(counts = data1, project = "GSM3304009")

# Loading P2 Normal processed data

data <- read.table("../input/GSM3304010_P2_Normal_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p2normal <- CreateSeuratObject(counts = data1, project = "GSM3304010")

# Loading P3 Tumor processed data

data <- read.table("../input/GSM3304011_P3_Tumor_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p3tumor <- CreateSeuratObject(counts = data1, project = "GSM3304011")

# Loading P3 Normal processed data

data <- read.table("../input/GSM3304012_P3_Normal_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p3normal <- CreateSeuratObject(counts = data1, project = "GSM3304012")

# Loading P4 Tumor processed data

data <- read.table("../input/GSM3304013_P4_Tumor_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p4tumor <- CreateSeuratObject(counts = data1, project = "GSM3304013")

# Loading P4 Normal processed data

data <- read.table("../input/GSM3304014_P4_Normal_processed_data.txt",row.names = 1,header = T)
data1 <- data[-1]
p4normal <- CreateSeuratObject(counts = data1, project = "GSM3304014")
```

Add percentage of reads that map to mitochondrial genome

```
p1tumor[["percent.mt"]] <- PercentageFeatureSet(p1tumor, pattern = "^MT-")
p1normal[["percent.mt"]] <- PercentageFeatureSet(p1normal, pattern = "^MT-")
p2tumor[["percent.mt"]] <- PercentageFeatureSet(p2tumor, pattern = "^MT-")
p2normal[["percent.mt"]] <- PercentageFeatureSet(p2normal, pattern = "^MT-")
p3tumor[["percent.mt"]] <- PercentageFeatureSet(p3tumor, pattern = "^MT-")
p3normal[["percent.mt"]] <- PercentageFeatureSet(p3normal, pattern = "^MT-")
p4tumor[["percent.mt"]] <- PercentageFeatureSet(p4tumor, pattern = "^MT-")
p4normal[["percent.mt"]] <- PercentageFeatureSet(p4normal, pattern = "^MT-")
```

Visualise QC metrics as violin plots

```
VlnPlot(p1tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p1normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p2tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p2normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p3tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p3normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p4tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(p4normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

Now, filter cells. The first one has been done for you.

```
p1tumorFilteredObj <- subset(p1tumor, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &  percent.mt < 10)
```

---
#### Normalization, Scaling and Dimensional Reduction
---

Once the samples are filtered, then run the following steps for each sample in the iteration:

```
# Define the sample names
samples <- c("p1tumor", "p1normal", "p2tumor", "p2normal", "p3tumor", "p3normal", "p4tumor", "p4normal")

# Loop through each sample
for (sample in samples) {
  
  # Get the original Seurat object dynamically
  seurat_obj <- get(paste0(sample, "FilteredObj"))
  
  normalized_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  normalized2_obj <- FindVariableFeatures(normalized_obj, selection.method = "vst", nfeatures = 2000)
  
  all_genes <- rownames(normalized2_obj)
  scaled_obj <- ScaleData(object = normalized2_obj, features = all_genes)
  
  scaled2_obj <- RunPCA(object = scaled_obj, features = VariableFeatures(object = scaled_obj), 
                        ndims.print = 1:5, nfeatures.print = 5)
  
  ElbowPlot(scaled2_obj, ndims = 50)
  
  neighbors_obj <- FindNeighbors(scaled2_obj, reduction = "pca", dims = 1:25, k.param = 30)
  
  clusters_obj <- FindClusters(object = neighbors_obj)
  
  umap_obj <- RunUMAP(clusters_obj, reduction = "pca", dims = 1:25, n.neighbors = 30L, min.dist = 0.3, seed.use = 123456L)
  
  # Store results dynamically in the environment
  assign(paste0(sample, "NormalisedObj"), normalized_obj)
  assign(paste0(sample, "Normalised2Obj"), normalized2_obj)
  assign(paste0(sample, "ScaledObj"), scaled_obj)
  assign(paste0(sample, "Scaled2Obj"), scaled2_obj)
  assign(paste0(sample, "FindNeighbors"), neighbors_obj)
  assign(paste0(sample, "FindClusters"), clusters_obj)
  assign(paste0(sample, "UMAP"), umap_obj)
}
```

The results for each sample will be stored in its corresponding sample name + UMAP. For example, if you would like to access the object for p3tumor, the results is stored in p3tumorUMAP. For p4normal, the object is p4normalUMAP, and so on.

---
### Doublet removal
---

---
#### Doublet identification
---

- The next step is to identify doublets in the datasets. For this we will be using `DoubletFinder` tool in R.
- Follow the steps as mentioned to identify and label the cells as `Singlet` or `Doublet`.

The first sample has been done for you:

```
sweep.res.list <- paramSweep(p1tumorFilteredObj, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])

annotations <- p1tumorFilteredObj@active.ident
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*length(p1tumorFilteredObj@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seuratDoublet <- doubletFinder(p1tumorFilteredObj, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
pANN_String <- paste("pANN_0.25_",pk,"_",list(nExp_poi),sep="")

p1tumorFilteredDoublet <- doubletFinder(seuratDoublet, PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = pANN_String)
```

```
# Use this format: DF.classifications_pN_pk_nExp_poi

classificationsCol <- p1tumorFilteredDoublet@meta.data["DF.classifications_0.25_0.24_134"]

p1tumorFilteredDoublet@meta.data[,"DF_hi.lo"] <- classificationsCol[1]

```

Follow the steps to generate doublet marked objects for rest of the samples.

---
#### Singlet selection and doublet removal
---

- In the next step, we use `subset` function in R to select Singlet and filter Doublet cells from the dataset.
- This is iterated for all the samples.

```
p1tumorSinglets <- subset(p1tumorFilteredDoublet, subset = DF_hi.lo == "Singlet")
p1normalSinglets <- subset(p1normalFilteredDoublet, subset = DF_hi.lo == "Singlet")
p2tumorSinglets <- subset(p2tumorFilteredDoublet, subset = DF_hi.lo == "Singlet")
p2normalSinglets <- subset(p2normalFilteredDoublet, subset = DF_hi.lo == "Singlet")
p3tumorSinglets <- subset(p3tumorFilteredDoublet, subset = DF_hi.lo == "Singlet")
p3normalSinglets <- subset(p3normalFilteredDoublet, subset = DF_hi.lo == "Singlet")
p4tumorSinglets <- subset(p4tumorFilteredDoublet, subset = DF_hi.lo == "Singlet")
p4normalSinglets <- subset(p4normalFilteredDoublet, subset = DF_hi.lo == "Singlet")
```

---
### Normalization and finding variable features
---

- In the next step, simply create a list called `sample_list` which stores all the samples.

```
sample_list <- list()
sample_list[["p1tumorSinglets"]] <- p1tumor
sample_list[["p1normalSinglets"]] <- p1normal
sample_list[["p2tumorSinglets"]] <- p2tumor
sample_list[["p2normalSinglets"]] <- p2normal
sample_list[["p3tumorSinglets"]] <- p3tumor
sample_list[["p3normalSinglets"]] <- p3normal
sample_list[["p4tumorSinglets"]] <- p4tumor
sample_list[["p4normalSinglets"]] <- p4normal
```

- Once the samples are added to the list, run `NormalizeData` and `FindVariableFeatures` functions to perform normalization and finding variable features using VST and `nfeatures=2000`.
- Interate this step over all the samples in the `sample_list`. 

```
for (i in 1:length(sample_list)) {
  sample_list[[i]] <- NormalizeData(sample_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sample_list[[i]] <- FindVariableFeatures(sample_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = T)
}
```

---
### Perform Integration of samples
---

- We will not integrate data with 8 samples (2 conditions for each Patient).
- First step is to find integration anchors.

#### Step 1: Finding Integration Anchors

```
sample_anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:30)
```

- This step identifies anchors (shared biological structures) across multiple datasets in sample_list to help correct batch effects.
- It first selects highly variable genes (default = 2000 features).
- It then runs canonical correlation analysis (CCA) to find common biological signals across datasets.
- The parameter dims = 1:30 means that the first 30 principal components (PCs) are used to define the anchors.

#### Step 2: Integrating the Data

```
allsamples_seurat <- IntegrateData(anchorset = sample_anchors, dims = 1:30)
```

- This step applies batch correction using the integration anchors found in Step 1.
- It merges all datasets into a single Seurat object (allsamples_seurat), where only the 2000 most variable features are stored in the "integrated" assay.
- The batch-corrected data is now stored in the "integrated" assay.

#### Step 3: Switching Back to RNA Assay

```
DefaultAssay(allsamples_seurat) <- "RNA"
```

- Changes the active assay from "integrated" to "RNA".
- This is necessary for normalization, scaling, and differential expression analysis because "integrated" data should not be used for these steps.

Why It’s Needed:

- The "integrated" assay is only useful for dimensionality reduction (PCA, UMAP), clustering, and visualization.
- Gene expression comparisons (e.g., `FindMarkers()`) should be done on the "RNA" assay to retain true biological signal.

#### Step 4: Merging Layers

```
allsamples_seurat1 <- JoinLayers(allsamples_seurat)
```

- This step merges different data layers within the Seurat object to ensure all layers are accessible.
- Specifically, it merges:
  - Raw expression data (RNA)
  - Integrated (batch-corrected) data
  - Any other assays or metadata layers
- This is useful for downstream analysis where you might need access to both raw and batch-corrected data.

---
### Perform Data Scaling
---

- We will be using `all.genes` to capture all genes for data scaling.

```
all.genes <- rownames(allsamples_seurat1)
allsamplesScaleData <- ScaleData(allsamples_seurat1, features = all.genes)
```

---
### Perform linear dimensional reduction
---

- Next we perform PCA on the scaled data. 
- We use variable features by default, in which case only the previously determined variable features are used as input.

```
allsamplesPCA <- RunPCA(object = allsamplesScaleData, features = VariableFeatures(object = allsamplesScaleData), do.print = TRUE, ndims.print = 1:5, nfeatures.print = 5)
```

#### Plot ElbowPlot to calculate number of PCs

- Next, we plot Elbow plot by running `ElbowPlot` function to determine number of PCs which define highest variance in the datasets.
- We create this plot with `ndims = 50`

```
ElbowPlot(allsamplesPCA, ndims = 50)
```

<img src="/assets/img/ElbowPlot_cancerData.png" alt="gui1" width="1200"/>

In this tutorial, I have chosen PC 20. Beyond this number, the variance seems to be stabilised and should be good enough to perform clustering.

---
### Clustering
---

- Next, we perform clustering by running `FindNeighbors()` and `FindClusters()` functions using the dimensionality determined by plotting ElbowPlot().
- We take first 20 PCs as an input in this case.

```
allsamplesFN <- FindNeighbors(allsamplesPCA, reduction = "pca", dims = 1:20)
allsamplesFC <- FindClusters(object = allsamplesFN)
```

---
### Non-linear dimensional reduction
---

- In the next step, we run the `RunUMAP()` function to group the cells together constructing a UMAP plot.
- Again, we use first 20 PCs.

```
allsamplesUMAP <- RunUMAP(allsamplesFC, reduction = "pca", dims = 1:20, verbose = T)
```

---
### Generate UMAP plot
---

- Before we generate the UMAP plot, it is important to change the default assay layer to integrated.

```
DefaultAssay(allsamplesUMAP) <- "integrated"
DimPlot(allsamplesUMAP, reduction = "umap")
```

<img src="/assets/img/UMAP_cancerData.png" alt="gui1" width="1200"/>

#### UMAP plot by sample name

- To plot the UMAP by sampleName information stored in `meta.data`, use the following code:

```
DimPlot(allsamplesUMAP, reduction = "umap", group.by = "sampleName")
```

<img src="/assets/img/UMAP_sampleName.png" alt="gui1" width="1200"/>


#### UMAP plot by sample type (Tumor/Normal)

- To plot the UMAP by sampleType information stored in `meta.data`, use the following code:

```
DimPlot(allsamplesUMAP, reduction = "umap", group.by = "sampleType")
```

<img src="/assets/img/UMAP_sampleType.png" alt="gui1" width="1200"/>

#### UMAP plot splitted by sample type (Tumor/Normal)

- To plot the UMAP split-by sampleType information stored in `meta.data`, use the following code:

```
DimPlot(allsamplesUMAP, reduction = "umap", split.by = "sampleType")
```

<img src="/assets/img/UMAP_splitbySampleType.png" alt="gui1" width="1200"/>

- Now, turn on the labels for each cluster.
- **How do you turn-on the labels ?**
- **Which clusters are expressed in Tumor and Normal cells ?**


<br>

---
## Differential Gene Expression (DGE)
---

- In Seurat, Differential Expression (DE) analysis is used to identify genes that are differentially expressed between different cell groups (clusters, conditions, or sample types).
- Seurat provides several methods for performing DGE, primarily through the `FindMarkers()` and `FindAllMarkers()` functions.
- For Example: If you would like to find all the markers for clusters 2, then you would run the following command:

```
cluster2.markers <- FindMarkers(allsamplesUMAP, ident.1 = 2)
```

Now, if you want to find all markers distinguishing cluster 5 from clusters 0 and 3, run the following command:

```
cluster5.markers <- FindMarkers(allsamplesUMAP, ident.1 = 5, ident.2 = c(0, 3))
```

List the top 5 markers which are expressed significantly in cluster 2.

```
head(cluster2.markers, n = 5)
```

If you want to find markers for every cluster compared to all remaining cells, then run the following command:

```
allsamples.markers <- FindAllMarkers(allsamplesUMAP)
```

To plot the expression of genes expressed in cluster 5, run the following command:

```
VlnPlot(allsamplesUMAP, features = c("TYROBP", "KLRF1", "SPON2", "FCGR3A"), ncol = 2)
```

<img src="/assets/img/ViolinPlot_clus5_markers.png" alt="gui1" width="1200"/>

Now, to plot the expression of cluster 2 and 5 genes as feature plots, run the following command:

```
FeaturePlot(allsamplesUMAP, features = c("TYROBP", "KLRF1", "SPON2", "FCGR3A", "FGFBP2", "FCN1", "LST1", "AIF1", "C5AR1"))
```

<img src="/assets/img/FeaturePlot_cancerData.png" alt="gui1" width="1200"/>


- To plot the Feature plot by sampleType i.e Tumor/Normal for *JUN1* and *MALAT1* genes which are upregulated in the intratumoral epithelial cell, run the following command:

```
FeaturePlot(allsamplesUMAP, features = c("JUN", "MALAT1"), split.by = "sampleType")
```

<img src="/assets/img/FeaturePlot_JUN1_MALAT1.png" alt="gui1" width="1200"/>

<br>

---
## Cell Type Identification using CellMarker
---

- CellMarker is a manually curated database human and mouse species.
- Its manually curated over 100,000 published papers, 4124 entries including the cell marker information.
- **13,605** cell markers of **467** cell types in **158** human tissues/sub-tissues and **9148** cell makers of **389** cell types in **81** mouse tissues/sub-tissues were collected and deposited in CellMarker.
- **Paper Reference**: [CellMarker: a manually curated resource of cell markers in human and mouse.](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky900/5115823)

---
---
- For this dataset, we will be downloading the "Human cell markers" from CellMarker database.
- **Human Cell markers:** [download](http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx) or visit the download page: [download page](http://117.50.127.228/CellMarker/CellMarker_download.html)
- Other cell markers are as follows:
  - Mouse cell markers: [download]([http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt](http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Mouse.xlsx))
  - Single cell markers: [download]([http://xteam.xbio.top/CellMarker/download/Single_cell_markers.txt](http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx))
  - All cell markers: [download]([http://xteam.xbio.top/CellMarker/download/all_cell_markers.txt](http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_All.xlsx))


