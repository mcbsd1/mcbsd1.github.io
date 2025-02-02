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
- Create a folder on your local machine, called `scrnaseq-nextflow`
- Now, Drag and Drop the files from the `results` folder in the remote server to your `scrnaseq-nextflow` folder on your lopcal machine.
- Once the files are copied, these should be visible on your local machine under `scrnaseq-nextflow` folder.

---
<img src="/assets/img/filezilla-login.png" alt="gui1" width="1200"/>


---
## Processing scrnaseq data (1k_mouse_kidney_CNIK_3pv3) in Seurat 
---
---

- Once the counts have been downloaded into the `scrnaseq-nextflow` folder on your lopcal machine, start **RStudio** on Mac/Windows.

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
---

- Copy the line in your markdown document:

```{r}
# load libraries

library(Seurat)
```

<img src="/assets/img/rmarkdown-1.png" alt="gui1" width="800"/>



