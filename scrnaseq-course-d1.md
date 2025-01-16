---
title: "Single Cell RNAseq Data Processing & Analysing Course"
author: "Sumukh Deshpande"
#full-width: true
---

---
# Day 1

---
---

We will cover:

- Introduction to course.
- Basic Unix.
- Introduction to HAWK.
- Basic Unix Continued.
- Fetch data with nf-core/fetchngs pipeline.

---

<br>

---
## Course Introduction
---
---

You should have one of three questions:

- Have you found a nice bulk-RNAseq dataset from a paper and want to download and use it for your research?
- Have you recently generated bulk-RNAseq data and need to process and analyse it?
- You are planning a bulk-RNAseq experiment and would like to know how to process the data when the time comes?

---

<br>

---
## What you will learn on this course:
---
---

- Downloading scrnaseq datasets and genome files.
- Setting up directory for the data processing and analysis.
- How to run Nextflow pipeline on HPCs.
- How to process 10x cellranger pipeline.
- How to process and analyse scrnaseq count data in R. 
- How to visualise the analysed data and perform additional downstream analysis.

---

<br>

---
## Introduction to Single-Cell RNA-seq
---
---

- RNA-seq allows profiling the transcripts in a sample in an efficient and cost-effective way.
- Part of its success is due to the fact that RNA-seq allows for an unbiased sampling of all transcripts in a sample, rather than being limited to a pre-determined set of transcripts (as in microarrays or RT-qPCR).
- Typically, RNA-seq has been used in samples composed of a mixture of cells, referred to as **bulk RNA-seq**, and has many applications.
- For example, it can be used to characterise expression signatures between tissues in healthy/diseased, wild-type/mutant or control/treated samples.
- However, with bulk RNA-seq we can only estimate the **average expression level** for each gene across a population of cells, without regard for the heterogeneity in gene expression across individual cells of that sample.

<img src="/assets/img/singlecell-v.-bulk-image.jpg" alt="singlecell" width="1200"/> 


- Unlike with the bulk approach, with scRNA-seq we can estimate a **distribution of expression levels** for each gene across a population of cells.

- This allows us to answer new biological questions where **cell-specific changes in the transcriptome** are important. For example discovering new or rare cell types, identifying differential cell composition between healthy/diseased tissues or understanding cell differentiation during development.

---

<br>

---
## Sample Preparation Protocols
---
---

Broadly speaking, a typical scRNA-seq protocol consists of the following steps (illustrated in the figure below):

- Tissue dissection and cell dissociating to obtain a suspension of cells.

- Optionally cells may be selected 

- Capture single cells into individual reaction containers (e.g. wells or oil droplets).

- Extracting the RNA from each cell.

- Reverse-transcribing the RNA to more stable cDNA.

- Amplifying the cDNA (either by in vitro transcription or by PCR).

- Preparing the sequencing library with adequate molecular adapters.

- Sequencing, usually with paired-end Illumina protocols.

- Processing the raw data to obtain a count matrix of genes-by-cells

- Carrying several downstream analysis (the focus of this course).



---

<br>

---
## Basic Unix: Learning Objectives
---
---

- Learn the concept of using the command line.
- Learn how to navigate and manipulate files and data.
- Learn how to run and manage programs.

---

<br>

---
## Basic Unix: Brief History
---
---

- UNIX is a suite of programs that make up an operating system (like Windows and Mac).
- First developed in 1960's and has been in constant development ever since.
- It's a stable, multi-use, multi-tasking system for servers, desktops, and laptops.
- UNIX systems also have a grahpical user interface (GUI), providing an easy to use Windows-like icon-based environment.
- Linux is a clone of UNIX (they're the same thing). UNIX is a commercial product, whereas Linux is open-source (you can download, install, and use it).
- We will be using Linux on this course.

---

<br>

---
## Basic Unix: Graphical User Interface (GUI) & Command Line/Shell
---
---

- All Windows/Mac/Linux PC's use a GUI to allow users to easily navigate and use the PC.
- The GUI is what allows us to point and click on things, which in turn opens the respective programs etc.
- These operations can also be performed using the command line through the use of a shell.



This is a shell. We can use this to type commands etc.

---

<br>

---
## Introduction to HAWK
---
---

- HAWK is Cardiff and Bangor universities High Performance Compute (HPC) system.
- Swansea and Aberystwyth universities use SunBird - both systems are run by Supercomputing Wales (SCW).
- Both HPC's are the same, but set-up slightly different - we will be working with HAWK.
- If we need to run analyses/softwares/code that requires a lot of computing power, we need to use HPC's.
- The way we interact with HAWK is through the command line via Bash.
- Think of HAWK as a computer that's located in the cloud. We can connect to it via out Unix shells - Terminal, MobaXterm, and FileZilla.

---

<br>

---
