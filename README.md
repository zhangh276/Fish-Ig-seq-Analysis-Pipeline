# Fish Ig-seq Analysis Pipeline

## Overview

Fish Ig-seq Analysis Pipeline is a reproducible workflow for processing and analysing fish immunoglobulin sequencing data. It supports UMI-based read correction, clonotype assembly, V(D)J gene assignment, IGHV usage profiling, somatic hypermutation analysis, and comparative repertoire analysis across samples, tissues, immune conditions, and species.

This pipeline was developed for the analysis of newly generated zebrafish Ig-seq datasets and previously published Ig-seq datasets from humans, mice, and rainbow trout.

## Data sources

Ig-seq datasets generated in this study from head kidney samples of AB and TU zebrafish are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA1027976.

Previously published Ig-seq datasets from humans, mice, and rainbow trout were obtained from the NCBI SRA and reanalysed in this study.

IGHV sequences and annotations were retrieved from the IMGT database.

## System requirements

### Operating systems

The pipeline has been tested on:

- Ubuntu 20.04 or later
- Ubuntu 22.04

### Software dependencies

The following software is required:

- MiXCR v4.7.0
- MIGEC v1.2.9
- R v4.4.1
- Python v3.10 or later

Required R packages include:

- tidyverse
- data.table
- ggplot2
- reshape2
- vegan

Required Python packages include:

- pandas
- numpy
- biopython

No non-standard hardware is required. The demo dataset can be run on a standard desktop or laptop computer.

## Installation

Clone this repository:

```bash
git clone https://github.com/zhangh276/Fish-Ig-seq-Analysis-Pipeline.git
cd Fish-Ig-seq-Analysis-Pipeline
