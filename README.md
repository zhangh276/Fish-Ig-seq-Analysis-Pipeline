# Fish Ig-seq Analysis Pipeline

## Overview

Fish Ig-seq Analysis Pipeline is a reproducible workflow for processing and analysing immunoglobulin sequencing (Ig-seq) data from fish and mammals. The pipeline was developed for the analysis of newly generated zebrafish Ig-seq datasets and the reanalysis of previously published Ig-seq datasets from humans, mice, and rainbow trout.

The workflow supports UMI-based read correction, clonotype assembly, V(D)J gene assignment, IGHV usage profiling, somatic hypermutation analysis, repertoire diversity estimation, and comparative repertoire analysis across samples, tissues, immune conditions, strains, species, and datasets.

## Repository structure

```text
Fish-Ig-seq-Analysis-Pipeline/
├── README.md
├── LICENSE
├── scripts/
│   ├── preprocessing/
│   ├── mixcr_analysis/
│   ├── ighv_usage/
│   ├── shm_analysis/
│   ├── diversity_analysis/
│   └── comparative_analysis/
├── example_data/
│   ├── metadata.tsv
│   ├── fastq/
│   └── reference/
├── example_output/
├── run_demo.sh
└── run_pipeline.sh

The exact structure of the repository may be updated as the analysis workflow is further improved and documented.
