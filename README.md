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

Data sources

Ig-seq datasets generated in this study from head kidney samples of AB and TU zebrafish are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA1027976.

Previously published Ig-seq datasets from humans, mice, and rainbow trout were obtained from the NCBI SRA and reanalysed in this study. These datasets include human PBMC Ig-seq datasets, mouse PBMC, spleen, and bone marrow Ig-seq datasets, and rainbow trout kidney and spleen Ig-seq datasets.

IGHV sequences and annotations were retrieved from the IMGT database.

System requirements
Operating systems

The pipeline was developed and tested on Linux-based systems.

Tested operating system:

Ubuntu 22.04

The pipeline is expected to run on other Linux distributions with the required software dependencies installed. macOS may also be supported for some analysis steps, although it has not been fully tested.

Required software

The following software is required:

MiXCR v4.7.0
MIGEC
R v4.4.1
Python v3.10 or later
Bash

Depending on the specific analysis module used, additional command-line tools may be required, such as:

FastQC
cutadapt
seqkit
Required R packages

The following R packages are used in downstream repertoire analyses and visualization:

tidyverse
data.table
ggplot2
reshape2
dplyr
stringr
vegan
pheatmap
ape
Required Python packages

The following Python packages may be required for data formatting and sequence processing:

pandas
numpy
biopython
Hardware requirements

No non-standard hardware is required. The demo dataset can be run on a standard desktop or laptop computer.

For large-scale Ig-seq datasets, higher memory and multi-core CPUs are recommended. Actual requirements depend on the size of the FASTQ files and the number of samples analysed.

Installation

Clone this repository:

git clone https://github.com/zhangh276/Fish-Ig-seq-Analysis-Pipeline.git
cd Fish-Ig-seq-Analysis-Pipeline

Install the required dependencies according to the official documentation of MiXCR, MIGEC, R, and Python.

If a conda environment file is provided, the environment can be created using:

conda env create -f environment.yml
conda activate fish_igseq

Alternatively, install required Python packages using:

pip install -r requirements.txt

R packages can be installed within R using:

install.packages(c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "reshape2",
  "dplyr",
  "stringr",
  "vegan",
  "pheatmap",
  "ape"
))

Typical installation time is approximately 10–30 minutes on a standard desktop computer, depending on internet speed, package manager configuration, and software availability.

Demo

A small demo dataset is provided in the example_data/ directory to demonstrate the basic functionality of the pipeline.

To run the demo:

bash run_demo.sh

The demo processes example Ig-seq input files and generates output files in the example_output/ directory.

Expected demo output

The demo is expected to generate the following types of output files:

example_output/
├── processed_reads/
├── clonotype_tables/
├── vdj_assignment/
├── ighv_usage/
├── shm_summary/
├── diversity_summary/
└── plots/

Expected output includes:

Processed read files
Clonotype tables
V(D)J gene assignment summaries
IGHV usage tables
Somatic hypermutation summary tables
Repertoire diversity summaries
Example repertoire comparison plots
Expected demo run time

The demo dataset is expected to run within approximately 5–20 minutes on a standard desktop or laptop computer.

Actual run time may vary depending on hardware, software configuration, and the size of the demo input files.

Instructions for use

To analyse user-provided Ig-seq data, prepare the following input files:

Paired-end FASTQ files
Sample metadata table
Germline V(D)J reference sequences
Primer and barcode information, if applicable

The general workflow includes the following steps:

1. Preprocessing of raw sequencing reads
2. UMI-guided read correction where applicable
3. Clonotype assembly
4. V(D)J gene assignment
5. IGHV usage quantification
6. Somatic hypermutation analysis
7. Repertoire diversity analysis
8. Repertoire similarity and correlation analysis
9. Cross-sample, cross-condition, and cross-species comparative analysis

A typical command for running the pipeline is:

bash run_pipeline.sh example_data/metadata.tsv example_data/reference/

Users should modify the metadata file and reference sequence files according to their own datasets.

Input file format
Metadata file

The metadata file should be a tab-delimited file containing sample information and paths to input FASTQ files.

Example:

sample_id    species        strain    tissue       condition    R1_fastq                         R2_fastq
ZF_AB_01     Danio_rerio    AB        head_kidney  resting      fastq/ZF_AB_01_R1.fastq.gz        fastq/ZF_AB_01_R2.fastq.gz
ZF_TU_01     Danio_rerio    TU        head_kidney  resting      fastq/ZF_TU_01_R1.fastq.gz        fastq/ZF_TU_01_R2.fastq.gz

Recommended columns include:

sample_id
species
strain
tissue
condition
isotype
R1_fastq
R2_fastq

Additional columns can be included as needed, such as:

dataset
BioProject
BioSample
SRA_accession
treatment
time_point
replicate
Reference files

Reference files should include germline V(D)J gene sequences in FASTA format.

Example:

reference/
├── IGHV.fasta
├── IGHD.fasta
├── IGHJ.fasta
└── IGHC.fasta

For zebrafish and other fish species, custom curated V(D)J references may be required because public annotations are incomplete for some loci and species.

Output files

The main output files include:

UMI-corrected read files
MiXCR clonotype tables
V(D)J assignment summaries
IGHV usage matrices
Somatic hypermutation summary tables
Diversity estimates
Repertoire correlation matrices
Cross-sample comparison results
Summary plots

Example output directories:

processed_reads/
clonotype_tables/
vdj_assignment/
ighv_usage/
shm_summary/
diversity_summary/
comparative_analysis/
plots/
Workflow description
1. Raw-read preprocessing

Raw paired-end FASTQ files are processed to remove low-quality reads and extract sample-specific barcode and UMI information where applicable.

2. UMI-guided read correction

For UMI-containing libraries, UMI-guided error correction is performed to reduce sequencing and PCR errors before clonotype assembly.

3. Clonotype assembly

Processed reads are assembled into clonotypes using MiXCR or other appropriate repertoire assembly tools. Clonotypes are defined according to V(D)J assignment and CDR3 sequence information.

4. V(D)J gene assignment

Reads and clonotypes are assigned to germline V, D, and J genes using species-specific or custom-curated reference sequences.

5. IGHV usage analysis

IGHV gene usage is quantified for each sample. The pipeline can generate IGHV usage profiles at the read, UMI, clonotype, or V3J-type level, depending on the dataset and analysis design.

6. Somatic hypermutation analysis

Somatic hypermutation levels are estimated by comparing assembled sequences with assigned germline IGHV sequences. Mutation burden can be summarized across samples, groups, tissues, conditions, species, or datasets.

7. Diversity and similarity analysis

Repertoire diversity and similarity are analysed using clonotype abundance, IGHV usage, VJ usage, and other repertoire-level features. Pairwise repertoire correlations and group-level comparisons can be performed.

8. Comparative repertoire analysis

The pipeline supports comparative analyses across individuals, strains, tissues, immune conditions, species, and independent datasets. The same standardized workflow is applied to newly generated and previously published Ig-seq datasets to improve comparability.

Reproducibility

The same standardized computational workflow was applied to all newly generated and previously published Ig-seq datasets analysed in this study.

The repository provides code, example input files, demo instructions, and expected output files to help users test and reproduce the analysis workflow.

For full reproduction of the manuscript analyses, users should download the relevant datasets from NCBI SRA using the BioProject accession numbers described above and follow the analysis steps provided in this repository.

Notes on large datasets

Raw FASTQ files from NCBI SRA are not included in this repository because of file size limitations. Users can download the raw sequencing data directly from NCBI SRA.

Example command using the SRA Toolkit:

prefetch SRR_ACCESSION
fasterq-dump SRR_ACCESSION --split-files

The resulting FASTQ files should be placed in the appropriate input directory and listed in the metadata file.

License

This project is released under the MIT License. See the LICENSE file for details.

Citation

If you use this pipeline, please cite the associated manuscript.

Suggested citation format will be added after publication.

Contact

For questions, please contact the corresponding authors or open an issue in this GitHub repository.

Repository:

https://github.com/zhangh276/Fish-Ig-seq-Analysis-Pipeline/
