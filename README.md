# Fish Ig-seq Analysis Pipeline

## Overview

Fish Ig-seq Analysis Pipeline is an analysis workflow specifically designed for fish Ig-seq data and is tailored to our fish repertoire library preparation strategy based on 5′ RACE and UMI-containing reads. This pipeline supports UMI-based read correction, custom MiXCR reference library construction, clonotype assembly, and V(D)J gene assignment, followed by downstream characterization of repertoire features, including IGHV gene usage and VJ combination analysis, Pearson correlation analysis, CDR3 clonotype network analysis, and V-region somatic hypermutation analysis.

This project uses zebrafish Ig-seq data as an example to demonstrate our analysis workflow. The pipeline can also be extended to Ig-seq datasets from other fish species and mammalian species. This standardized analytical framework facilitates comparative analyses of antibody repertoire features across species.

This repository provides custom scripts, demo data, reference files, and instructions for running the analysis pipeline.

## Repository structure

```text
Fish-Ig-seq-Analysis-Pipeline/
├── README.md
├── LICENSE
├── metadata.csv
├── run_demo.sh
├── scripts/
│   ├── 00_build_mixcr_library.sh
│   ├── 01_migec_umi_correction.sh
│   ├── 02_mixcr_alignment_assemble_cdr3_clones.sh
│   ├── 03_export_cdr3_clones_tsv.sh
│   ├── 04_mixcr_alignment_assemble_vdjregion_clones.sh
│   ├── 05_export_vdjregion_clones_and_mutation_tsv.sh
│   ├── 06_ighv_vj_usage_analysis.R
│   ├── 07_pairwise_correlation_analysis.R
│   ├── 08_cdr3_clonotype_network_analysis.R
│   ├── 09_v_region_mutation_analysis.R
│   └── install_R_packages.R
├── example_data/
│   ├── fastq/
│   │   ├── demo1_R1.fastq.gz
│   │   ├── demo1_R2.fastq.gz
│   │   ├── demo2_R1.fastq.gz
│   │   ├── demo2_R2.fastq.gz
│   │   ├── demo3_R1.fastq.gz
│   │   └── demo3_R2.fastq.gz
│   └── reference/
│       ├── v-genes.IGH.fasta
│       ├── d-genes.IGH.fasta
│       ├── j-genes.IGH.fasta
│       └── c-genes.IGH.fasta
└── results/
```

The `results/` directory is generated automatically when the demo workflow is run. Existing `results/` files are removed at the beginning of each demo run.

## Data sources

Ig-seq datasets generated in this study from head kidney samples of AB and TU zebrafish are available in the NCBI Sequence Read Archive (SRA) under BioProject accession number PRJNA1027976. Three demo datasets from individual zebrafish IgM repertoires are provided in the `example_data/fastq/`.

Zebrafish IGHV, IGHD, IGHJ, and IGHC reference sequences used in the demo workflow are provided in `example_data/reference/`. 

## System requirements

### Operating system

The workflow was developed and tested on a Linux system.

Tested system:

```text
CentOS Linux
```

The pipeline is expected to run on other Linux distributions with the required dependencies installed. macOS may support some steps but has not been fully tested.

### Required software

The following software is required:

```text
Java
MiXCR (tested version: v4.7.0)
MIGEC (tested version: v1.2.9)
seqtk
R
```

The workflow uses shell scripts and R scripts. Bash-compatible shell is required.

### Required R packages

The following R packages are required for the downstream R analyses:

```text
stringdist
igraph
dplyr
tidyr
data.table
stringr
```

These packages can be installed using:

```bash
Rscript scripts/install_R_packages.R
```

Alternatively, they can be installed manually in R:

```r
install.packages(c(
  "stringdist",
  "igraph",
  "dplyr",
  "tidyr",
  "data.table",
  "stringr"
), repos = "https://cloud.r-project.org")
```

### Hardware requirements

No non-standard hardware is required for the demo workflow.

The demo dataset can be run on a standard desktop or laptop computer. For full-scale Ig-seq datasets, higher memory and multi-core CPUs are recommended because MiXCR analysis and CDR3 network construction can be computationally intensive.

## Installation

Clone this repository:

```bash
git clone https://github.com/zhangh276/Fish-Ig-seq-Analysis-Pipeline.git
cd Fish-Ig-seq-Analysis-Pipeline
```

Make the shell scripts executable:

```bash
chmod +x run_demo.sh
chmod +x scripts/*.sh
```

Install required R packages:

```bash
Rscript scripts/install_R_packages.R
```

MiXCR and MIGEC should be installed separately according to their official documentation. After installation, specify their paths in `run_demo.sh` or export the paths before running the demo.

Example:

```bash
export MIXCR_HOME=/path/to/mixcr-4.7.0
export MIGEC_JAR=/path/to/migec-1.2.9.jar
```

In this workflow, `MIXCR_HOME` should point to the MiXCR installation directory, not the `mixcr` executable itself. For example:

```text
/path/to/mixcr-4.7.0
```

The MiXCR executable is expected to be:

```text
/path/to/mixcr-4.7.0/mixcr
```

The MiXCR library directory is expected to be:

```text
/path/to/mixcr-4.7.0/libraries
```

The custom reference library generated in Step 00 will be copied into this `libraries/` directory.

## Typical installation time

Typical installation time is approximately 10-20 minutes on a standard desktop computer, depending on internet speed, package manager configuration, and whether MiXCR and MIGEC have already been installed.

## Demo dataset

A small demo dataset is provided in `example_data/fastq/`, including three paired-end Ig-seq samples:

```text
demo1_R1.fastq.gz
demo1_R2.fastq.gz
demo2_R1.fastq.gz
demo2_R2.fastq.gz
demo3_R1.fastq.gz
demo3_R2.fastq.gz
```

The sample metadata file is provided as:

```text
metadata.csv
```

Example metadata format:

```csv
sample_id,strain,condition,R1_fastq,R2_fastq
demo1,AB,demo,example_data/fastq/demo1_R1.fastq.gz,example_data/fastq/demo1_R2.fastq.gz
demo2,AB,demo,example_data/fastq/demo2_R1.fastq.gz,example_data/fastq/demo2_R2.fastq.gz
demo3,AB,demo,example_data/fastq/demo3_R1.fastq.gz,example_data/fastq/demo3_R2.fastq.gz
```

The demo data are intended only to demonstrate and test the workflow. They should not be used for biological interpretation.

## Running the demo

Before running the demo, set the MiXCR and MIGEC paths.

Example:

```bash
export MIXCR_HOME=/path/to/mixcr-4.7.0
export MIGEC_JAR=/path/to/migec-1.2.9.jar
```

Then run:

```bash
bash run_demo.sh metadata.csv
```

The demo workflow automatically removes any previous `results/` directory and generates a new one.

## Expected demo run time

Expected run time for the demo dataset is approximately 10 minutes on a standard desktop or laptop computer.

Run time depends on the number of reads in the demo FASTQ files and the number of CDR3 clonotypes used for the network analysis.

## Workflow description

The demo workflow contains the following steps.

### Step 00: Build custom MiXCR reference library

Script:

```text
scripts/00_build_mixcr_library.sh
```

This step builds a custom MiXCR IGH reference library from V, D, J, and C FASTA files.

Input files:

```text
example_data/reference/v-genes.IGH.fasta
example_data/reference/d-genes.IGH.fasta
example_data/reference/j-genes.IGH.fasta
example_data/reference/c-genes.IGH.fasta
```

Output:

```text
results/reference_library/Dr-AB-IGH.json.gz
```

The generated library is also copied to:

```text
${MIXCR_HOME}/libraries/Dr-AB-IGH.json.gz
```

Downstream MiXCR commands use the library name:

```text
Dr-AB-IGH
```

### Step 01: UMI-guided read correction

Script:

```text
scripts/01_migec_umi_correction.sh
```

This step uses MIGEC to process UMI-containing reads and generate corrected consensus reads.

Input:

```text
metadata.csv
example_data/fastq/*.fastq.gz
```

Output:

```text
results/01_migec/<sample_id>/consensus/<sample_id>_R1.t2.fastq
results/01_migec/<sample_id>/consensus/<sample_id>_R2.t2.fastq
results/01_migec/<sample_id>/consensus/<sample_id>_R2.t2.rc.fastq
```

The reverse-complemented R2 file is used in downstream MiXCR analysis.

### Step 02: MiXCR alignment and CDR3 clonotype assembly

Script:

```text
scripts/02_mixcr_alignment_assemble_cdr3_clones.sh
```

This step uses `mixcr analyze generic-amplicon` to perform alignment and assemble CDR3-based clonotypes.

Main MiXCR parameters:

```text
--species zebrafish
--library Dr-AB-IGH
--rna
--rigid-left-alignment-boundary
--floating-right-alignment-boundary C
--split-clones-by V
--split-clones-by J
--assemble-clonotypes-by CDR3
```

Output:

```text
results/02_cdr3_clones/<sample_id>_CDR3.clns
```

### Step 03: Export CDR3 clonotypes as TSV

Script:

```text
scripts/03_export_cdr3_clones_tsv.sh
```

This step exports CDR3 clonotypes from MiXCR `.clns` files to TSV format.

Output:

```text
results/03_cdr3_tsv/<sample_id>_CDR3_clones_IGH.tsv
```

These files are used for IGHV/VJ usage analysis, repertoire correlation analysis, and CDR3 network analysis.

### Step 04: MiXCR alignment and VDJRegion clonotype assembly

Script:

```text
scripts/04_mixcr_alignment_assemble_vdjregion_clones.sh
```

This step uses `mixcr analyze generic-amplicon` to assemble VDJRegion-based clonotypes.

Main MiXCR parameters:

```text
--species zebrafish
--library Dr-AB-IGH
--rna
--rigid-left-alignment-boundary
--floating-right-alignment-boundary C
--split-clones-by V
--split-clones-by J
--assemble-clonotypes-by VDJRegion
```

Output:

```text
results/04_vdjregion_clones/<sample_id>_VDJRegion.clns
```

### Step 05: Export VDJRegion clonotypes and mutation tables

Script:

```text
scripts/05_export_vdjregion_clones_and_mutation_tsv.sh
```

This step exports VDJRegion clonotypes and mutation-related fields from MiXCR `.clns` files.

Output:

```text
results/05_vdjregion_mutation_tsv/<sample_id>_VDJRegion_clones.tsv
results/05_vdjregion_mutation_tsv/<sample_id>_mutation_IGH.tsv
```

Mutation-related exported fields include FR/CDR region lengths and substitution mutation counts.

### Step 06: IGHV and VJ usage analysis

Script:

```text
scripts/06_ighv_vj_usage_analysis.R
```

This step calculates IGHV gene usage and VJ combination usage from CDR3 clonotype tables.

Input:

```text
results/03_cdr3_tsv/<sample_id>_CDR3_clones_IGH.tsv
example_data/reference/v-genes.IGH.fasta
example_data/reference/j-genes.IGH.fasta
```

The script extracts all V and J gene names from the reference FASTA files. Gene allele information is removed by keeping the part before the first `*`.

For example:

```text
IGHV1-1*01 -> IGHV1-1
IGHJ2*01   -> IGHJ2
```

Two weighting schemes are used:

```text
V3J_weight: each CDR3 clonotype contributes one count
reads_weight: each CDR3 clonotype is weighted by readCount
```

Output:

```text
results/06_ighv_vj_usage/clone_stat.tsv
results/06_ighv_vj_usage/IGHV_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/IGHV_reads_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_reads_weight_matrix.tsv
```

### Step 07: Pairwise repertoire correlation analysis

Script:

```text
scripts/07_pairwise_correlation_analysis.R
```

This step calculates pairwise Pearson correlation coefficients between samples based on the four usage frequency matrices generated in Step 06.

Input:

```text
results/06_ighv_vj_usage/IGHV_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/IGHV_reads_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_reads_weight_matrix.tsv
```

Output:

```text
results/07_repertoire_correlation/IGHV_V3J_weight_correlation_matrix.tsv
results/07_repertoire_correlation/IGHV_reads_weight_correlation_matrix.tsv
results/07_repertoire_correlation/VJ_V3J_weight_correlation_matrix.tsv
results/07_repertoire_correlation/VJ_reads_weight_correlation_matrix.tsv
```

### Step 08: CDR3 clonotype network analysis

Script:

```text
scripts/08_cdr3_clonotype_network_analysis.R
```

This step constructs CDR3 clonotype networks using Levenshtein distance.

Network definition:

```text
Nodes: unique V3J clonotypes, defined by V gene + nSeqCDR3 + J gene
Edges: pairs of clonotypes with Levenshtein distance = 1 between nSeqCDR3 sequences
```

Input:

```text
results/03_cdr3_tsv/<sample_id>_CDR3_clones_IGH.tsv
```

Output:

```text
results/08_cdr3_network/<sample_id>_CDR3_network_nodes.tsv
results/08_cdr3_network/<sample_id>_CDR3_network_edges.tsv
results/08_cdr3_network/<sample_id>_network_LD1_CDR3_nSeq.pdf
results/08_cdr3_network/clone_network_stat.tsv
```

The default maximum number of clonotypes used for network construction is 3000 per sample. This value is set in `run_demo.sh` and can be modified.

### Step 09: V-region mutation analysis

Script:

```text
scripts/09_v_region_mutation_analysis.R
```

This step summarizes V-region mutation burden based on the mutation tables generated in Step 05.

Input:

```text
results/05_vdjregion_mutation_tsv/<sample_id>_mutation_IGH.tsv
```

Rows containing `region_not_covered` are removed before analysis.

The script calculates read-weighted mutation statistics:

```text
cloneCount
allReads
unmutatedReads
unmutatedReadsFraction
mutatedReads
mutatedReadsFraction
meanMutationsPerRead
```

Output:

```text
results/09_v_region_mutation/<sample_id>_V_region_mutation_filtered.tsv
results/09_v_region_mutation/v_region_mutation_stat.tsv
```

## Expected output

After running:

```bash
bash run_demo.sh metadata.csv
```

the following output directories should be generated:

```text
results/
├── reference_library/
├── 01_migec/
├── 02_cdr3_clones/
├── 03_cdr3_tsv/
├── 04_vdjregion_clones/
├── 05_vdjregion_mutation_tsv/
├── 06_ighv_vj_usage/
├── 07_repertoire_correlation/
├── 08_cdr3_network/
└── 09_v_region_mutation/
```

Key final output files include:

```text
results/06_ighv_vj_usage/IGHV_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/IGHV_reads_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_V3J_weight_matrix.tsv
results/06_ighv_vj_usage/VJ_reads_weight_matrix.tsv
results/07_repertoire_correlation/IGHV_V3J_weight_correlation_matrix.tsv
results/07_repertoire_correlation/IGHV_reads_weight_correlation_matrix.tsv
results/07_repertoire_correlation/VJ_V3J_weight_correlation_matrix.tsv
results/07_repertoire_correlation/VJ_reads_weight_correlation_matrix.tsv
results/08_cdr3_network/clone_network_stat.tsv
results/09_v_region_mutation/v_region_mutation_stat.tsv
```

## Running individual steps

Each step can also be run separately.

### Step 06 only

```bash
Rscript scripts/06_ighv_vj_usage_analysis.R metadata.csv results/03_cdr3_tsv example_data/reference results/06_ighv_vj_usage
```

### Step 07 only

```bash
Rscript scripts/07_pairwise_correlation_analysis.R results/06_ighv_vj_usage results/07_repertoire_correlation
```

### Step 08 only

```bash
Rscript scripts/08_cdr3_clonotype_network_analysis.R metadata.csv results/03_cdr3_tsv results/08_cdr3_network 3000
```

### Step 09 only

```bash
Rscript scripts/09_v_region_mutation_analysis.R metadata.csv results/05_vdjregion_mutation_tsv results/09_v_region_mutation
```

## Input file formats

### Metadata file

The metadata file should be a comma-separated file with at least the following columns:

```text
sample_id
strain
condition
R1_fastq
R2_fastq
```

Example:

```csv
sample_id,strain,condition,R1_fastq,R2_fastq
demo1,AB,demo,example_data/fastq/demo1_R1.fastq.gz,example_data/fastq/demo1_R2.fastq.gz
demo2,TU,demo,example_data/fastq/demo2_R1.fastq.gz,example_data/fastq/demo2_R2.fastq.gz
demo3,AB,demo,example_data/fastq/demo3_R1.fastq.gz,example_data/fastq/demo3_R2.fastq.gz
```

Additional columns can be included and will be retained in downstream summary files.

### Reference files

The reference directory should contain:

```text
v-genes.IGH.fasta
d-genes.IGH.fasta
j-genes.IGH.fasta
c-genes.IGH.fasta
```

These files are used to build a custom MiXCR reference library in Step 00.

### FASTQ files

The workflow expects paired-end FASTQ files listed in `metadata.csv`.

Example:

```text
example_data/fastq/demo1_R1.fastq.gz
example_data/fastq/demo1_R2.fastq.gz
```

## Notes on full-scale analysis

The raw FASTQ files for the full manuscript datasets are not included in this repository because of file size limitations. Users can download the full datasets directly from the NCBI SRA using the BioProject and SRR accession numbers provided in the manuscript and supplementary tables.

Example command using SRA Toolkit:

```bash
prefetch SRR_ACCESSION
fasterq-dump SRR_ACCESSION --split-files
```

The resulting FASTQ files should be listed in the metadata file and placed in the appropriate input directory.

## Reproducibility

The same standardized computational workflow was applied to the zebrafish Ig-seq data and to the previously published human/mouse/trout Ig-seq datasets reanalysed in the study.

This repository provides:

```text
custom analysis scripts
demo input data
reference FASTA files
demo metadata
expected output structure
instructions for running the full demo workflow
```

The demo workflow is intended to verify that the code and software environment are functional. Full-scale reproduction of all manuscript analyses requires downloading the corresponding raw sequencing datasets from NCBI SRA.

## Troubleshooting

### MiXCR library error

If MiXCR reports that the library name cannot be a path, make sure that the custom library has been copied into the MiXCR library directory:

```text
${MIXCR_HOME}/libraries/Dr-AB-IGH.json.gz
```

Downstream commands should use:

```text
--library Dr-AB-IGH
```

not:

```text
--library results/reference_library/Dr-AB-IGH.json.gz
```

### R packages not found

Install required R packages:

```bash
Rscript scripts/install_R_packages.R
```

### Rscript not found

Install R and ensure `Rscript` is available in the system PATH.

Check with:

```bash
Rscript --version
```

### seqtk not found

Install `seqtk` and ensure it is available in the system PATH.

Check with:

```bash
seqtk
```

### MIGEC jar not found

Set the path to the MIGEC jar file:

```bash
export MIGEC_JAR=/path/to/migec-1.2.9.jar
```

### MiXCR executable not found

Set the MiXCR installation directory:

```bash
export MIXCR_HOME=/path/to/mixcr-4.7.0
```

The workflow expects the executable at:

```text
${MIXCR_HOME}/mixcr
```

## License

This project is released under the MIT License. See the `LICENSE` file for details.

## Citation

If you use this pipeline, please cite the associated manuscript.

Suggested citation information will be added after publication.

## Contact

For questions or issues, please contact the corresponding authors or open an issue in this GitHub repository.

Repository:

```text
https://github.com/zhangh276/Fish-Ig-seq-Analysis-Pipeline/
```
