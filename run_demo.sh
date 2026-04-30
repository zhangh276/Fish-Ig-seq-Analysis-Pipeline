#!/usr/bin/env bash
set -euo pipefail

# Run the full demo workflow for the Fish Ig-seq Analysis Pipeline.

# ============================================================
# Software paths
#
# Please modify these paths according to your local installation.
#
# MIXCR_HOME should be the MiXCR installation directory, not
# the mixcr executable itself.
#
# Example:
#   export MIXCR_HOME=/path/to/mixcr-4.7.0
#   export MIGEC_JAR=/path/to/migec-1.2.9.jar
#   bash run_demo.sh metadata.csv
# ============================================================

export MIXCR_HOME="${MIXCR_HOME:-/path/to/mixcr-4.7.0}"
export MIXCR_BIN="${MIXCR_BIN:-${MIXCR_HOME}/mixcr}"
export MIGEC_JAR="${MIGEC_JAR:-/path/to/migec-1.2.9.jar}"

# ============================================================
# MiXCR reference settings
#
# Step 00 builds Dr-AB-IGH.json.gz and copies it into:
#   ${MIXCR_HOME}/libraries/
#
# Downstream MiXCR commands use the library name, not a path.
# ============================================================

export MIXCR_SPECIES="${MIXCR_SPECIES:-zebrafish}"
export MIXCR_LIBRARY="${MIXCR_LIBRARY:-Dr-AB-IGH}"

# ============================================================
# Input metadata
# ============================================================

METADATA="${1:-metadata.csv}"

# ============================================================
# Check required software and input files
# ============================================================

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: metadata file not found: ${METADATA}"
  exit 1
fi

if [[ ! -d "${MIXCR_HOME}" ]]; then
  echo "ERROR: MiXCR installation directory not found:"
  echo "${MIXCR_HOME}"
  echo "Please modify MIXCR_HOME in run_demo.sh or export it before running."
  exit 1
fi

if [[ ! -x "${MIXCR_BIN}" ]]; then
  echo "ERROR: MiXCR executable not found or not executable:"
  echo "${MIXCR_BIN}"
  echo "Please modify MIXCR_BIN in run_demo.sh or export it before running."
  exit 1
fi

if [[ ! -d "${MIXCR_HOME}/libraries" ]]; then
  echo "ERROR: MiXCR library directory not found:"
  echo "${MIXCR_HOME}/libraries"
  echo "Please check MIXCR_HOME. This workflow will not create the MiXCR library directory automatically."
  exit 1
fi

if [[ ! -f "${MIGEC_JAR}" ]]; then
  echo "ERROR: MIGEC jar file not found:"
  echo "${MIGEC_JAR}"
  echo "Please modify MIGEC_JAR in run_demo.sh or export it before running."
  exit 1
fi

if ! command -v java >/dev/null 2>&1; then
  echo "ERROR: java is not available in PATH."
  exit 1
fi

if ! command -v seqtk >/dev/null 2>&1; then
  echo "ERROR: seqtk is not available in PATH."
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available in PATH."
  exit 1
fi

echo "Using metadata: ${METADATA}"
echo "Using MiXCR home: ${MIXCR_HOME}"
echo "Using MiXCR binary: ${MIXCR_BIN}"
echo "Using MiXCR library directory: ${MIXCR_HOME}/libraries"
echo "Using MIGEC: ${MIGEC_JAR}"
echo "Using MiXCR species: ${MIXCR_SPECIES}"
echo "Using MiXCR library name: ${MIXCR_LIBRARY}"

# ============================================================
# Install/check required R packages
# ============================================================

echo "============================================================"
echo "Checking required R packages"
echo "============================================================"
Rscript scripts/install_R_packages.R

# ============================================================
# Clean previous demo results
# ============================================================

RESULTS_DIR="results"

if [[ -d "${RESULTS_DIR}" ]]; then
  echo "Removing previous results directory: ${RESULTS_DIR}"
  rm -rf "${RESULTS_DIR}"
fi

mkdir -p "${RESULTS_DIR}"

# ============================================================
# Run workflow
# ============================================================

echo "============================================================"
echo "Step 00: Build custom MiXCR reference library"
echo "============================================================"
bash scripts/00_build_mixcr_library.sh \
  example_data/reference \
  results/reference_library

echo "============================================================"
echo "Step 01: MIGEC UMI correction"
echo "============================================================"
bash scripts/01_migec_umi_correction.sh \
  "${METADATA}" \
  results/01_migec

echo "============================================================"
echo "Step 02: MiXCR alignment and CDR3 clonotype assembly"
echo "============================================================"
bash scripts/02_mixcr_alignment_assemble_cdr3_clones.sh \
  "${METADATA}" \
  results/01_migec \
  results/02_cdr3_clones

echo "============================================================"
echo "Step 03: Export CDR3 clonotypes as TSV"
echo "============================================================"
bash scripts/03_export_cdr3_clones_tsv.sh \
  "${METADATA}" \
  results/02_cdr3_clones \
  results/03_cdr3_tsv

echo "============================================================"
echo "Step 04: MiXCR alignment and VDJRegion clonotype assembly"
echo "============================================================"
bash scripts/04_mixcr_alignment_assemble_vdjregion_clones.sh \
  "${METADATA}" \
  results/01_migec \
  results/04_vdjregion_clones

echo "============================================================"
echo "Step 05: Export VDJRegion clonotypes and mutation tables as TSV"
echo "============================================================"
bash scripts/05_export_vdjregion_clones_and_mutation_tsv.sh \
  "${METADATA}" \
  results/04_vdjregion_clones \
  results/05_vdjregion_mutation_tsv

echo "============================================================"
echo "Step 06: IGHV and VJ usage analysis"
echo "============================================================"
Rscript scripts/06_ighv_vj_usage_analysis.R \
  "${METADATA}" \
  results/03_cdr3_tsv \
  example_data/reference \
  results/06_ighv_vj_usage

echo "============================================================"
echo "Step 07: Pairwise repertoire correlation analysis"
echo "============================================================"
Rscript scripts/07_pairwise_correlation_analysis.R \
  results/06_ighv_vj_usage \
  results/07_repertoire_correlation

echo "============================================================"
echo "Step 08: CDR3 clonotype network analysis"
echo "============================================================"
Rscript scripts/08_cdr3_clonotype_network_analysis.R \
  "${METADATA}" \
  results/03_cdr3_tsv \
  results/08_cdr3_network \
  3000

echo "============================================================"
echo "Step 09: V-region mutation analysis"
echo "============================================================"
Rscript scripts/09_v_region_mutation_analysis.R \
  "${METADATA}" \
  results/05_vdjregion_mutation_tsv \
  results/09_v_region_mutation

echo "============================================================"
echo "Demo workflow finished successfully."
echo "Output files are available in the results/ directory."
echo "============================================================"
