#!/usr/bin/env bash
set -euo pipefail

# Step 04: MiXCR alignment and VDJRegion-based clonotype assembly.
# This step uses mixcr analyze generic-amplicon.
# Input:
#   MIGEC-corrected reads from Step 01
# Output:
#   results/04_vdjregion_clones/<sample_id>_VDJRegion.clns

METADATA="${1:-metadata.csv}"
MIGEC_OUT_DIR="${2:-results/01_migec}"
OUT_DIR="${3:-results/04_vdjregion_clones}"

MIXCR_BIN="${MIXCR_BIN:-mixcr}"
MIXCR_SPECIES="${MIXCR_SPECIES:-zebrafish}"
MIXCR_LIBRARY="${MIXCR_LIBRARY:-Dr-AB-IGH}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: metadata file not found: ${METADATA}"
  exit 1
fi

tail -n +2 "${METADATA}" | while IFS=',' read -r sample_id strain condition R1_fastq R2_fastq
do
  echo "Running MiXCR VDJRegion clonotype analysis for sample: ${sample_id}"

  CONSENSUS_DIR="${MIGEC_OUT_DIR}/${sample_id}/consensus"

  R1_CORRECTED="${CONSENSUS_DIR}/${sample_id}_R1.t2.fastq"
  R2_RC="${CONSENSUS_DIR}/${sample_id}_R2.t2.rc.fastq"

  if [[ ! -f "${R1_CORRECTED}" ]]; then
    echo "ERROR: corrected R1 file not found: ${R1_CORRECTED}"
    exit 1
  fi

  if [[ ! -f "${R2_RC}" ]]; then
    echo "ERROR: reverse-complemented corrected R2 file not found: ${R2_RC}"
    exit 1
  fi

  "${MIXCR_BIN}" analyze generic-amplicon \
    --species "7955" \
    --library "${MIXCR_LIBRARY}" \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    --split-clones-by V \
    --split-clones-by J \
    --assemble-clonotypes-by VDJRegion \
    "${R1_CORRECTED}" \
    "${R2_RC}" \
    "${OUT_DIR}/${sample_id}_VDJRegion"

  echo "VDJRegion clonotype analysis finished for sample: ${sample_id}"
done

echo "Step 04 finished: MiXCR VDJRegion clonotype analysis completed."
