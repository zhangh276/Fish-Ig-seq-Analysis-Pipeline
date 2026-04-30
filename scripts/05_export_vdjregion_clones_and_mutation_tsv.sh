#!/usr/bin/env bash
set -euo pipefail

# Step 05: Export VDJRegion clonotypes and mutation information as TSV files.
# Input:
#   results/04_vdjregion_clones/<sample_id>_VDJRegion.clns
# Output:
#   results/05_vdjregion_mutation_tsv/<sample_id>_VDJRegion_clones.tsv
#   results/05_vdjregion_mutation_tsv/<sample_id>_mutation.tsv

METADATA="${1:-metadata.csv}"
VDJ_DIR="${2:-results/04_vdjregion_clones}"
OUT_DIR="${3:-results/05_vdjregion_mutation_tsv}"

MIXCR_BIN="${MIXCR_BIN:-mixcr}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: metadata file not found: ${METADATA}"
  exit 1
fi

tail -n +2 "${METADATA}" | while IFS=',' read -r sample_id strain condition R1_fastq R2_fastq
do
  echo "Exporting VDJRegion clonotype and mutation tables for sample: ${sample_id}"

  CLNS="${VDJ_DIR}/${sample_id}_VDJRegion.clns"

  VDJ_TSV="${OUT_DIR}/${sample_id}_VDJRegion_clones.tsv"
  MUT_TSV="${OUT_DIR}/${sample_id}_mutation.tsv"

  if [[ ! -f "${CLNS}" ]]; then
    echo "ERROR: VDJRegion clns file not found: ${CLNS}"
    exit 1
  fi

  "${MIXCR_BIN}" exportClones \
    -f \
    "${CLNS}" \
    "${VDJ_TSV}"

  "${MIXCR_BIN}" exportClones \
    -f \
    -nLength FR1 \
    -nLength CDR1 \
    -nLength FR2 \
    -nLength CDR2 \
    -nLength FR3 \
    -nMutationsCount FR1 substitutions \
    -nMutationsCount CDR1 substitutions \
    -nMutationsCount FR2 substitutions \
    -nMutationsCount CDR2 substitutions \
    -nMutationsCount FR3 substitutions \
    "${CLNS}" \
    "${MUT_TSV}"

  echo "VDJRegion clonotype TSV generated: ${VDJ_TSV}"
  echo "Mutation TSV generated: ${MUT_TSV}"
done

echo "Step 05 finished: VDJRegion and mutation TSV export completed."
