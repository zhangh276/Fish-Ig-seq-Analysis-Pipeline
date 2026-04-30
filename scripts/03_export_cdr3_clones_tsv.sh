#!/usr/bin/env bash
set -euo pipefail

# Step 03: Export CDR3-based clonotypes as TSV files.
# Input:
#   results/02_cdr3_clones/<sample_id>_CDR3.clns
# Output:
#   results/03_cdr3_tsv/<sample_id>_CDR3_clones.tsv

METADATA="${1:-metadata.csv}"
CDR3_DIR="${2:-results/02_cdr3_clones}"
OUT_DIR="${3:-results/03_cdr3_tsv}"

MIXCR_BIN="${MIXCR_BIN:-mixcr}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: metadata file not found: ${METADATA}"
  exit 1
fi

tail -n +2 "${METADATA}" | while IFS=',' read -r sample_id strain condition R1_fastq R2_fastq
do
  echo "Exporting CDR3 clonotype table for sample: ${sample_id}"

  CLNS="${CDR3_DIR}/${sample_id}_CDR3.clns"
  TSV="${OUT_DIR}/${sample_id}_CDR3_clones.tsv"

  if [[ ! -f "${CLNS}" ]]; then
    echo "ERROR: CDR3 clns file not found: ${CLNS}"
    exit 1
  fi

  "${MIXCR_BIN}" exportClones \
    -f \
    "${CLNS}" \
    "${TSV}"

  echo "CDR3 clonotype TSV generated: ${TSV}"
done

echo "Step 03 finished: CDR3 clonotype TSV export completed."
