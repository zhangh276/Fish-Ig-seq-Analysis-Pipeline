#!/usr/bin/env bash
set -euo pipefail

# Step 01: UMI-guided read correction using MIGEC.
# Input:
#   metadata.csv
# Output:
#   results/01_migec/<sample_id>/consensus/<sample_id>_R1.t2.fastq
#   results/01_migec/<sample_id>/consensus/<sample_id>_R2.t2.fastq
#   results/01_migec/<sample_id>/consensus/<sample_id>_R2.t2.rc.fastq

METADATA="${1:-metadata.csv}"
OUT_DIR="${2:-results/01_migec}"

MIGEC_JAR="${MIGEC_JAR:-migec-1.2.9.jar}"
UMI_PATTERN="${UMI_PATTERN:-NNNNNNNNNNNNGTACGG}"

mkdir -p "${OUT_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  echo "ERROR: metadata file not found: ${METADATA}"
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

if [[ ! -f "${MIGEC_JAR}" ]]; then
  echo "ERROR: MIGEC jar file not found: ${MIGEC_JAR}"
  echo "Please set MIGEC_JAR, for example:"
  echo "export MIGEC_JAR=/path/to/migec-1.2.9.jar"
  exit 1
fi

tail -n +2 "${METADATA}" | while IFS=',' read -r sample_id strain condition R1_fastq R2_fastq
do
  echo "Processing sample with MIGEC: ${sample_id}"

  SAMPLE_DIR="${OUT_DIR}/${sample_id}"
  BARCODE_FILE="${SAMPLE_DIR}/barcode.txt"
  CHECKOUT_DIR="${SAMPLE_DIR}/Checkout"
  CONSENSUS_DIR="${SAMPLE_DIR}/consensus"

  mkdir -p "${SAMPLE_DIR}" "${CHECKOUT_DIR}" "${CONSENSUS_DIR}"

  if [[ ! -f "${R1_fastq}" ]]; then
    echo "ERROR: R1 FASTQ file not found: ${R1_fastq}"
    exit 1
  fi

  if [[ ! -f "${R2_fastq}" ]]; then
    echo "ERROR: R2 FASTQ file not found: ${R2_fastq}"
    exit 1
  fi

  # MIGEC barcode file format used in the original analysis:
  # sample_id    UMI_pattern    empty_column    R2_fastq    R1_fastq
  printf "%s\t%s\t\t%s\t%s\n" \
    "${sample_id}" \
    "${UMI_PATTERN}" \
    "${R2_fastq}" \
    "${R1_fastq}" \
    > "${BARCODE_FILE}"

  java -jar "${MIGEC_JAR}" CheckoutBatch \
    -u \
    -t \
    -o "${BARCODE_FILE}" \
    "${CHECKOUT_DIR}/"

  java -jar "${MIGEC_JAR}" Assemble \
    -m 2 \
    "${CHECKOUT_DIR}/${sample_id}_R1.fastq" \
    "${CHECKOUT_DIR}/${sample_id}_R2.fastq" \
    "${CONSENSUS_DIR}/"

  R2_CORRECTED="${CONSENSUS_DIR}/${sample_id}_R2.t2.fastq"
  R2_RC="${CONSENSUS_DIR}/${sample_id}_R2.t2.rc.fastq"

  if [[ ! -f "${R2_CORRECTED}" ]]; then
    echo "ERROR: MIGEC corrected R2 file not found: ${R2_CORRECTED}"
    exit 1
  fi

  # Generate reverse-complemented R2 file for MiXCR analysis.
  # The original R2.t2.fastq is kept unchanged.
  seqtk seq -r "${R2_CORRECTED}" > "${R2_RC}"

  echo "MIGEC correction finished for sample: ${sample_id}"
done

echo "Step 01 finished: MIGEC UMI correction completed."
