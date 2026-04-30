et -euo pipefail

# Step 00: Build a custom MiXCR reference library from V/D/J/C FASTA files.
#
# Input:
#   example_data/reference/v-genes.IGH.fasta
#   example_data/reference/d-genes.IGH.fasta
#   example_data/reference/j-genes.IGH.fasta
#   example_data/reference/c-genes.IGH.fasta
#
# Temporary output:
#   results/reference_library/Dr-AB-IGH.json.gz
#
# Final library location:
#   ${MIXCR_HOME}/libraries/Dr-AB-IGH.json.gz
#
# Downstream MiXCR command should use:
#   --library Dr-AB-IGH

# MiXCR installation directory and executable.
MIXCR_HOME="${MIXCR_HOME:-/RAID5_32T/zh/software/MIXCR/mixcr-4.7.0}"
MIXCR_BIN="${MIXCR_BIN:-${MIXCR_HOME}/mixcr}"

# Input reference directory and temporary output directory.
REF_DIR="${1:-example_data/reference}"
OUT_DIR="${2:-results/reference_library}"

# Reference FASTA files.
V_FASTA="${REF_DIR}/v-genes.IGH.fasta"
D_FASTA="${REF_DIR}/d-genes.IGH.fasta"
J_FASTA="${REF_DIR}/j-genes.IGH.fasta"
C_FASTA="${REF_DIR}/c-genes.IGH.fasta"

# Library name and output file.
LIBRARY_NAME="Dr-AB-IGH"
LIBRARY_OUT="${OUT_DIR}/${LIBRARY_NAME}.json.gz"

# Existing MiXCR library directory.
MIXCR_LIBRARY_DIR="${MIXCR_HOME}/libraries"
FINAL_LIBRARY="${MIXCR_LIBRARY_DIR}/${LIBRARY_NAME}.json.gz"

mkdir -p "${OUT_DIR}"

# Check MiXCR executable.
if [[ ! -x "${MIXCR_BIN}" ]]; then
  echo "ERROR: MiXCR executable not found or not executable:"
  echo "${MIXCR_BIN}"
  echo "Please set MIXCR_HOME or MIXCR_BIN correctly."
  exit 1
fi

# Check existing MiXCR library directory.
if [[ ! -d "${MIXCR_LIBRARY_DIR}" ]]; then
  echo "ERROR: MiXCR library directory does not exist:"
  echo "${MIXCR_LIBRARY_DIR}"
  echo "Please check MIXCR_HOME. This script will not create the MiXCR library directory automatically."
  exit 1
fi

# Check required FASTA files.
if [[ ! -f "${V_FASTA}" ]]; then
  echo "ERROR: V gene FASTA file not found:"
  echo "${V_FASTA}"
  exit 1
fi

if [[ ! -f "${J_FASTA}" ]]; then
  echo "ERROR: J gene FASTA file not found:"
  echo "${J_FASTA}"
  exit 1
fi

# Build MiXCR library command.
CMD=(
  "${MIXCR_BIN}" buildLibrary
  --force-overwrite
  --debug
  --v-genes-from-fasta "${V_FASTA}"
  --v-gene-feature VRegion
  --j-genes-from-fasta "${J_FASTA}"
  --chain IGH
  --taxon-id 7955
  --species Dr_AB_IGH
)

# D and C FASTA files are optional.
if [[ -f "${D_FASTA}" ]]; then
  CMD+=(--d-genes-from-fasta "${D_FASTA}")
else
  echo "WARNING: D gene FASTA file not found, skipping D genes:"
  echo "${D_FASTA}"
fi

if [[ -f "${C_FASTA}" ]]; then
  CMD+=(--c-genes-from-fasta "${C_FASTA}")
else
  echo "WARNING: C gene FASTA file not found, skipping C genes:"
  echo "${C_FASTA}"
fi

CMD+=("${LIBRARY_OUT}")

echo "Building custom MiXCR library:"
echo "${CMD[@]}"

"${CMD[@]}"

echo "Custom MiXCR library generated:"
echo "${LIBRARY_OUT}"

# Copy the custom library to the existing MiXCR library directory.
cp "${LIBRARY_OUT}" "${FINAL_LIBRARY}"

echo "Custom MiXCR library copied to:"
echo "${FINAL_LIBRARY}"

echo "Use the following MiXCR library name in downstream steps:"
echo "${LIBRARY_NAME}"
