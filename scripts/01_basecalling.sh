#!/bin/bash

# Exit on error
set -e

# Check to make sure that the runtime_config is there.
if [ ! -f "scripts/runtime_config.sh" ]; then
  log_error "ERROR: runtime_config.sh not found. Please run 'setup' step first..."
  exit 1
fi

# Source the logging utility and runtime_config
# Assumes that logging.sh is in the utils folder of the current folder.
source "$(dirname "$0")/utils/logging.sh"
source "scripts/runtime_config.sh"

CONFIG_FILE=$1
DORADO_MODEL_NAME=$(yq e '.parameters.basecalling.base_model_name' "${CONFIG_FILE}")
MODEL_SPEED=$(yq e '.parameters.basecalling.model_speed' "${CONFIG_FILE}")
BASECALLING_MODIFICATIONS=$(yq e '.parameters.basecalling.basecalling_modifications' "${CONFIG_FILE}")
BATCHSIZE=$(yq e '.parameters.basecalling.batch_size' "${CONFIG_FILE}")
POD5_DIR=$(yq e '.paths.pod5_dir' "${CONFIG_FILE}")
POD5_NAME=$(yq e '.paths.pod5_name' "${CONFIG_FILE}")
BASECALLING_DIR=$(yq e '.paths.basecalled_output_dir' "${CONFIG_FILE}")
BASECALLED_FILE=$(yq e '.paths.unaligned_bam_name' "${CONFIG_FILE}")
MODEL_DIR="models/"

check_vars \
  "DORADO_MODEL_NAME" \
  "MODEL_SPEED" \
  "BASECALLING_MODIFICATIONS" \
  "BATCHSIZE" \
  "POD5_DIR" \
  "POD5_NAME" \
  "BASECALLING_DIR" \
  "BASECALLED_FILE"

# --- Main script ---
log_info "Using Dorado executable at ${DORADO_EXECUTABLE}" >&2

download_dorado_model() {
  #dna_r9.4.1_e8_fast@v3.4_5mCG_5hmCG@v0
  log_info "Downloading Dorado model ${DORADO_MODEL_NAME} into ${MODEL_DIR}"
  "${DORADO_EXECUTABLE}" download --model dna_r9.4.1_e8_fast@v3.4 --models-directory "${MODEL_DIR}"
  #"${DORADO_EXECUTABLE}" download --model dna_r9.4.1_e8_fast@v3.4_5mCG_5hmCG@v0 --models-directory "${MODEL_DIR}"
}

basecalling_pod5() {
  local POD5_INPUT="${POD5_DIR}${POD5_NAME}"
  local BASECALLED_OUTPUT="${BASECALLING_DIR}${BASECALLED_FILE}"

  #  download_dorado_model
  # ---- You could also make it so that the alignment and basecalling is done at once. That should be the next thing to test.

  local dorado_cmd=(
    "${DORADO_EXECUTABLE}"
    basecaller
    "${MODEL_SPEED},${BASECALLING_MODIFICATIONS}"
    "${POD5_INPUT}"
    --batchsize "${BATCHSIZE}"
  )

  log_info "Executing command: ${dorado_cmd[*]} > ${BASECALLED_OUTPUT}"

  "${dorado_cmd[@]}" > "${BASECALLED_OUTPUT}"



  #"${DORADO_EXECUTABLE}" basecaller "${MODEL_SPEED},${BASECALLING_MODIFICATIONS}" --modified-bases "${BASECALLING_MODIFICATIONS}" --batchsize "${BATCHSIZE}" data/pod5_output/all_reads.pod5 > data/basecalled_output/newcalls.bam
  #"${DORADO_EXECUTABLE}" basecaller "${MODEL_SPEED},${BASECALLING_MODIFICATIONS}" "${POD5_INPUT}" --batchsize "${BATCHSIZE}"  > data/basecalled_output/newcalls2.bam
  #"${DORADO_EXECUTABLE}" basecaller fast --mpd5mCG_5hmCG --batchsize 64 data/pod5_output/all_reads.pod5 > data/basecalled_output/newcalls.bam

  # Output bam file can be checked out with samtools - if we want to see whether there's any methylation in the bam file use
  # samtools view path/to/file | grep 'MM:Z:' | head -n 5

}


#download_dorado_model
basecalling_pod5