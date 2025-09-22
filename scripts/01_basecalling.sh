#!/bin/bash

set -e

CONFIG_FILE=$1

if [ ! -f "scripts/runtime_config.sh" ]; then
  echo "ERROR: runtime_config.sh not found. Please run 'setup' step first..."
  exit 1
fi
source "scripts/runtime_config.sh"

DORADO_MODEL_NAME=$(yq e '.parameters.basecalling.base_model_name' "${CONFIG_FILE}")
MODEL_SPEED=$(yq e '.parameters.basecalling.model_speed' "${CONFIG_FILE}")
BASECALLING_MODIFICATIONS=$(yq e '.parameters.basecalling.basecalling_modifications' "${CONFIG_FILE}")
BATCHSIZE=$(yq e '.parameters.analysis.batch_size' "${CONFIG_FILE}")
MODEL_DIR="models/"

# --- Main script ---
echo "INFO: Using Dorado executable at ${DORADO_EXECUTABLE}" >&2


check_device_for_dorado() {
  if command -v nvidia-smi &> /dev/null; then
    echo "NVIDIA GPU found, so we'll use that for basecalling"

    # Tell dorado to tun with --device flag cuda:0
    export DORADO_DEVICE_FLAG="--device cuda:0"
    # Default batchsize is 128
    export DORADO_BATCHSIZE_FLAG="--batchsize 128"

  else
    echo "No NVIDIA GPU found, so we'll run with the CPU only"

    # Tell dorado to run with --device flag cpu
    export DORADO_DEVICE_FLAG="--device cpu"
    # Default batchsize is 128
    export DORADO_BATCHSIZE_FLAG=""
  fi

  echo "Dorado device flags set: ${DORADO_DEVICE_FLAG} ${DORADO_BATCHSIZE_FLAG}"
}

download_dorado_model() {
  #dna_r9.4.1_e8_fast@v3.4_5mCG_5hmCG@v0
  echo "Downloading Dorado model"
  "${DORADO_EXECUTABLE}" download --model dna_r9.4.1_e8_fast@v3.4 --models-directory "${MODEL_DIR}"
  "${DORADO_EXECUTABLE}" download --model dna_r9.4.1_e8_fast@v3.4_5mCG_5hmCG@v0 --models-directory "${MODEL_DIR}"

}

basecalling_pod5() {
  echo "Adding dorado to the path"
  local POD5_INPUT="data/pod5_output/all_reads.pod5"
  local BASECALLED_OUTPUT="data/basecalled_output/calls.bam"

  # ---- You need to add the download_dorado_model here
  #  download_dorado_model
  # ---- You could also make it so that the alignment and basecalling is done at once. That should be the next thing to test.

  #"${DORADO_EXECUTABLE}" basecaller "${MODEL_SPEED},${BASECALLING_MODIFICATIONS}" --modified-bases "${BASECALLING_MODIFICATIONS}" --batchsize "${BATCHSIZE}" data/pod5_output/all_reads.pod5 > data/basecalled_output/newcalls.bam
  "${DORADO_EXECUTABLE}" basecaller fast,5mCG_5hmCG data/pod5_output/all_reads.pod5 > data/basecalled_output/newcalls.bam
  #"${DORADO_EXECUTABLE}" basecaller fast --mpd5mCG_5hmCG --batchsize 64 data/pod5_output/all_reads.pod5 > data/basecalled_output/newcalls.bam

  # Output bam file can be checked out with samtools - if we want to see whether there's any methylation in the bam file use
  # samtools view path/to/file | grep 'Mm:Z:' | head -n 5

}


download_dorado_model
basecalling_pod5
