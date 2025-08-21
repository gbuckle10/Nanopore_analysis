#!/bin/bash

# --- Config ---



# --- Dorado installation ---

install_dorado() {
  local VERSION=$1
  echo "--- Downloading and setting up Dorado ---"
  echo "--- Checking for existing Dorado installation ---"
  if [ -d "tools/dorado_v$VERSION" ]; then
    echo -e "Dorado folder already exists"
    return 0
  fi

  echo "Dorado not found, so I'll install it now."

  echo "Creating installation directory in tools"
  mkdir -p "tools"

  echo "Downloading dorado from https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-win64.zip"
  curl -o tools/dorado.zip  "https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-win64.zip"

  echo "Downloading dorado version 0.9.6 from https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.6-win64.zip"
  echo "Download complete. Unzipping and setting up."

  unzip -q tools/dorado.zip -d tools
  mv "tools/dorado-$VERSION-win64"* "tools/dorado_v$VERSION"
  rm tools/dorado.zip

  echo "-------- Dorado was successfully installed ---------------"

  echo "-------- Setting Dorado exe to path --------------"
  LOCAL_DORADO_BIN="tools/dorado_v$VERSION/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"


  echo "Installing pod5 package"
  pip install pod5
}

download_fast5_data() {
  alias mkdir='/usr/bin/mkdir'
  set -e

  echo "----Starting Data Download Script----"

  S3_BASE_PATH="s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"

  echo "--- Making fast5 input and basecalled output folders ---"
  mkdir -p data
  mkdir -p data/fast5_input
  mkdir -p data/basecalled_output

  echo "--- Downloading 20 fast5 files ---"
  #aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request | head -n 20 | awk '{print $4}' | xargs -I {} aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/{} data/fast5_input/ --no-sign-request

  aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request \
  | head -n 20 \
  | awk '{print $4}' \
  | while read -r filename; do
      aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"$filename" data/fast5_input/"$filename" --no-sign-request
  done
}

convert_fast5_to_pod5() {
  echo "Converting fast5 to pod5"

  mkdir data/pod5_output/all_reads.pod5

  pod5 convert fast5 data/fast5_input --output data/pod5_output/all_reads.pod5
}

basecalling_pod5() {
  # This line doesn't work on the old fast5 data I found online...it will work for new ones but the data I downloaded is too old for this version of dorado
  dorado basecaller hac data/pod5_output/all_reads.pod5/output.pod5 > data/basecalled_output/calls.bam
}

run_script(){
  DORADO_VERSION="1.1.1"
  install_dorado "$DORADO_VERSION"
}

run_script