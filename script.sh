#!/bin/bash

# --- Config ---



# --- Dorado installation ---

install_dorado() {
  local VERSION=$1
  echo "--- Downloading and setting up Dorado ---"

  echo "-- Detecting OS --"
  case "$(uname -s)" in
    Linux*)
      os_type="linux-x64"
      archive_filename="dorado.tar.gz"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-${os_type}.tar.gz"
      ;;
    CYGWIN*|MINGW*|MSYS*)
      os_type="win64"
      archive_filename="dorado.zip"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-${os_type}.zip"
      ;;
    *)
      echo "I don't have a clue what operating system you're trying to use..."
      return 1
      ;;
  esac
  echo "The OS has been found, it's $os_type"

  echo "--- Checking for existing Dorado installation ---"
  expected_folder=tools/dorado-$VERSION-${os_type}
  if [ -d "tools/dorado-$VERSION-${os_type}" ]; then
    echo -e "Dorado folder already exists"
    return 0
  fi

  echo "Dorado folder ${expected_folder} not found, so I'll install it now."

  echo "Creating installation directory in tools"
  mkdir -p "tools"

  # Include a check here for whether the archive_filename is there to save more time.
  echo "Downloading dorado version $VERSION for ${os_type} from ${download_url}"
  curl -o tools/${archive_filename}  $download_url

  echo "Download complete. Unzipping and setting up."

  case ${os_type} in
    "linux-x64")
      tar -xzvf "tools/${archive_filename}" -C "tools/"
      ;;
    "win64")
      unzip -q "tools/${archive_filename}" -d "tools/"
      ;;
  esac

  mv "tools/dorado-$VERSION-${os_type}"* "tools/dorado_v$VERSION-${os_type}"
  rm tools/${archive_filename}

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
  DORADO_VERSION="0.9.6"
  install_dorado "$DORADO_VERSION"
}

run_script