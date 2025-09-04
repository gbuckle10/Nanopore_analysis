#!/bin/bash

set -e

DORADO_VERSION=$1
FAST5_DOWNLOAD_URL=$2
FAST5_DESTINATION_DIR=$3
NUM_FAST5_FILES=$4
POD5_DIR=$5
SHOULD_DOWNLOAD_FAST5_FILES=$6
SHOULD_CONVERT_FAST5_TO_POD5=$7

RUNTIME_CONFIG="scripts/runtime_config.sh"

make_directories() {
  mkdir -p data
  mkdir -p logs
  mkdir -p models
  mkdir -p reference_genomes
  mkdir -p scripts
  mkdir -p tools
}

install_dorado() {
  echo "--- Downloading and setting up Dorado ---" >&2

  echo "-- Detecting OS --" >&2
  case "$(uname -s)" in
    Linux*)
      os_type="linux-x64"
      archive_filename="dorado.tar.gz"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-${os_type}.tar.gz"
      ;;
    CYGWIN*|MINGW*|MSYS*)
      os_type="win64"
      archive_filename="dorado.zip"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-${os_type}.zip"
      ;;
    *)
      echo "I don't have a clue what operating system you're trying to use..." >&2
      return 1
      ;;
  esac
  echo "The OS has been found, it's $os_type" >&2

  DORADO_DIR=tools/dorado-${DORADO_VERSION}-${os_type}
  echo "--- Checking for existing Dorado installation in ${DORADO_DIR}---"
  if [ -d "${DORADO_DIR}" ]; then
    echo "Dorado folder already exists" >&2

  else 
    echo "Dorado folder ${DORADO_DIR} not found, so I'll install it now." >&2

    # Include a check here for whether the archive_filename is there to save more time.
    echo "Downloading dorado version ${DORADO_VERSION} for ${os_type} from ${download_url}" >&2
    curl -o tools/${archive_filename}  "$download_url"
  
    echo "Download complete. Unzipping and setting up." >&2
  
    case ${os_type} in
      "linux-x64")
        tar -xzvf "tools/${archive_filename}" -C "tools/"
        ;;
      "win64")
        unzip -q "tools/${archive_filename}" -d "tools/"
        ;;
    esac
  
    folder_to_remove=tools/${archive_filename}
    echo "Removing ${folder_to_remove}" >&2
  
    rm "${folder_to_remove}"
    
    
  fi
  
  LOCAL_DORADO_BIN="${DORADO_DIR}/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"

  DORADO_EXECUTABLE_PATH=$(realpath "${LOCAL_DORADO_BIN}/dorado")
  echo "INFO: Saving Dorado path to ${RUNTIME_CONFIG}" >&2
  echo "export DORADO_EXECUTABLE=\"${DORADO_EXECUTABLE_PATH}\"" > "${RUNTIME_CONFIG}"


  echo "-------- Dorado was successfully installed ---------------" >&2
}

download_fast5_data() {
  #alias mkdir='/usr/bin/mkdir'


  echo "----Starting Data Download Script----"


  echo "--- Making fast5 input folder at ${FAST5_DESTINATION_DIR}---"
  mkdir -p "${FAST5_DESTINATION_DIR}"

  echo "--- Downloading ${NUM_FAST5_FILES} fast5 files into '${FAST5_DESTINATION_DIR}' ---"
  echo "--- Downloading from '${FAST5_DOWNLOAD_URL}'"

  aws s3 ls "${FAST5_DOWNLOAD_URL}" --no-sign-request \
  | head -n "${NUM_FAST5_FILES}" \
  | awk '{print $4}' \
  | while read -r filename; do
      aws s3 cp "${FAST5_DOWNLOAD_URL}${filename}" "${FAST5_DESTINATION_DIR}""$filename" --no-sign-request || true
  done

  echo "--- Fast5 input downloaded ---"

}

convert_fast5_to_pod5() {
  echo "--- Converting fast5 to pod5 ---"

  mkdir -p "${POD5_DIR}"

  pod5 convert fast5 "${FAST5_DESTINATION_DIR}" --output "${POD5_DIR}/all_reads.pod5" --force-overwrite
}

download_methylation_atlas_and_illumina_manifest(){
  mkdir -p "data/atlas/"

  if ! [ -f "data/atlas/illumina_manifest.csv" ]; then
    echo "The illumina manifest doesn't already exist so we will download it."
    wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv -O "data/atlas/illumina_manifest.csv"
  fi
  if ! [ -f "data/atlas/full_atlas.csv" ]; then
    echo "The methylation atlas doesn't already exist so we will download it."
    wget https://github.com/nloyfer/meth_atlas/raw/refs/heads/master/full_atlas.csv.gz -O "data/atlas/full_atlas.csv.gz"
    gunzip -v "data/atlas/full_atlas.csv.gz"
    #rm "data/atlas/full_atlas.csv.gz"
  fi

  echo "Downloading Illumina manifest."
  wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv -O "data/atlas/illumina_manifest.csv"
  echo "Illumina manifest successfully downloaded."
}

make_directories
install_dorado
download_methylation_atlas_and_illumina_manifest
if [[ "$SHOULD_DOWNLOAD_FAST5_FILES" == "true" ]]; then
  download_fast5_data
fi
if [[ "$SHOULD_CONVERT_FAST5_TO_POD5" == "true" ]]; then
  convert_fast5_to_pod5
fi

