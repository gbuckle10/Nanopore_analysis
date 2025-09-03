#!/bin/bash


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

  dorado_folder=tools/dorado-${VERSION}-${os_type}
  echo "--- Checking for existing Dorado installation in ${dorado_folder}---"
  if [ -d "${dorado_folder}" ]; then
    echo "Dorado folder already exists"
    LOCAL_DORADO_BIN="${dorado_folder}/bin"
    export PATH="$LOCAL_DORADO_BIN:$PATH"
    echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"
    return 0
  fi

  echo "Dorado folder ${dorado_folder} not found, so I'll install it now."

  echo "Creating installation directory in tools"
  mkdir -p "tools"

  # Include a check here for whether the archive_filename is there to save more time.
  echo "Downloading dorado version $VERSION for ${os_type} from ${download_url}"
  curl -o tools/${archive_filename}  "$download_url"

  echo "Download complete. Unzipping and setting up."

  case ${os_type} in
    "linux-x64")
      tar -xzvf "tools/${archive_filename}" -C "tools/"
      ;;
    "win64")
      unzip -q "tools/${archive_filename}" -d "tools/"
      ;;
  esac

  folder_to_remove=tools/${archive_filename}
  echo "Removing ${folder_to_remove}"

  rm "${folder_to_remove}"

  echo "-------- Dorado was successfully installed ---------------"

  echo "-------- Setting Dorado exe to path --------------"
  LOCAL_DORADO_BIN="${dorado_folder}/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"

}

install_dorado "0.9.6"