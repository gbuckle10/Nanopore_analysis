#!/bin/bash

# scripts/utils/logging.sh

# Print a formatted log message to stderr
log_info() {
  #local message="$1"
  #echo -e "$(date '+%Y-%m-%d %H:%M:%S') - ${GREEN}INFO${NC}     - ${message}" >&2
  echo "[Bash] INFO: $1" >&2
}

log_warn() {
  echo "[Bash] WARNING: $1" >&2
}

log_error() {
  echo "[Bash] ERROR: $1" >&2
}

handle_error() {
  local exit_code=$?
  local line_number=$1
  local failed_command="${BASH_COMMAND}"

  log_error "Command failed on line ${line_number} with exit code ${exit_code}"
  log_error "Command that failed: '${failed_command}'"
}
trap 'handle_error $LINENO' ERR

check_vars() {
  for var in "$@"; do
    if [[ -z "${!var}" ]]; then
      log_error "Required variable '${var_name}' is not set or is empty. Please check config.yaml"
      exit 1
    fi
  done
}