#!/bin/bash

# scripts/utils/logging.sh

# ANSI colour codes
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

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
