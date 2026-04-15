#!/usr/bin/env bash
# Wrapper script for download_gfs_archive_python.py
# Provides consistent interface with other download tools

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/download_gfs_archive_python.py"

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Error: $PYTHON_SCRIPT not found" >&2
  exit 1
fi

python3 "$PYTHON_SCRIPT" "$@"
