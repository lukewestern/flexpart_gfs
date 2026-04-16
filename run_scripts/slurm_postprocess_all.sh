#!/bin/bash
#SBATCH -J flexpart_postprocess_all
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p fdr
#SBATCH --mem=8G
#SBATCH -o slurm-postprocess-all-%j.out
#SBATCH -e slurm-postprocess-all-%j.err

set -euo pipefail

# Optional environment variables:
#   OUTROOT: root directory containing FLEXPART run folders (default: /net/fs06/d2/$USER/flexpart_outs)
#   FINAL_DIR: final directory for postprocessed netCDF files (default: OUTROOT)
#   POSTPROCESS_DRIVER_PYTHON: python for postprocess_all_outputs.py (default: python3)
#   POSTPROCESS_PYTHON_CMD: python for postprocess_footprint.py subprocess (default: /home/$USER/.conda/envs/flexpart/bin/python)
#   POSTPROCESS_LOWEST_MAGL: low-level footprint cutoff (default: 100)
#   POSTPROCESS_SOURCE_LAYER_THICKNESS_M: conversion thickness for SRR units (default: 100)
#   POSTPROCESS_OVERWRITE: set to 1 to overwrite existing final files (default: 0)
#   KEEP_RUN_DIRS: set to 1 to keep per-run directories after successful move (default: 0)
#   LIMIT: optional max number of grid_time files to process (default: 0 = all)
#   DRY_RUN: set to 1 for preview only (default: 0)
#   FLEXPART_REPO_ROOT: absolute path to repository root (auto-resolved by launcher)

OUTROOT="${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs}"
FINAL_DIR="${FINAL_DIR:-${OUTROOT}}"
POSTPROCESS_DRIVER_PYTHON="${POSTPROCESS_DRIVER_PYTHON:-python3}"
POSTPROCESS_PYTHON_CMD="${POSTPROCESS_PYTHON_CMD:-/home/${USER}/.conda/envs/flexpart/bin/python}"
POSTPROCESS_LOWEST_MAGL="${POSTPROCESS_LOWEST_MAGL:-100}"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="${POSTPROCESS_SOURCE_LAYER_THICKNESS_M:-100}"
POSTPROCESS_OVERWRITE="${POSTPROCESS_OVERWRITE:-0}"
KEEP_RUN_DIRS="${KEEP_RUN_DIRS:-0}"
LIMIT="${LIMIT:-0}"
DRY_RUN="${DRY_RUN:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${FLEXPART_REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
POSTPROCESS_ALL_SCRIPT="${REPO_ROOT}/run_scripts/postprocess_all_outputs.py"

if [[ ! -f "${POSTPROCESS_ALL_SCRIPT}" ]]; then
  echo "ERROR: postprocess_all_outputs.py not found: ${POSTPROCESS_ALL_SCRIPT}"
  exit 2
fi

if ! command -v "${POSTPROCESS_DRIVER_PYTHON}" >/dev/null 2>&1; then
  echo "ERROR: POSTPROCESS_DRIVER_PYTHON not found: ${POSTPROCESS_DRIVER_PYTHON}"
  exit 2
fi

if [[ ! -x "${POSTPROCESS_PYTHON_CMD}" ]]; then
  echo "ERROR: POSTPROCESS_PYTHON_CMD is not executable: ${POSTPROCESS_PYTHON_CMD}"
  exit 2
fi

echo "Root dir: ${OUTROOT}"
echo "Final dir: ${FINAL_DIR}"
echo "Driver python: $("${POSTPROCESS_DRIVER_PYTHON}" -c 'import sys; print(sys.executable)')"
echo "Postprocess python: ${POSTPROCESS_PYTHON_CMD}"

cmd=("${POSTPROCESS_DRIVER_PYTHON}" "${POSTPROCESS_ALL_SCRIPT}"
  --root-dir "${OUTROOT}"
  --final-dir "${FINAL_DIR}"
  --python "${POSTPROCESS_PYTHON_CMD}"
  --postprocess-lowest-magl "${POSTPROCESS_LOWEST_MAGL}"
  --postprocess-source-layer-thickness-m "${POSTPROCESS_SOURCE_LAYER_THICKNESS_M}")

if [[ "${POSTPROCESS_OVERWRITE}" == "1" ]]; then
  cmd+=(--overwrite)
fi
if [[ "${KEEP_RUN_DIRS}" == "1" ]]; then
  cmd+=(--keep-run-dirs)
fi
if [[ "${DRY_RUN}" == "1" ]]; then
  cmd+=(--dry-run)
fi
if [[ "${LIMIT}" =~ ^[0-9]+$ ]] && [[ "${LIMIT}" -gt 0 ]]; then
  cmd+=(--limit "${LIMIT}")
fi

"${cmd[@]}"