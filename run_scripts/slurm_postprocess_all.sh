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
#   FINAL_DIR: final directory for postprocessed hourly netCDF files (default: OUTROOT/<receptor-lowercase>_hourly when RECEPTOR is set)
#   POSTPROCESS_DRIVER_PYTHON: python for postprocess_all_outputs.py (default: python3)
#   POSTPROCESS_PYTHON_CMD: python for postprocess_footprint.py subprocess (default: /home/$USER/.conda/envs/flexpart-post/bin/python)
#   POSTPROCESS_LOWEST_MAGL: low-level footprint cutoff (default: 100)
#   POSTPROCESS_SOURCE_LAYER_THICKNESS_M: conversion thickness for SRR units (default: 100)
#   POSTPROCESS_OVERWRITE: set to 1 to overwrite existing final files (default: 0)
#   KEEP_RUN_DIRS: set to 1 to keep per-run directories after successful move (default: 0)
#   LIMIT: optional max number of grid_time files to process (default: 0 = all)
#   WRITE_MONTHLY: set to 1 to write monthly aggregated NetCDF files (default: 1)
#   MONTHLY_DIR: directory for monthly NetCDF outputs (default: OUTROOT/<receptor-lowercase> when RECEPTOR is set)
#   RECEPTOR: optional site code used to derive MONTHLY_DIR when unset
#   DRY_RUN: set to 1 for preview only (default: 0)
#   FLEXPART_REPO_ROOT: absolute path to repository root (auto-resolved by launcher)

OUTROOT="${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs}"
if [[ -n "${FINAL_DIR:-}" ]]; then
  FINAL_DIR="${FINAL_DIR}"
elif [[ -n "${RECEPTOR:-}" ]]; then
  FINAL_DIR="${OUTROOT}/$(printf '%s' "${RECEPTOR}" | tr '[:upper:]' '[:lower:]')_hourly"
else
  FINAL_DIR="${OUTROOT}"
fi
POSTPROCESS_DRIVER_PYTHON="${POSTPROCESS_DRIVER_PYTHON:-python3}"
POSTPROCESS_PYTHON_CMD="${POSTPROCESS_PYTHON_CMD:-/home/${USER}/.conda/envs/flexpart-post/bin/python}"
POSTPROCESS_LOWEST_MAGL="${POSTPROCESS_LOWEST_MAGL:-100}"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="${POSTPROCESS_SOURCE_LAYER_THICKNESS_M:-100}"
POSTPROCESS_OVERWRITE="${POSTPROCESS_OVERWRITE:-0}"
KEEP_RUN_DIRS="${KEEP_RUN_DIRS:-0}"
LIMIT="${LIMIT:-0}"
WRITE_MONTHLY="${WRITE_MONTHLY:-1}"
if [[ -n "${MONTHLY_DIR:-}" ]]; then
  MONTHLY_DIR="${MONTHLY_DIR}"
elif [[ -n "${RECEPTOR:-}" ]]; then
  MONTHLY_DIR="${OUTROOT}/$(printf '%s' "${RECEPTOR}" | tr '[:upper:]' '[:lower:]')"
else
  MONTHLY_DIR="${OUTROOT}"
fi
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

# xarray/pandas/numpy stacks on this cluster are not stable on Python 3.13+ yet.
if ! "${POSTPROCESS_PYTHON_CMD}" -c 'import sys; raise SystemExit(0 if sys.version_info < (3, 13) else 1)' >/dev/null 2>&1; then
  echo "ERROR: POSTPROCESS_PYTHON_CMD must be Python <= 3.12 to avoid numpy/pandas segfaults on this system."
  echo "       Current: $("${POSTPROCESS_PYTHON_CMD}" -c 'import sys; print(sys.version.split()[0])') (${POSTPROCESS_PYTHON_CMD})"
  echo "       Recommended: /home/${USER}/.conda/envs/flexpart-post/bin/python"
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
if [[ "${WRITE_MONTHLY}" == "1" ]]; then
  cmd+=(--write-monthly --monthly-dir "${MONTHLY_DIR}")
fi

"${cmd[@]}"