#!/bin/bash
#SBATCH -J flexpart_postprocess_monthly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH -t 24:00:00
#SBATCH -p fdr
#SBATCH --mem=8G
#SBATCH -o slurm-postprocess-monthly-%A_%a.out
#SBATCH -e slurm-postprocess-monthly-%A_%a.err

set -euo pipefail

# Required environment variables:
#   MONTH_LIST: whitespace-separated YYYYMM values (comma-separated is also accepted)
# Optional environment variables:
#   OUTROOT: root directory containing FLEXPART run folders
#   FINAL_DIR: final directory for postprocessed hourly netCDF files
#   POSTPROCESS_DRIVER_PYTHON: python for postprocess_all_outputs.py
#   POSTPROCESS_PYTHON_CMD: python for postprocess_footprint.py subprocess
#   POSTPROCESS_FOOTPRINT_OUTHEIGHT_M: cumulative top outheight for footprint integration
#   POSTPROCESS_LOWEST_MAGL: deprecated alias for POSTPROCESS_FOOTPRINT_OUTHEIGHT_M
#   POSTPROCESS_SOURCE_LAYER_THICKNESS_M: conversion thickness for SRR units
#   POSTPROCESS_WORKERS: concurrent NetCDF postprocess workers within each month task
#   POSTPROCESS_OVERWRITE: set to 1 to overwrite existing final files
#   KEEP_RUN_DIRS: set to 1 to keep per-run directories after successful move
#   LIMIT: optional max number of grid_time files to process (0 = all)
#   WRITE_MONTHLY: set to 1 to write monthly aggregated NetCDF files
#   MONTHLY_DIR: directory for monthly NetCDF outputs
#   RECEPTOR: optional site code used to derive output directories
#   DRY_RUN: set to 1 for preview only
#   FLEXPART_REPO_ROOT: absolute path to repository root (auto-resolved by launcher)

if [[ -z "${MONTH_LIST:-}" ]]; then
  echo "ERROR: MONTH_LIST is required for monthly array jobs."
  exit 2
fi

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
POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_FOOTPRINT_OUTHEIGHT_M:-100}"
if [[ -n "${POSTPROCESS_LOWEST_MAGL:-}" ]]; then
  POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_LOWEST_MAGL}"
fi
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="${POSTPROCESS_SOURCE_LAYER_THICKNESS_M:-100}"
POSTPROCESS_WORKERS="${POSTPROCESS_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
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

if [[ "${POSTPROCESS_PYTHON_CMD}" == *"/.conda/envs/flexpart-post/"* ]]; then
  offending_pypi_lines=""
  if command -v conda >/dev/null 2>&1; then
    offending_pypi_lines="$(conda list -n flexpart-post 2>/dev/null | grep -E '^(numpy|pandas|xarray|netcdf4)\s+.+\s+pypi_0\s+pypi$' || true)"
  fi
  if [[ -n "${offending_pypi_lines}" ]]; then
    echo "ERROR: flexpart-post scientific stack was installed from pip (pypi), which is incompatible on older compute-node glibc."
    echo "       Offending package rows from 'conda list -n flexpart-post':"
    while IFS= read -r row; do
      echo "         ${row}"
    done <<< "${offending_pypi_lines}"
    echo "       Repair with conda packages before running:"
    echo "         conda install -n flexpart-post -c conda-forge 'numpy<2.0' 'pandas<2.2' 'xarray<2025' netcdf4"
    exit 2
  fi
fi

if ! "${POSTPROCESS_PYTHON_CMD}" - <<'PY' >/dev/null 2>&1
import importlib.util
import sys

required = ["numpy", "xarray", "pandas", "netCDF4"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    print(
        "Missing Python package(s) in POSTPROCESS_PYTHON_CMD environment: "
        + ", ".join(missing),
        file=sys.stderr,
    )
    raise SystemExit(1)
PY
then
  echo "ERROR: POSTPROCESS_PYTHON_CMD is missing required packages (numpy/xarray/pandas/netCDF4)."
  echo "       Interpreter: ${POSTPROCESS_PYTHON_CMD}"
  echo "       Install with one of:"
  echo "         ${POSTPROCESS_PYTHON_CMD} -m pip install numpy xarray pandas netCDF4"
  echo "         conda install -n flexpart-post numpy xarray pandas netcdf4"
  exit 2
fi

if ! "${POSTPROCESS_PYTHON_CMD}" -c 'import sys; raise SystemExit(0 if sys.version_info < (3, 13) else 1)' >/dev/null 2>&1; then
  echo "ERROR: POSTPROCESS_PYTHON_CMD must be Python <= 3.12 to avoid numpy/pandas segfaults on this system."
  echo "       Current: $("${POSTPROCESS_PYTHON_CMD}" -c 'import sys; print(sys.version.split()[0])') (${POSTPROCESS_PYTHON_CMD})"
  echo "       Recommended: /home/${USER}/.conda/envs/flexpart-post/bin/python"
  exit 2
fi

# Slurm --export uses commas to separate variables, so the launcher exports
# months as whitespace-delimited values. Keep comma parsing for backward compatibility.
months_raw="${MONTH_LIST}"
if [[ "${months_raw}" == *,* ]]; then
  IFS=',' read -r -a months <<< "${months_raw}"
else
  read -r -a months <<< "${months_raw}"
fi
count="${#months[@]}"
idx="${SLURM_ARRAY_TASK_ID:-1}"

if [[ ! "${idx}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID must be numeric; got ${idx}"
  exit 2
fi

if (( idx < 1 || idx > count )); then
  echo "INFO: SLURM_ARRAY_TASK_ID=${idx} outside month range 1-${count}; skipping."
  exit 0
fi

target_month="${months[$((idx - 1))]}"
if [[ ! "${target_month}" =~ ^[0-9]{6}$ ]]; then
  echo "ERROR: invalid month token at index ${idx}: ${target_month}"
  exit 2
fi

echo "Root dir: ${OUTROOT}"
echo "Final dir: ${FINAL_DIR}"
echo "Monthly dir: ${MONTHLY_DIR}"
echo "Target month: ${target_month} (${idx}/${count})"
echo "Driver python: $("${POSTPROCESS_DRIVER_PYTHON}" -c 'import sys; print(sys.executable)')"
echo "Postprocess python: ${POSTPROCESS_PYTHON_CMD}"

cmd=("${POSTPROCESS_DRIVER_PYTHON}" "${POSTPROCESS_ALL_SCRIPT}"
  --root-dir "${OUTROOT}"
  --final-dir "${FINAL_DIR}"
  --python "${POSTPROCESS_PYTHON_CMD}"
  --postprocess-footprint-outheight-m "${POSTPROCESS_FOOTPRINT_OUTHEIGHT_M}"
  --workers "${POSTPROCESS_WORKERS}"
  --postprocess-source-layer-thickness-m "${POSTPROCESS_SOURCE_LAYER_THICKNESS_M}"
  --month "${target_month}")

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
