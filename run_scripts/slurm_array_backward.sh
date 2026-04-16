#!/bin/bash
#SBATCH -J run_flexpart_backward
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0:30:00
#SBATCH -p fdr
#SBATCH --mem=10G
#SBATCH -a 1-1%1
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err

set -euo pipefail

# Array note:
# - This line is only a fallback for direct sbatch usage.
# - Normal usage should go through submit_slurm_array_backward.sh, which passes
#   a dynamic --array value computed from START_DATE/END_DATE.

# Required environment variables (set by submit helper or sbatch --export):
#   START_DATE: YYYYMMDDHH (inclusive)
#   END_DATE: YYYYMMDDHH (inclusive)
# Optional:
#   DOMAIN: domain name for run_backward_batch.py (default: EASTASIA)
#   RECEPTOR: site code (default: GSN)
#   STEP_HOURS: spacing between array timestamps (default: 1)
#   BACKWARD_DAYS: backward duration in days (default: 20)
#   NUM_PARTICLES: number of particles (default: 20000)
#   LSUBGRID: set COMMAND LSUBGRID (default: 0; use 1 for FLEXINVERT-style behavior)
#   OUTROOT: directory for per-task run folders (default: /net/fs06/d2/$USER/flexpart_outs)
#   POSTPROCESS_LOWEST_MAGL: low-level footprint cutoff (default: 100)
#   POSTPROCESS_SOURCE_LAYER_THICKNESS_M: conversion thickness for SRR units (default: 100)
#   DISABLE_AUTO_POSTPROCESS: set to 1 to skip postprocess in this job (default: 0)
#   PRUNE_TO_GRID_FILES: set to 1 to keep only output/grid_time_*.nc after FLEXPART (default: 0)
#   FLEXPART_REPO_ROOT: absolute path to repository root (recommended; auto-set by submit helper)
#   PYTHON_CMD: Python 3 executable to use when no venv is found (default: python3)
#   POSTPROCESS_PYTHON_CMD: Python executable for postprocessing (default: PYTHON_CMD)
#   USE_PROJECT_VENV: set to 1 to activate <repo>/.venv after conda (default: 0)
#   FLEXPART_EXE: optional absolute path to FLEXPART executable to use in runs
#   DEBUG_ENV: set to 1 to print environment/linker diagnostics

DOMAIN="${DOMAIN:-EASTASIA}"
RECEPTOR="${RECEPTOR:-GSN}"
STEP_HOURS="${STEP_HOURS:-1}"
BACKWARD_DAYS="${BACKWARD_DAYS:-20}"
NUM_PARTICLES="${NUM_PARTICLES:-20000}"
LSUBGRID="${LSUBGRID:-0}"
POSTPROCESS_LOWEST_MAGL="${POSTPROCESS_LOWEST_MAGL:-100}"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="${POSTPROCESS_SOURCE_LAYER_THICKNESS_M:-100}"
DISABLE_AUTO_POSTPROCESS="${DISABLE_AUTO_POSTPROCESS:-0}"
PRUNE_TO_GRID_FILES="${PRUNE_TO_GRID_FILES:-0}"
PYTHON_CMD="${PYTHON_CMD:-}"
POSTPROCESS_PYTHON_CMD="${POSTPROCESS_PYTHON_CMD:-}"
USE_PROJECT_VENV="${USE_PROJECT_VENV:-0}"
FLEXPART_EXE="${FLEXPART_EXE:-}"
DEBUG_ENV="${DEBUG_ENV:-0}"

if [[ -z "${START_DATE:-}" || -z "${END_DATE:-}" ]]; then
  echo "ERROR: START_DATE and END_DATE must be set (YYYYMMDDHH)."
  exit 2
fi

if [[ ! "${START_DATE}" =~ ^[0-9]{10}$ || ! "${END_DATE}" =~ ^[0-9]{10}$ ]]; then
  echo "ERROR: START_DATE/END_DATE must match YYYYMMDDHH."
  exit 2
fi

if ! [[ "${STEP_HOURS}" =~ ^[0-9]+$ ]] || [[ "${STEP_HOURS}" -lt 1 ]]; then
  echo "ERROR: STEP_HOURS must be an integer >= 1."
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${FLEXPART_REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
OUTROOT="${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs}"

if [[ ! -f "${REPO_ROOT}/run_scripts/run_backward_batch.py" ]]; then
  echo "ERROR: run_backward_batch.py not found under REPO_ROOT=${REPO_ROOT}"
  echo "Hint: submit through run_scripts/submit_slurm_array_backward.sh so FLEXPART_REPO_ROOT is exported."
  exit 2
fi

# Increase stack limit to prevent segfaults in deep simulation loops
# (FLEXPART has large temporary arrays in nested subroutine calls)
ulimit -s unlimited
echo "Stack size (soft): $(ulimit -s)"
echo "Stack size (hard): $(ulimit -Hs)"

# System-specific environment setup (Svante):
#   source ~/.bashrc
#   initialize conda shell hook
#   conda activate flexpart
if [[ -f "${HOME}/.bashrc" ]]; then
  # .bashrc (and /etc/bashrc) may reference PS1 and other interactive vars.
  # Temporarily disable nounset to avoid failures in non-interactive Slurm shells.
  set +u
  # shellcheck disable=SC1090
  source "${HOME}/.bashrc"
  set -u
fi

if command -v conda >/dev/null 2>&1; then
  # In non-interactive shells, conda activation requires the shell hook.
  # Some conda activate.d scripts assume unset vars are allowed; run activation
  # with nounset disabled and restore strict mode immediately afterwards.
  set +u
  conda_activate_ok=1
  # shellcheck disable=SC1090
  eval "$(conda shell.bash hook)" || conda_activate_ok=0
  if [[ "${conda_activate_ok}" -eq 1 ]]; then
    conda activate flexpart || conda_activate_ok=0
  fi
  set -u
  if [[ "${conda_activate_ok}" -ne 1 ]]; then
    echo "ERROR: failed to run 'conda activate flexpart'."
    exit 2
  fi
else
  echo "WARNING: conda not found in PATH after sourcing ~/.bashrc; continuing without conda activation."
fi

start_epoch="$(date -u -d "${START_DATE:0:8} ${START_DATE:8:2}:00:00" +%s)"
end_epoch="$(date -u -d "${END_DATE:0:8} ${END_DATE:8:2}:00:00" +%s)"
step_seconds="$(( STEP_HOURS * 3600 ))"

if [[ "${end_epoch}" -lt "${start_epoch}" ]]; then
  echo "ERROR: END_DATE is earlier than START_DATE."
  exit 2
fi

span_seconds="$(( end_epoch - start_epoch ))"
if (( span_seconds % step_seconds != 0 )); then
  echo "ERROR: date range is not divisible by STEP_HOURS=${STEP_HOURS}."
  exit 2
fi

total_tasks="$(( span_seconds / step_seconds + 1 ))"
idx="${SLURM_ARRAY_TASK_ID:-1}"

if (( idx < 1 || idx > total_tasks )); then
  echo "INFO: SLURM_ARRAY_TASK_ID=${idx} outside computed range 1-${total_tasks}; skipping."
  exit 0
fi

target_epoch="$(( start_epoch + (idx - 1) * step_seconds ))"
END_TIME="$(date -u -d "@${target_epoch}" +%Y%m%d%H)"

RUN_DIR="${OUTROOT}/${DOMAIN}_${RECEPTOR}_${END_TIME}"
mkdir -p "${OUTROOT}"

cd "${REPO_ROOT}"
if [[ "${USE_PROJECT_VENV}" == "1" ]]; then
  if [[ -f "${REPO_ROOT}/.venv/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${REPO_ROOT}/.venv/bin/activate"
  else
    echo "WARNING: USE_PROJECT_VENV=1 but ${REPO_ROOT}/.venv/bin/activate not found."
  fi
fi

pick_python_cmd() {
  local candidate
  for candidate in "$@"; do
    [[ -z "${candidate}" ]] && continue
    if [[ "${candidate}" == /* ]]; then
      [[ -x "${candidate}" ]] || continue
    else
      command -v "${candidate}" >/dev/null 2>&1 || continue
    fi
    if "${candidate}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 6) else 1)' >/dev/null 2>&1; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

PYTHON_CMD="$(pick_python_cmd "${PYTHON_CMD}" python python3 "${REPO_ROOT}/.venv/bin/python" || true)"
if [[ -z "${PYTHON_CMD}" ]]; then
  echo "ERROR: no Python >= 3.6 interpreter found."
  echo "Hint: conda env activation may have failed, or set PYTHON_CMD to a compatible interpreter."
  exit 2
fi

if [[ "${DEBUG_ENV}" == "1" ]]; then
  echo "=== Debug Info ==="
  command -v conda || true
  conda info --envs 2>/dev/null | grep flexpart || true
  echo "CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
  echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<empty>}"
  if [[ -n "${FLEXPART_EXE}" && -x "${FLEXPART_EXE}" ]]; then
    ldd "${FLEXPART_EXE}" 2>&1 | head -n 40 || true
  fi
  echo "==================="
fi

echo "Python executable: $("${PYTHON_CMD}" -c 'import sys; print(sys.executable)')"
echo "Python version: $("${PYTHON_CMD}" -c 'import sys; print(sys.version.split()[0])')"
if [[ -n "${POSTPROCESS_PYTHON_CMD}" ]]; then
  echo "Postprocess Python: ${POSTPROCESS_PYTHON_CMD}"
fi

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task=${idx}/${total_tasks} end_time=${END_TIME} receptor=${RECEPTOR} domain=${DOMAIN}"

cmd=("${PYTHON_CMD}" "${REPO_ROOT}/run_scripts/run_backward_batch.py" \
  --domain "${DOMAIN}" \
  --receptor "${RECEPTOR}" \
  --end-time "${END_TIME}" \
  --days "${BACKWARD_DAYS}" \
  --num-particles "${NUM_PARTICLES}" \
  --lsubgrid "${LSUBGRID}" \
  --outdir "${RUN_DIR}" \
  --postprocess-lowest-magl "${POSTPROCESS_LOWEST_MAGL}" \
  --postprocess-source-layer-thickness-m "${POSTPROCESS_SOURCE_LAYER_THICKNESS_M}")

if [[ "${DISABLE_AUTO_POSTPROCESS}" == "1" ]]; then
  cmd+=(--no-postprocess)
  echo "Automatic postprocess disabled for this FLEXPART job (DISABLE_AUTO_POSTPROCESS=1)."
fi

if [[ -n "${FLEXPART_EXE}" ]]; then
  cmd+=(--executable "${FLEXPART_EXE}")
fi

if [[ -n "${POSTPROCESS_PYTHON_CMD}" ]]; then
  cmd+=(--postprocess-python "${POSTPROCESS_PYTHON_CMD}")
fi

"${cmd[@]}"

# Keep only the final postprocessed footprint netCDF for this release hour.
OUTPUT_DIR="${RUN_DIR}/output"
if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${PRUNE_TO_GRID_FILES}" == "1" ]]; then
    shopt -s nullglob
    keep_files=("${OUTPUT_DIR}"/grid_time_*.nc)
    if (( ${#keep_files[@]} == 0 )); then
      echo "WARNING: PRUNE_TO_GRID_FILES=1 but no grid_time_*.nc found in ${OUTPUT_DIR}."
    else
      # Prune everything in output/ except grid_time_*.nc
      for p in "${OUTPUT_DIR}"/*; do
        keep=false
        for k in "${keep_files[@]}"; do
          if [[ "${p}" == "${k}" ]]; then
            keep=true
            break
          fi
        done
        if [[ "${keep}" == false ]]; then
          rm -rf "${p}"
        fi
      done
      
      # Delete all input/config files from RUN_DIR (keep only output/ subdirectory)
      for p in "${RUN_DIR}"/*; do
        if [[ "${p}" != "${OUTPUT_DIR}" ]]; then
          rm -rf "${p}"
        fi
      done
      
      echo "Pruned: kept ${#keep_files[@]} grid_time_*.nc file(s); deleted all other run-dir contents (AGECLASSES, AVAILABLE, COMMAND, IGBP_int1.dat, SPECIES, oh_fields, etc.)"
    fi
    shopt -u nullglob
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] completed end_time=${END_TIME}"
    exit 0
  fi

  shopt -s nullglob
  keep_files=("${OUTPUT_DIR}"/*_FLEXPART_GFS_"${DOMAIN}"_inert_"${END_TIME}".nc)
  if (( ${#keep_files[@]} == 0 )); then
    echo "WARNING: no final footprint file found in ${OUTPUT_DIR} for ${END_TIME}."
  else
    for p in "${OUTPUT_DIR}"/*; do
      keep=false
      for k in "${keep_files[@]}"; do
        if [[ "${p}" == "${k}" ]]; then
          keep=true
          break
        fi
      done
      if [[ "${keep}" == false ]]; then
        rm -rf "${p}"
      fi
    done
  fi
  shopt -u nullglob
fi

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] completed end_time=${END_TIME}"
