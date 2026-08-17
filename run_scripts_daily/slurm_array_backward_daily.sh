#!/bin/bash
#SBATCH -J flexpart_daily
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -p edr,fdr
#SBATCH --mem=60G
#SBATCH -a 1-1%1
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err

set -euo pipefail

# Daily FLEXPART backward array worker.
#
# Each array task runs all 24 hourly releases for one calendar day in a single
# FLEXPART simulation, producing 24 grid_time_*.nc files.
#
# Required environment variables (set by submit helper or sbatch --export):
#   START_DATE: YYYYMMDD (inclusive, first day)
#   END_DATE:   YYYYMMDD (inclusive, last day)
# Optional:
#   DOMAIN, RECEPTOR, BACKWARD_DAYS, NUM_PARTICLES, IPOUT, NXSHIFT
#   LSUBGRID, LINIT_COND
#   POSTPROCESS_FOOTPRINT_OUTHEIGHT_M, POSTPROCESS_SOURCE_LAYER_THICKNESS_M
#   DISABLE_AUTO_POSTPROCESS: 1 = skip per-run postprocess (default: 0)
#   PRUNE_OUTPUTS: 1 = after postprocessing, delete raw FLEXPART outputs, keep footprints/exit-points (default: 1)
#   SKIP_IF_RUN_DIR_EXISTS: 1 = skip task when ${OUTROOT}/${DOMAIN}_${RECEPTOR}_${DATE} already exists (default: 0)
#   OUTROOT, FLEXPART_EXE, PYTHON_CMD, POSTPROCESS_PYTHON_CMD
#   USE_PROJECT_VENV, DEBUG_ENV, FLEXPART_REPO_ROOT

DOMAIN="${DOMAIN:-EASTASIA}"
RECEPTOR="${RECEPTOR:-GSN}"
BACKWARD_DAYS="${BACKWARD_DAYS:-20}"
NUM_PARTICLES="${NUM_PARTICLES:-20000}"
IPOUT="${IPOUT:-1}"
NXSHIFT="${NXSHIFT:-}"
LSUBGRID="${LSUBGRID:-1}"
LINIT_COND="${LINIT_COND:-1}"
POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_FOOTPRINT_OUTHEIGHT_M:-100}"
if [[ -n "${POSTPROCESS_LOWEST_MAGL:-}" ]]; then
  POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_LOWEST_MAGL}"
fi
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="${POSTPROCESS_SOURCE_LAYER_THICKNESS_M:-100}"
DISABLE_AUTO_POSTPROCESS="${DISABLE_AUTO_POSTPROCESS:-0}"
PRUNE_OUTPUTS="${PRUNE_OUTPUTS:-1}"
SKIP_IF_POSTPROCESSED_EXISTS="${SKIP_IF_POSTPROCESSED_EXISTS:-1}"
SKIP_IF_RUN_DIR_EXISTS="${SKIP_IF_RUN_DIR_EXISTS:-0}"
PYTHON_CMD="${PYTHON_CMD:-}"
POSTPROCESS_PYTHON_CMD="${POSTPROCESS_PYTHON_CMD:-}"
USE_PROJECT_VENV="${USE_PROJECT_VENV:-0}"
FLEXPART_EXE="${FLEXPART_EXE:-}"
DEBUG_ENV="${DEBUG_ENV:-0}"

if [[ -z "${START_DATE:-}" || -z "${END_DATE:-}" ]]; then
  echo "ERROR: START_DATE and END_DATE must be set (YYYYMMDD)."
  exit 2
fi

if [[ ! "${START_DATE}" =~ ^[0-9]{8}$ || ! "${END_DATE}" =~ ^[0-9]{8}$ ]]; then
  echo "ERROR: START_DATE/END_DATE must match YYYYMMDD."
  exit 2
fi

if [[ -n "${IPOUT}" && ! "${IPOUT}" =~ ^[0-2]$ ]]; then
  echo "ERROR: IPOUT must be one of 0, 1, 2 when set."
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${FLEXPART_REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
OUTROOT="${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs}"

if [[ ! -f "${REPO_ROOT}/run_scripts_daily/run_backward_batch_daily.py" ]]; then
  echo "ERROR: run_backward_batch_daily.py not found under REPO_ROOT=${REPO_ROOT}"
  exit 2
fi

ulimit -s unlimited

# Environment setup (mirrors slurm_array_backward.sh)
if [[ -f "${HOME}/.bashrc" ]]; then
  set +u
  # shellcheck disable=SC1090
  source "${HOME}/.bashrc"
  set -u
fi

if command -v conda >/dev/null 2>&1; then
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
  echo "WARNING: conda not found; continuing without conda activation."
fi

# Compute DATE for this array task
start_epoch="$(date -u -d "${START_DATE}" +%s)"
end_epoch="$(date -u -d "${END_DATE}" +%s)"
total_tasks="$(( (end_epoch - start_epoch) / 86400 + 1 ))"
idx="${SLURM_ARRAY_TASK_ID:-1}"

if (( idx < 1 || idx > total_tasks )); then
  echo "INFO: SLURM_ARRAY_TASK_ID=${idx} outside range 1-${total_tasks}; skipping."
  exit 0
fi

target_epoch="$(( start_epoch + (idx - 1) * 86400 ))"
DATE="$(date -u -d "@${target_epoch}" +%Y%m%d)"

RUN_DIR="${OUTROOT}/${DOMAIN}_${RECEPTOR}_${DATE}"
mkdir -p "${OUTROOT}"

if [[ "${SKIP_IF_RUN_DIR_EXISTS}" == "1" && -d "${RUN_DIR}" ]]; then
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] skip date=${DATE}: existing run directory found at ${RUN_DIR}"
  exit 0
fi

# Optional fast-skip: if final daily postprocessed files already exist for this
# date/site/domain, skip this task before launching FLEXPART.
if [[ "${SKIP_IF_POSTPROCESSED_EXISTS}" == "1" ]]; then
  shopt -s nullglob
  existing_daily_netcdf=(
    "${RUN_DIR}/output/${RECEPTOR}"_*_FLEXPART_CFSv2_"${DOMAIN}"_inert_"${DATE}"*.nc
    "${RUN_DIR}/output/${RECEPTOR}"_*_FLEXPART_GFS_"${DOMAIN}"_inert_"${DATE}"*.nc
  )
  shopt -u nullglob
  if (( ${#existing_daily_netcdf[@]} > 0 )); then
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] skip date=${DATE}: found ${#existing_daily_netcdf[@]} existing postprocessed netcdf file(s) in ${RUN_DIR}/output"
    exit 0
  fi
fi

cd "${REPO_ROOT}"
if [[ "${USE_PROJECT_VENV}" == "1" && -f "${REPO_ROOT}/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "${REPO_ROOT}/.venv/bin/activate"
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
    if "${candidate}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3,6) else 1)' >/dev/null 2>&1; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

PYTHON_CMD="$(pick_python_cmd "${PYTHON_CMD}" python python3 "${REPO_ROOT}/.venv/bin/python" || true)"
if [[ -z "${PYTHON_CMD}" ]]; then
  echo "ERROR: no Python >= 3.6 interpreter found."
  exit 2
fi

if [[ "${DEBUG_ENV}" == "1" ]]; then
  echo "=== Debug Info ==="
  echo "CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
  echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-<empty>}"
  if [[ -n "${FLEXPART_EXE}" && -x "${FLEXPART_EXE}" ]]; then
    ldd "${FLEXPART_EXE}" 2>&1 | head -n 20 || true
  fi
  echo "==================="
fi

echo "Python: $("${PYTHON_CMD}" -c 'import sys; print(sys.executable)')"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] task=${idx}/${total_tasks} date=${DATE} receptor=${RECEPTOR} domain=${DOMAIN}"

# Per-run AVAILABLE file scoped to this run's window
RUN_AVAILABLE="${RUN_DIR}/AVAILABLE"
rm -f "${RUN_AVAILABLE}"

cmd=("${PYTHON_CMD}" "${REPO_ROOT}/run_scripts_daily/run_backward_batch_daily.py"
  --domain    "${DOMAIN}"
  --receptor  "${RECEPTOR}"
  --date      "${DATE}"
  --days      "${BACKWARD_DAYS}"
  --num-particles "${NUM_PARTICLES}"
  --lsubgrid  "${LSUBGRID}"
  --linit-cond "${LINIT_COND}"
  --outdir    "${RUN_DIR}"
  --gfs-available "${RUN_AVAILABLE}"
  --postprocess-footprint-outheight-m   "${POSTPROCESS_FOOTPRINT_OUTHEIGHT_M}"
  --postprocess-source-layer-thickness-m "${POSTPROCESS_SOURCE_LAYER_THICKNESS_M}")

if [[ -n "${IPOUT}" ]]; then
  cmd+=(--ipout "${IPOUT}")
fi

if [[ -n "${NXSHIFT}" ]]; then
  cmd+=(--nxshift "${NXSHIFT}")
fi

if [[ "${DISABLE_AUTO_POSTPROCESS}" == "1" ]]; then
  cmd+=(--no-postprocess)
fi

if [[ -n "${FLEXPART_EXE}" ]]; then
  cmd+=(--executable "${FLEXPART_EXE}")
fi

if [[ -n "${POSTPROCESS_PYTHON_CMD}" ]]; then
  cmd+=(--postprocess-python "${POSTPROCESS_PYTHON_CMD}")
fi

"${cmd[@]}"

# ---------------------------------------------------------------------------
# Post-run cleanup
# ---------------------------------------------------------------------------
OUTPUT_DIR="${RUN_DIR}/output"
if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${PRUNE_OUTPUTS}" == "1" ]]; then
    shopt -s nullglob
    # Keep: grid_time_*.nc, daily postprocessed footprints, and exit-point CSVs
    # Delete: PARTFXTR, DRYDEP, DRYWET, partoutput_*.nc, and other raw outputs
    keep_files=(
      "${OUTPUT_DIR}"/grid_time_*.nc
      "${OUTPUT_DIR}"/*_FLEXPART_CFSv2_"${DOMAIN}"_inert_"${DATE}"*.nc
      "${OUTPUT_DIR}"/*_FLEXPART_GFS_"${DOMAIN}"_inert_"${DATE}"*.nc
      "${OUTPUT_DIR}"/*_domain_exit_points.csv
    )
    if (( ${#keep_files[@]} < 2 )); then
      echo "WARNING: PRUNE_OUTPUTS=1 but no postprocessed daily footprint files found in ${OUTPUT_DIR}."
    else
      # Remove everything in output/ except postprocessed files
      for p in "${OUTPUT_DIR}"/*; do
        keep=false
        for k in "${keep_files[@]}"; do
          [[ "${p}" == "${k}" ]] && keep=true && break
        done
        [[ "${keep}" == false ]] && rm -rf "${p}"
      done
      # Remove all config/input files from RUN_DIR (keep only output/ subdir)
      for p in "${RUN_DIR}"/*; do
        [[ "${p}" != "${OUTPUT_DIR}" ]] && rm -rf "${p}"
      done
      echo "Pruned: kept ${#keep_files[@]} grid/postprocessed file(s); deleted all other run-dir contents."
    fi
    shopt -u nullglob
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] completed date=${DATE}"
    exit 0
  fi

  # Default: keep final postprocessed footprint NetCDFs for this date (all 24 hours)
  shopt -s nullglob
  keep_files=(
    "${OUTPUT_DIR}"/*_FLEXPART_CFSv2_"${DOMAIN}"_inert_"${DATE}"*.nc
    "${OUTPUT_DIR}"/*_FLEXPART_GFS_"${DOMAIN}"_inert_"${DATE}"*.nc
  )
  if (( ${#keep_files[@]} == 0 )); then
    echo "WARNING: no final footprint files found in ${OUTPUT_DIR} for ${DATE}."
  else
    for p in "${OUTPUT_DIR}"/*; do
      keep=false
      for k in "${keep_files[@]}"; do
        [[ "${p}" == "${k}" ]] && keep=true && break
      done
      [[ "${keep}" == false ]] && rm -rf "${p}"
    done
  fi
  shopt -u nullglob
fi

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] completed date=${DATE}"
