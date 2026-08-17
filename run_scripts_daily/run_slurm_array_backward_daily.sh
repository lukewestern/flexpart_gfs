#!/bin/bash
set -euo pipefail

# User-facing launcher for daily FLEXPART backward Slurm array runs.
#
# How to run:
#   1) Edit run_scripts_daily/slurm_array_config_daily.sh (or pass a custom config as $1).
#   2) From repo root:
#        ./run_scripts_daily/run_slurm_array_backward_daily.sh
#
# Each array task runs all 24 hourly releases for one calendar day in a single
# FLEXPART simulation — roughly 20x less met-field I/O than 24 separate hourly runs.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

resolve_submit_script() {
  local candidate
  for candidate in \
    "${SCRIPT_DIR}/submit_slurm_array_backward_daily.sh" \
    "${REPO_ROOT}/run_scripts_daily/submit_slurm_array_backward_daily.sh"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

SUBMIT_SCRIPT="$(resolve_submit_script || true)"
if [[ -z "${SUBMIT_SCRIPT}" ]]; then
  echo "ERROR: could not locate submit_slurm_array_backward_daily.sh"
  exit 2
fi

# Optional config file override
CONFIG_FILE="${1:-${REPO_ROOT}/run_scripts_daily/slurm_array_config_daily.sh}"

# ---------------------------
# Default configuration
# ---------------------------
START_DATE=""
END_DATE=""
MAX_CONCURRENT="20"

DOMAIN="EASTASIA"
RECEPTOR="GSN"
BACKWARD_DAYS="20"
NUM_PARTICLES="20000"
IPOUT="0"
NXSHIFT=""
LSUBGRID="1"
LINIT_COND="1"
POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
DISABLE_AUTO_POSTPROCESS="1"
PRUNE_OUTPUTS="1"
SKIP_IF_POSTPROCESSED_EXISTS="1"
SKIP_IF_RUN_DIR_EXISTS="0"
OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"
FLEXPART_EXE="${REPO_ROOT}/src/FLEXPART"
PYTHON_CMD=""
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
USE_PROJECT_VENV="0"
DEBUG_ENV="0"
POSTPROCESS_LOWEST_MAGL=""
LINIT_COND="1"
KEEP_RUN_DIRS="1"

# Load user overrides
if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG_FILE}"
  echo "Loaded config overrides from ${CONFIG_FILE}"
fi

if [[ -z "${START_DATE}" || -z "${END_DATE}" ]]; then
  echo "ERROR: START_DATE and END_DATE must be set in the config file or environment (YYYYMMDD)."
  exit 2
fi

# Export for slurm worker
export DOMAIN RECEPTOR BACKWARD_DAYS NUM_PARTICLES IPOUT NXSHIFT LSUBGRID LINIT_COND
export POSTPROCESS_FOOTPRINT_OUTHEIGHT_M POSTPROCESS_LOWEST_MAGL POSTPROCESS_SOURCE_LAYER_THICKNESS_M
export DISABLE_AUTO_POSTPROCESS PRUNE_OUTPUTS SKIP_IF_POSTPROCESSED_EXISTS SKIP_IF_RUN_DIR_EXISTS
export OUTROOT FLEXPART_EXE PYTHON_CMD POSTPROCESS_PYTHON_CMD USE_PROJECT_VENV DEBUG_ENV

echo "Submitting daily FLEXPART array run"
echo "  START_DATE=${START_DATE}"
echo "  END_DATE=${END_DATE}"
echo "  MAX_CONCURRENT=${MAX_CONCURRENT}"
echo "  DOMAIN=${DOMAIN}"
echo "  RECEPTOR=${RECEPTOR}"
echo "  OUTROOT=${OUTROOT}"
echo "  FLEXPART_EXE=${FLEXPART_EXE}"
echo "  IPOUT=${IPOUT}"
echo "  NXSHIFT=${NXSHIFT:-<auto>}"
echo "  BACKWARD_DAYS=${BACKWARD_DAYS}"
echo "  NUM_PARTICLES=${NUM_PARTICLES}  (x24 = $((NUM_PARTICLES * 24)) total per run)"
echo "  LINIT_COND=${LINIT_COND}"
echo "  DISABLE_AUTO_POSTPROCESS=${DISABLE_AUTO_POSTPROCESS}"
echo "  SKIP_IF_POSTPROCESSED_EXISTS=${SKIP_IF_POSTPROCESSED_EXISTS}"
echo "  SKIP_IF_RUN_DIR_EXISTS=${SKIP_IF_RUN_DIR_EXISTS}"
echo "  PRUNE_OUTPUTS=${PRUNE_OUTPUTS}"
"${SUBMIT_SCRIPT}" "${START_DATE}" "${END_DATE}" "${MAX_CONCURRENT}"
