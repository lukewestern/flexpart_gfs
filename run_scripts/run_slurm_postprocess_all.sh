#!/bin/bash
set -euo pipefail

# User-facing launcher for single-job postprocessing over all outputs in OUTROOT.
#
# Usage:
#   ./run_scripts/run_slurm_postprocess_all.sh [CONFIG_FILE]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SUBMIT_SCRIPT="${SCRIPT_DIR}/submit_slurm_postprocess_all.sh"
CONFIG_FILE="${1:-${REPO_ROOT}/run_scripts/slurm_array_config.sh}"

OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"
FINAL_DIR=""
POSTPROCESS_DRIVER_PYTHON="python3"
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
POSTPROCESS_LOWEST_MAGL="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
POSTPROCESS_OVERWRITE="0"
KEEP_RUN_DIRS="0"
LIMIT="0"
WRITE_MONTHLY="1"
MONTHLY_DIR=""
DRY_RUN="0"

if [[ -f "${CONFIG_FILE}" ]]; then
	# shellcheck disable=SC1090
	source "${CONFIG_FILE}"
	echo "Loaded config overrides from ${CONFIG_FILE}"
fi

if [[ -z "${FINAL_DIR}" ]]; then
	if [[ -n "${RECEPTOR:-}" ]]; then
		FINAL_DIR="${OUTROOT}/$(printf '%s' "${RECEPTOR}" | tr '[:upper:]' '[:lower:]')_hourly"
	else
		FINAL_DIR="${OUTROOT}"
	fi
fi

if [[ -z "${MONTHLY_DIR}" ]]; then
	if [[ -n "${RECEPTOR:-}" ]]; then
		MONTHLY_DIR="${OUTROOT}/$(printf '%s' "${RECEPTOR}" | tr '[:upper:]' '[:lower:]')"
	else
		MONTHLY_DIR="${OUTROOT}"
	fi
fi

if [[ ! -x "${SUBMIT_SCRIPT}" ]]; then
  echo "ERROR: submit script not executable: ${SUBMIT_SCRIPT}"
  exit 2
fi

export OUTROOT
export FINAL_DIR
export POSTPROCESS_DRIVER_PYTHON
export POSTPROCESS_PYTHON_CMD
export POSTPROCESS_LOWEST_MAGL
export POSTPROCESS_SOURCE_LAYER_THICKNESS_M
export POSTPROCESS_OVERWRITE
export KEEP_RUN_DIRS
export LIMIT
export WRITE_MONTHLY
export MONTHLY_DIR
export DRY_RUN

echo "Submitting postprocess-all Slurm job"
echo "  OUTROOT=${OUTROOT}"
echo "  FINAL_DIR=${FINAL_DIR}"
echo "  POSTPROCESS_DRIVER_PYTHON=${POSTPROCESS_DRIVER_PYTHON}"
echo "  POSTPROCESS_PYTHON_CMD=${POSTPROCESS_PYTHON_CMD}"
echo "  POSTPROCESS_OVERWRITE=${POSTPROCESS_OVERWRITE}"
echo "  KEEP_RUN_DIRS=${KEEP_RUN_DIRS}"
echo "  LIMIT=${LIMIT}"
echo "  WRITE_MONTHLY=${WRITE_MONTHLY}"
echo "  MONTHLY_DIR=${MONTHLY_DIR}"
echo "  DRY_RUN=${DRY_RUN}"

"${SUBMIT_SCRIPT}"
