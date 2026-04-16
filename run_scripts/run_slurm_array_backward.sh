#!/bin/bash
set -euo pipefail

# User-facing launcher for FLEXPART backward Slurm array runs.
#
# How to run:
#   1) Edit the configuration block below.
#   2) From repo root, run:
#        ./run_scripts/run_slurm_array_backward.sh
#
# What it does:
#   - Exports the configured options as environment variables.
#   - Calls submit_slurm_array_backward.sh with START_DATE, END_DATE,
#     and MAX_CONCURRENT.
#   - Submits the Slurm array job via sbatch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# If this script is submitted via sbatch, BASH_SOURCE points to a temp copy under
# /tmp/slurmd/... and sibling scripts are not available there. Resolve helper paths
# from the original submit directory/repo when possible.
resolve_submit_script() {
	local candidate
	for candidate in \
		"${SCRIPT_DIR}/submit_slurm_array_backward.sh" \
		"${REPO_ROOT}/run_scripts/submit_slurm_array_backward.sh" \
		"${SLURM_SUBMIT_DIR:-}/run_scripts/submit_slurm_array_backward.sh"; do
		if [[ -n "${candidate}" && -x "${candidate}" ]]; then
			echo "${candidate}"
			return 0
		fi
	done
	return 1
}

SUBMIT_SCRIPT="$(resolve_submit_script || true)"
if [[ -z "${SUBMIT_SCRIPT}" ]]; then
	echo "ERROR: could not locate run_scripts/submit_slurm_array_backward.sh"
	echo "       Tried SCRIPT_DIR=${SCRIPT_DIR} and SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}"
	echo "       Run from repo root as: ./run_scripts/run_slurm_array_backward.sh"
	exit 2
fi

# Optional config file to override defaults below.
# Usage:
#   ./run_scripts/run_slurm_array_backward.sh [CONFIG_FILE]
CONFIG_FILE="${1:-${REPO_ROOT}/run_scripts/slurm_array_config.sh}"

# ---------------------------
# Configuration
# ---------------------------
START_DATE="2018020100"   # inclusive, YYYYMMDDHH
END_DATE="2018020123"     # inclusive, YYYYMMDDHH
MAX_CONCURRENT="30"       # Slurm array throttle (%N)

DOMAIN="EASTASIA"
RECEPTOR="GSN"
STEP_HOURS="1"
BACKWARD_DAYS="20"
NUM_PARTICLES="20000"
POSTPROCESS_LOWEST_MAGL="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
DISABLE_AUTO_POSTPROCESS="0"
PRUNE_TO_GRID_FILES="0"
OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"
FLEXPART_EXE="${REPO_ROOT}/src/FLEXPART"
PYTHON_CMD=""
POSTPROCESS_PYTHON_CMD="/home/lwestern/.conda/envs/flexpart/bin/python"
USE_PROJECT_VENV="0"
DEBUG_ENV="0"

# Load user overrides if provided.
if [[ -f "${CONFIG_FILE}" ]]; then
	# shellcheck disable=SC1090
	source "${CONFIG_FILE}"
	echo "Loaded config overrides from ${CONFIG_FILE}"
fi

# ---------------------------
# Submit
# ---------------------------
export DOMAIN
export RECEPTOR
export STEP_HOURS
export BACKWARD_DAYS
export NUM_PARTICLES
export POSTPROCESS_LOWEST_MAGL
export POSTPROCESS_SOURCE_LAYER_THICKNESS_M
export DISABLE_AUTO_POSTPROCESS
export PRUNE_TO_GRID_FILES
export OUTROOT
export FLEXPART_EXE
export PYTHON_CMD
export POSTPROCESS_PYTHON_CMD
export USE_PROJECT_VENV
export DEBUG_ENV

echo "Submitting configured FLEXPART array run"
echo "  START_DATE=${START_DATE}"
echo "  END_DATE=${END_DATE}"
echo "  MAX_CONCURRENT=${MAX_CONCURRENT}"
echo "  DOMAIN=${DOMAIN}"
echo "  RECEPTOR=${RECEPTOR}"
echo "  OUTROOT=${OUTROOT}"
echo "  FLEXPART_EXE=${FLEXPART_EXE}"
echo "  PYTHON_CMD=${PYTHON_CMD:-<auto>}"
echo "  POSTPROCESS_PYTHON_CMD=${POSTPROCESS_PYTHON_CMD:-<same as PYTHON_CMD>}"
echo "  USE_PROJECT_VENV=${USE_PROJECT_VENV}"
echo "  DISABLE_AUTO_POSTPROCESS=${DISABLE_AUTO_POSTPROCESS}"
echo "  PRUNE_TO_GRID_FILES=${PRUNE_TO_GRID_FILES}"

"${SUBMIT_SCRIPT}" "${START_DATE}" "${END_DATE}" "${MAX_CONCURRENT}"
