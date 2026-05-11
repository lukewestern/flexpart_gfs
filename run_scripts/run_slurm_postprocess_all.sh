#!/bin/bash
set -euo pipefail

# User-facing launcher for single-job postprocessing over all outputs in OUTROOT.
#
# Usage:
#   ./run_scripts/run_slurm_postprocess_all.sh [CONFIG_FILE]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# If this script is submitted via sbatch, BASH_SOURCE may point to a temp copy
# under /tmp/slurmd/... where sibling scripts are unavailable.
resolve_job_script() {
	local candidate
	for candidate in \
		"${SCRIPT_DIR}/slurm_postprocess_all.sh" \
		"${REPO_ROOT}/run_scripts/slurm_postprocess_all.sh" \
		"${SLURM_SUBMIT_DIR:-}/run_scripts/slurm_postprocess_all.sh"; do
		if [[ -n "${candidate}" && -f "${candidate}" ]]; then
			echo "${candidate}"
			return 0
		fi
	done
	return 1
}

JOB_SCRIPT="$(resolve_job_script || true)"
if [[ -z "${JOB_SCRIPT}" ]]; then
	echo "ERROR: job script not found: run_scripts/slurm_postprocess_all.sh"
	echo "       Tried SCRIPT_DIR=${SCRIPT_DIR} and SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}"
	echo "       Run from repo root as: ./run_scripts/run_slurm_postprocess_all.sh"
	exit 2
fi

# Normalize REPO_ROOT from the resolved job script location.
REPO_ROOT="$(cd "$(dirname "${JOB_SCRIPT}")/.." && pwd)"
CONFIG_FILE="${1:-${REPO_ROOT}/run_scripts/slurm_array_config.sh}"

OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"
FINAL_DIR=""
POSTPROCESS_DRIVER_PYTHON="python3"
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
POSTPROCESS_WORKERS="8"
POSTPROCESS_OVERWRITE="0"
KEEP_RUN_DIRS="1"
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

# Backward-compatible alias for older configs.
if [[ -n "${POSTPROCESS_LOWEST_MAGL:-}" ]]; then
	POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_LOWEST_MAGL}"
fi

export OUTROOT
export FINAL_DIR
export POSTPROCESS_DRIVER_PYTHON
export POSTPROCESS_PYTHON_CMD
export POSTPROCESS_FOOTPRINT_OUTHEIGHT_M
export POSTPROCESS_SOURCE_LAYER_THICKNESS_M
export POSTPROCESS_WORKERS
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
echo "  POSTPROCESS_FOOTPRINT_OUTHEIGHT_M=${POSTPROCESS_FOOTPRINT_OUTHEIGHT_M}"
echo "  POSTPROCESS_WORKERS=${POSTPROCESS_WORKERS}"
echo "  POSTPROCESS_OVERWRITE=${POSTPROCESS_OVERWRITE}"
echo "  KEEP_RUN_DIRS=${KEEP_RUN_DIRS}"
echo "  LIMIT=${LIMIT}"
echo "  WRITE_MONTHLY=${WRITE_MONTHLY}"
echo "  MONTHLY_DIR=${MONTHLY_DIR}"
echo "  DRY_RUN=${DRY_RUN}"

sbatch \
	--export=ALL,FLEXPART_REPO_ROOT="${REPO_ROOT}" \
	"${JOB_SCRIPT}"
