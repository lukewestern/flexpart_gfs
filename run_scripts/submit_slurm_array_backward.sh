#!/bin/bash
set -euo pipefail

# Helper to submit run_scripts/slurm_array_backward.sh with a date-derived array range.
#
# Usage:
#   run_scripts/submit_slurm_array_backward.sh START_DATE END_DATE [MAX_CONCURRENT]
#
# Notes:
#   - The script computes task count from START_DATE..END_DATE and STEP_HOURS.
#   - It submits sbatch with a dynamic --array=1-N%MAX_CONCURRENT.
#
# Example:
#   run_scripts/submit_slurm_array_backward.sh 2018020100 2018022800 20
#
# Full example with all optional inputs:
#   DOMAIN=EASTASIA \
#   RECEPTOR=GSN \
#   STEP_HOURS=1 \
#   BACKWARD_DAYS=20 \
#   NUM_PARTICLES=20000 \
#   POSTPROCESS_LOWEST_MAGL=100 \
#   POSTPROCESS_SOURCE_LAYER_THICKNESS_M=100 \
#   OUTROOT=/net/fs06/d2/$USER/flexpart_outs \
#   run_scripts/submit_slurm_array_backward.sh 2018020100 2018022800 20
#
# Optional env overrides:
#   DOMAIN, RECEPTOR, STEP_HOURS, BACKWARD_DAYS, NUM_PARTICLES,
#   POSTPROCESS_LOWEST_MAGL, POSTPROCESS_SOURCE_LAYER_THICKNESS_M, OUTROOT,
#   PYTHON_CMD, USE_PROJECT_VENV, FLEXPART_EXE

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 START_DATE END_DATE [MAX_CONCURRENT]"
  exit 2
fi

START_DATE="$1"
END_DATE="$2"
MAX_CONCURRENT="${3:-20}"

if [[ ! "${START_DATE}" =~ ^[0-9]{10}$ || ! "${END_DATE}" =~ ^[0-9]{10}$ ]]; then
  echo "ERROR: START_DATE/END_DATE must match YYYYMMDDHH."
  exit 2
fi

if ! [[ "${MAX_CONCURRENT}" =~ ^[0-9]+$ ]] || [[ "${MAX_CONCURRENT}" -lt 1 ]]; then
  echo "ERROR: MAX_CONCURRENT must be an integer >= 1."
  exit 2
fi

STEP_HOURS="${STEP_HOURS:-1}"
if ! [[ "${STEP_HOURS}" =~ ^[0-9]+$ ]] || [[ "${STEP_HOURS}" -lt 1 ]]; then
  echo "ERROR: STEP_HOURS must be an integer >= 1."
  exit 2
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

n_tasks="$(( span_seconds / step_seconds + 1 ))"
array_spec="1-${n_tasks}%${MAX_CONCURRENT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/slurm_array_backward.sh"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ ! -f "${JOB_SCRIPT}" ]]; then
  echo "ERROR: job script not found: ${JOB_SCRIPT}"
  exit 2
fi

echo "Submitting array job: ${array_spec}"
echo "  START_DATE=${START_DATE}"
echo "  END_DATE=${END_DATE}"
echo "  STEP_HOURS=${STEP_HOURS}"
echo "  FLEXPART_REPO_ROOT=${REPO_ROOT}"
if [[ -n "${FLEXPART_EXE:-}" ]]; then
  echo "  FLEXPART_EXE=${FLEXPART_EXE}"
fi

sbatch \
  --array="${array_spec}" \
  --export=ALL,START_DATE="${START_DATE}",END_DATE="${END_DATE}",FLEXPART_REPO_ROOT="${REPO_ROOT}",FLEXPART_EXE="${FLEXPART_EXE:-}" \
  "${JOB_SCRIPT}"
