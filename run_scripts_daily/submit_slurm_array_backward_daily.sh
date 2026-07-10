#!/bin/bash
set -euo pipefail

# Helper: compute day-array spec and submit slurm_array_backward_daily.sh.
#
# Usage:
#   submit_slurm_array_backward_daily.sh START_DATE END_DATE [MAX_CONCURRENT]
#
# START_DATE and END_DATE: YYYYMMDD (inclusive on both ends)
# MAX_CONCURRENT: Slurm array throttle %N (default: 20)

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 START_DATE END_DATE [MAX_CONCURRENT]"
  exit 2
fi

START_DATE="$1"
END_DATE="$2"
MAX_CONCURRENT="${3:-20}"

if [[ ! "${START_DATE}" =~ ^[0-9]{8}$ || ! "${END_DATE}" =~ ^[0-9]{8}$ ]]; then
  echo "ERROR: START_DATE/END_DATE must match YYYYMMDD."
  exit 2
fi

if ! [[ "${MAX_CONCURRENT}" =~ ^[0-9]+$ ]] || [[ "${MAX_CONCURRENT}" -lt 1 ]]; then
  echo "ERROR: MAX_CONCURRENT must be an integer >= 1."
  exit 2
fi

start_epoch="$(date -u -d "${START_DATE}" +%s)"
end_epoch="$(date -u -d "${END_DATE}" +%s)"

if [[ "${end_epoch}" -lt "${start_epoch}" ]]; then
  echo "ERROR: END_DATE is earlier than START_DATE."
  exit 2
fi

n_days="$(( (end_epoch - start_epoch) / 86400 + 1 ))"
array_spec="1-${n_days}%${MAX_CONCURRENT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/slurm_array_backward_daily.sh"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ ! -f "${JOB_SCRIPT}" ]]; then
  echo "ERROR: job script not found: ${JOB_SCRIPT}"
  exit 2
fi

echo "Submitting daily array job: ${array_spec}"
echo "  START_DATE=${START_DATE}"
echo "  END_DATE=${END_DATE}"
echo "  n_days=${n_days}"
echo "  FLEXPART_REPO_ROOT=${REPO_ROOT}"
if [[ -n "${FLEXPART_EXE:-}" ]]; then
  echo "  FLEXPART_EXE=${FLEXPART_EXE}"
fi

sbatch \
  --array="${array_spec}" \
  --export=ALL,START_DATE="${START_DATE}",END_DATE="${END_DATE}",FLEXPART_REPO_ROOT="${REPO_ROOT}",FLEXPART_EXE="${FLEXPART_EXE:-}" \
  "${JOB_SCRIPT}"
