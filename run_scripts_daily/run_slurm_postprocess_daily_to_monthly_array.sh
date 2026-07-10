#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_scripts_daily/run_slurm_postprocess_daily_to_monthly_array.sh START_MONTH END_MONTH SITE DOMAIN [options]

Example:
  ./run_scripts_daily/run_slurm_postprocess_daily_to_monthly_array.sh 201801 201812 THD WESTUSA

Options:
  --root-dir DIR              Daily FLEXPART output root
                              default: /net/fs06/d2/${USER}/flexpart_outs_daily
  --python PYTHON             Python used to launch postprocess_daily_to_monthly.py
                              default: python3
  --postprocess-python PATH   Python used by postprocess subprocesses
                              default: /home/${USER}/.conda/envs/flexpart-post/bin/python
  --workers N                 Workers per Slurm task; must be an integer >= 1
                              default: 2
  --max-parallel N            Maximum array tasks running at once
                              default: 12
  --overwrite                 Overwrite existing daily/monthly products
  --dry-run                   Print work without writing or deleting
  --delete-run-dirs           After verification, delete DOMAIN_SITE_YYYYMMDD run dirs
  --delete-daily-files        After verification, delete daily postprocessed NetCDFs
  -h, --help                  Show this help

The job always includes exit diagnostics in the NetCDF output and never writes
exit CSV files.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/slurm_postprocess_daily_to_monthly_array.sh"

[[ -f "${JOB_SCRIPT}" ]] || die "job script not found: ${JOB_SCRIPT}"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ $# -ge 4 ]] || { usage; exit 2; }

START_MONTH="$1"
END_MONTH="$2"
SITE="$3"
DOMAIN="$4"
shift 4

OUTROOT="/net/fs06/d2/${USER}/flexpart_outs_daily"
POSTPROCESS_DRIVER_PYTHON="python3"
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
POSTPROCESS_WORKERS="1"
POSTPROCESS_MONTHLY_MAX_CONCURRENT="12"
POSTPROCESS_OVERWRITE="0"
DELETE_FLEXPART_RUN_DIRS="0"
DELETE_DAILY_POSTPROCESS_FILES="0"
DRY_RUN="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root-dir)
      [[ $# -ge 2 ]] || die "--root-dir requires a value"
      OUTROOT="$2"
      shift 2
      ;;
    --python)
      [[ $# -ge 2 ]] || die "--python requires a value"
      POSTPROCESS_DRIVER_PYTHON="$2"
      shift 2
      ;;
    --postprocess-python)
      [[ $# -ge 2 ]] || die "--postprocess-python requires a value"
      POSTPROCESS_PYTHON_CMD="$2"
      shift 2
      ;;
    --workers)
      [[ $# -ge 2 ]] || die "--workers requires a value"
      POSTPROCESS_WORKERS="$2"
      shift 2
      ;;
    --max-parallel)
      [[ $# -ge 2 ]] || die "--max-parallel requires a value"
      POSTPROCESS_MONTHLY_MAX_CONCURRENT="$2"
      shift 2
      ;;
    --overwrite)
      POSTPROCESS_OVERWRITE="1"
      shift
      ;;
    --dry-run)
      DRY_RUN="1"
      shift
      ;;
    --delete-run-dirs)
      DELETE_FLEXPART_RUN_DIRS="1"
      shift
      ;;
    --delete-daily-files)
      DELETE_DAILY_POSTPROCESS_FILES="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[[ "${START_MONTH}" =~ ^[0-9]{6}$ ]] || die "START_MONTH must match YYYYMM, got ${START_MONTH}"
[[ "${END_MONTH}" =~ ^[0-9]{6}$ ]] || die "END_MONTH must match YYYYMM, got ${END_MONTH}"
[[ "${START_MONTH}" -le "${END_MONTH}" ]] || die "START_MONTH must be <= END_MONTH"
[[ -n "${SITE}" ]] || die "SITE is required"
[[ -n "${DOMAIN}" ]] || die "DOMAIN is required"

if ! [[ "${POSTPROCESS_WORKERS}" =~ ^[0-9]+$ ]] || [[ "${POSTPROCESS_WORKERS}" -lt 1 ]]; then
  die "--workers must be an integer >= 1; got ${POSTPROCESS_WORKERS}"
fi

if ! [[ "${POSTPROCESS_MONTHLY_MAX_CONCURRENT}" =~ ^[0-9]+$ ]] || [[ "${POSTPROCESS_MONTHLY_MAX_CONCURRENT}" -lt 1 ]]; then
  die "--max-parallel must be an integer >= 1"
fi

command -v "${POSTPROCESS_DRIVER_PYTHON}" >/dev/null 2>&1 || die "python not found: ${POSTPROCESS_DRIVER_PYTHON}"
[[ -x "${POSTPROCESS_PYTHON_CMD}" ]] || die "postprocess python is not executable: ${POSTPROCESS_PYTHON_CMD}"

month_list="$(${POSTPROCESS_DRIVER_PYTHON} - "${START_MONTH}" "${END_MONTH}" <<'PY'
from __future__ import annotations

import sys

start = sys.argv[1]
end = sys.argv[2]

year = int(start[:4])
month = int(start[4:6])
end_year = int(end[:4])
end_month = int(end[4:6])

months = []
while (year, month) <= (end_year, end_month):
    months.append(f"{year:04d}{month:02d}")
    month += 1
    if month == 13:
        year += 1
        month = 1

print(" ".join(months))
PY
)"

read -r -a month_arr <<< "${month_list}"
month_count="${#month_arr[@]}"
array_limit="${POSTPROCESS_MONTHLY_MAX_CONCURRENT}"
if (( array_limit > month_count )); then
  array_limit="${month_count}"
fi
array_spec="1-${month_count}%${array_limit}"

export OUTROOT
export SITE
export DOMAIN
export POSTPROCESS_DRIVER_PYTHON
export POSTPROCESS_PYTHON_CMD
export POSTPROCESS_WORKERS
export POSTPROCESS_OVERWRITE
export DELETE_FLEXPART_RUN_DIRS
export DELETE_DAILY_POSTPROCESS_FILES
export DRY_RUN

echo "Submitting daily-to-monthly postprocess array job"
echo "  OUTROOT=${OUTROOT}"
echo "  SITE=${SITE}"
echo "  DOMAIN=${DOMAIN}"
echo "  MONTHS=${month_list}"
echo "  ARRAY=${array_spec}"
echo "  POSTPROCESS_DRIVER_PYTHON=${POSTPROCESS_DRIVER_PYTHON}"
echo "  POSTPROCESS_PYTHON_CMD=${POSTPROCESS_PYTHON_CMD}"
echo "  POSTPROCESS_WORKERS=${POSTPROCESS_WORKERS}"
echo "  DELETE_FLEXPART_RUN_DIRS=${DELETE_FLEXPART_RUN_DIRS}"
echo "  DELETE_DAILY_POSTPROCESS_FILES=${DELETE_DAILY_POSTPROCESS_FILES}"
echo "  DRY_RUN=${DRY_RUN}"

sbatch \
  --array="${array_spec}" \
  --cpus-per-task="${POSTPROCESS_WORKERS}" \
  --export=ALL,FLEXPART_REPO_ROOT="${REPO_ROOT}",MONTH_LIST="${month_list}" \
  "${JOB_SCRIPT}"
