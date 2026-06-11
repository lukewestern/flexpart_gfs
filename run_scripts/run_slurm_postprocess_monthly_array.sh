#!/bin/bash
set -euo pipefail

# User-facing launcher for month-sharded postprocessing over outputs in OUTROOT.
#
# Usage:
#   ./run_scripts/run_slurm_postprocess_monthly_array.sh [CONFIG_FILE]
#
# Optional month bounds from config/env:
#   START_DATE=YYYYMMDDHH
#   END_DATE=YYYYMMDDHH
# These are used only to restrict discovered months before submission.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

resolve_job_script() {
  local candidate
  for candidate in \
    "${SCRIPT_DIR}/slurm_postprocess_monthly_array.sh" \
    "${REPO_ROOT}/run_scripts/slurm_postprocess_monthly_array.sh" \
    "${SLURM_SUBMIT_DIR:-}/slurm_postprocess_monthly_array.sh" \
    "${SLURM_SUBMIT_DIR:-}/run_scripts/slurm_postprocess_monthly_array.sh"; do
    if [[ -n "${candidate}" && -f "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

JOB_SCRIPT="$(resolve_job_script || true)"
if [[ -z "${JOB_SCRIPT}" ]]; then
  echo "ERROR: job script not found: run_scripts/slurm_postprocess_monthly_array.sh"
  echo "       Tried SCRIPT_DIR=${SCRIPT_DIR} and SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}"
  echo "       Run from repo root as: ./run_scripts/run_slurm_postprocess_monthly_array.sh"
  exit 2
fi

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
POSTPROCESS_MONTHLY_MAX_CONCURRENT=""
START_DATE="${START_DATE:-}"
END_DATE="${END_DATE:-}"

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

if [[ -n "${POSTPROCESS_LOWEST_MAGL:-}" ]]; then
  POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="${POSTPROCESS_LOWEST_MAGL}"
fi

if [[ -z "${POSTPROCESS_MONTHLY_MAX_CONCURRENT}" ]]; then
  if [[ -n "${MAX_CONCURRENT:-}" ]]; then
    POSTPROCESS_MONTHLY_MAX_CONCURRENT="${MAX_CONCURRENT}"
  else
    POSTPROCESS_MONTHLY_MAX_CONCURRENT="6"
  fi
fi

if ! [[ "${POSTPROCESS_MONTHLY_MAX_CONCURRENT}" =~ ^[0-9]+$ ]] || [[ "${POSTPROCESS_MONTHLY_MAX_CONCURRENT}" -lt 1 ]]; then
  echo "ERROR: POSTPROCESS_MONTHLY_MAX_CONCURRENT must be an integer >= 1"
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
  offending_pypi_lines="$(conda list -n flexpart-post 2>/dev/null | grep -E '^(numpy|pandas|xarray|netcdf4)\s+.+\s+pypi_0\s+pypi$' || true)"
  if [[ -n "${offending_pypi_lines}" ]]; then
    echo "ERROR: flexpart-post scientific stack was installed from pip (pypi), which is incompatible on older compute-node glibc."
    echo "       Offending package rows from 'conda list -n flexpart-post':"
    while IFS= read -r row; do
      echo "         ${row}"
    done <<< "${offending_pypi_lines}"
    echo "       Repair with conda packages before submitting:"
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

month_list="$(${POSTPROCESS_DRIVER_PYTHON} - "${OUTROOT}" "${START_DATE}" "${END_DATE}" <<'PY'
from __future__ import annotations

from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
start = sys.argv[2].strip()
end = sys.argv[3].strip()

if not root.is_dir():
    print("")
    raise SystemExit(0)

re_grid = re.compile(r"^grid_time_(\d{14})\.nc$")
re_hour = re.compile(r"^\d{10}$")

start_month = ""
end_month = ""
if start:
    if not re_hour.match(start):
        raise SystemExit(f"START_DATE must match YYYYMMDDHH, got: {start}")
    start_month = start[:6]
if end:
    if not re_hour.match(end):
        raise SystemExit(f"END_DATE must match YYYYMMDDHH, got: {end}")
    end_month = end[:6]

months: set[str] = set()
for path in root.glob("**/output/grid_time_*.nc"):
    m = re_grid.match(path.name)
    if not m:
        continue
    month = m.group(1)[:6]
    if start_month and month < start_month:
        continue
    if end_month and month > end_month:
        continue
    months.add(month)

print("\n".join(sorted(months)))
PY
)"

if [[ -z "${month_list}" ]]; then
  echo "No months discovered under ${OUTROOT}; nothing to submit."
  exit 0
fi

mapfile -t month_arr <<< "${month_list}"
month_count="${#month_arr[@]}"
array_spec="1-${month_count}%${POSTPROCESS_MONTHLY_MAX_CONCURRENT}"
month_list_export="${month_list//$'\n'/ }"

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
export START_DATE
export END_DATE

echo "Submitting monthly postprocess array job"
echo "  OUTROOT=${OUTROOT}"
echo "  FINAL_DIR=${FINAL_DIR}"
echo "  MONTHLY_DIR=${MONTHLY_DIR}"
echo "  MONTHS=${month_list_export}"
echo "  ARRAY=${array_spec}"
echo "  POSTPROCESS_DRIVER_PYTHON=${POSTPROCESS_DRIVER_PYTHON}"
echo "  POSTPROCESS_PYTHON_CMD=${POSTPROCESS_PYTHON_CMD}"
echo "  POSTPROCESS_WORKERS=${POSTPROCESS_WORKERS}"
echo "  POSTPROCESS_OVERWRITE=${POSTPROCESS_OVERWRITE}"
echo "  KEEP_RUN_DIRS=${KEEP_RUN_DIRS}"
echo "  WRITE_MONTHLY=${WRITE_MONTHLY}"
echo "  DRY_RUN=${DRY_RUN}"

sbatch \
  --array="${array_spec}" \
  --export=ALL,FLEXPART_REPO_ROOT="${REPO_ROOT}",MONTH_LIST="${month_list_export}" \
  "${JOB_SCRIPT}"
