#!/bin/bash
#SBATCH -J flexpart_daily_to_monthly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH -t 6:00:00
#SBATCH -p fdr,edr
#SBATCH --mem=20G
#SBATCH -o slurm-daily-to-monthly-%A_%a.out
#SBATCH -e slurm-daily-to-monthly-%A_%a.err

set -euo pipefail

# Required environment variables:
#   MONTH_LIST: whitespace-separated YYYYMM values, one month per array task
# Required environment variables:
#   SITE: site/receptor code, e.g. THD
#   DOMAIN: domain name, e.g. WESTUSA
# Optional environment variables:
#   OUTROOT: root directory containing DOMAIN_SITE_YYYYMMDD run folders
#   POSTPROCESS_DRIVER_PYTHON: python used to launch postprocess_daily_to_monthly.py
#   POSTPROCESS_PYTHON_CMD: python used by postprocess subprocesses
#   POSTPROCESS_WORKERS: concurrent daily postprocess workers per array task
#   POSTPROCESS_OVERWRITE: set to 1 to overwrite existing products
#   DELETE_FLEXPART_RUN_DIRS: set to 1 to delete monthly FLEXPART run dirs after verification
#   DELETE_DAILY_POSTPROCESS_FILES: set to 1 to delete daily NetCDF intermediates after verification
#   DRY_RUN: set to 1 to preview commands without writing/deleting
#   FLEXPART_REPO_ROOT: absolute path to repository root

if [[ -z "${MONTH_LIST:-}" ]]; then
  echo "ERROR: MONTH_LIST is required for daily-to-monthly array jobs."
  exit 2
fi

OUTROOT="${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs_daily}"
SITE="${SITE:-}"
DOMAIN="${DOMAIN:-}"
POSTPROCESS_DRIVER_PYTHON="${POSTPROCESS_DRIVER_PYTHON:-python3}"
POSTPROCESS_PYTHON_CMD="${POSTPROCESS_PYTHON_CMD:-/home/${USER}/.conda/envs/flexpart-post/bin/python}"
POSTPROCESS_WORKERS="${POSTPROCESS_WORKERS:-2}"
POSTPROCESS_OVERWRITE="${POSTPROCESS_OVERWRITE:-0}"
DELETE_FLEXPART_RUN_DIRS="${DELETE_FLEXPART_RUN_DIRS:-0}"
DELETE_DAILY_POSTPROCESS_FILES="${DELETE_DAILY_POSTPROCESS_FILES:-0}"
DRY_RUN="${DRY_RUN:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${FLEXPART_REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
POSTPROCESS_SCRIPT="${REPO_ROOT}/run_scripts_daily/postprocess_daily_to_monthly.py"

if [[ -z "${SITE}" ]]; then
  echo "ERROR: SITE is required. Submit with run_slurm_postprocess_daily_to_monthly_array.sh START_MONTH END_MONTH SITE DOMAIN."
  exit 2
fi

if [[ -z "${DOMAIN}" ]]; then
  echo "ERROR: DOMAIN is required. Submit with run_slurm_postprocess_daily_to_monthly_array.sh START_MONTH END_MONTH SITE DOMAIN."
  exit 2
fi

if [[ ! -f "${POSTPROCESS_SCRIPT}" ]]; then
  echo "ERROR: postprocess script not found: ${POSTPROCESS_SCRIPT}"
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

if ! [[ "${POSTPROCESS_WORKERS}" =~ ^[0-9]+$ ]] || [[ "${POSTPROCESS_WORKERS}" -lt 1 ]]; then
  echo "ERROR: POSTPROCESS_WORKERS must be an integer >= 1; got ${POSTPROCESS_WORKERS}."
  exit 2
fi

slurm_cpus_per_task="${SLURM_CPUS_PER_TASK:-2}"
if ! [[ "${slurm_cpus_per_task}" =~ ^[0-9]+$ ]] || [[ "${slurm_cpus_per_task}" -lt "${POSTPROCESS_WORKERS}" ]]; then
  echo "ERROR: SLURM_CPUS_PER_TASK must be >= POSTPROCESS_WORKERS (${POSTPROCESS_WORKERS}); got ${slurm_cpus_per_task}."
  exit 2
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
  exit 2
fi

# Slurm --export uses commas to separate variables, so the launcher exports
# months as whitespace-delimited values. Keep comma parsing for direct use.
months_raw="${MONTH_LIST}"
if [[ "${months_raw}" == *,* ]]; then
  IFS=',' read -r -a months <<< "${months_raw}"
else
  read -r -a months <<< "${months_raw}"
fi
count="${#months[@]}"
idx="${SLURM_ARRAY_TASK_ID:-1}"

if [[ ! "${idx}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID must be numeric; got ${idx}"
  exit 2
fi

if (( idx < 1 || idx > count )); then
  echo "INFO: SLURM_ARRAY_TASK_ID=${idx} outside month range 1-${count}; skipping."
  exit 0
fi

target_month="${months[$((idx - 1))]}"
if [[ ! "${target_month}" =~ ^[0-9]{6}$ ]]; then
  echo "ERROR: invalid month token at index ${idx}: ${target_month}"
  exit 2
fi

echo "Root dir: ${OUTROOT}"
echo "Site: ${SITE}"
echo "Domain: ${DOMAIN}"
echo "Target month: ${target_month} (${idx}/${count})"
echo "Postprocess workers: ${POSTPROCESS_WORKERS}"
echo "Delete FLEXPART run dirs after verification: ${DELETE_FLEXPART_RUN_DIRS}"
echo "Delete daily postprocess files after verification: ${DELETE_DAILY_POSTPROCESS_FILES}"

cmd=("${POSTPROCESS_DRIVER_PYTHON}" "${POSTPROCESS_SCRIPT}"
  --root-dir "${OUTROOT}"
  --site "${SITE}"
  --domain "${DOMAIN}"
  --month "${target_month}"
  --workers "${POSTPROCESS_WORKERS}"
  --python "${POSTPROCESS_PYTHON_CMD}"
  --include-exits
  --no-exit-csv
  --complete-months-only)

if [[ "${POSTPROCESS_OVERWRITE}" == "1" ]]; then
  cmd+=(--overwrite)
fi
if [[ "${DRY_RUN}" == "1" ]]; then
  cmd+=(--dry-run)
fi

"${cmd[@]}"

if [[ "${DELETE_FLEXPART_RUN_DIRS}" != "1" && "${DELETE_DAILY_POSTPROCESS_FILES}" != "1" ]]; then
  echo "Cleanup switches are off; keeping FLEXPART run directories and daily postprocess files."
  exit 0
fi

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1; skipping verification cleanup."
  exit 0
fi

"${POSTPROCESS_PYTHON_CMD}" - \
  "${OUTROOT}" \
  "${SITE}" \
  "${DOMAIN}" \
  "${target_month}" \
  "${DELETE_FLEXPART_RUN_DIRS}" \
  "${DELETE_DAILY_POSTPROCESS_FILES}" <<'PY'
from __future__ import annotations

import calendar
import shutil
import sys
from pathlib import Path

import pandas as pd
import xarray as xr

root = Path(sys.argv[1]).resolve()
site = sys.argv[2]
domain = sys.argv[3]
month = sys.argv[4]
delete_run_dirs = sys.argv[5] == "1"
delete_daily_files = sys.argv[6] == "1"

year = int(month[:4])
month_num = int(month[4:6])
expected_days = calendar.monthrange(year, month_num)[1]
expected_times = expected_days * 24

daily_dir = root / f"{site}_daily"
monthly_dir = root / site
daily_files = sorted(daily_dir.glob(f"*-*_FLEXPART_*_{domain}_inert_{month}??.nc"))
monthly_files = sorted(monthly_dir.glob(f"*-*_FLEXPART_*_{domain}_inert_{month}.nc"))

if len(daily_files) != expected_days:
    raise SystemExit(
        f"Refusing cleanup: found {len(daily_files)} daily NetCDF files in {daily_dir}, "
        f"expected {expected_days} for {month}."
    )

missing_daily = []
bad_daily = []
expected_daily_dates = {f"{month}{day:02d}" for day in range(1, expected_days + 1)}
found_daily_dates = set()

for path in daily_files:
    day = path.stem[-8:]
    found_daily_dates.add(day)
    with xr.open_dataset(path, decode_times=False) as ds:
        ntime = int(ds.sizes.get("time", 0))
        if ntime != 24:
            bad_daily.append(f"{path.name}: time={ntime}, expected 24")

missing_daily = sorted(expected_daily_dates - found_daily_dates)
if missing_daily or bad_daily:
    details = []
    if missing_daily:
        details.append("missing daily dates: " + ", ".join(missing_daily))
    if bad_daily:
        details.append("bad daily time counts: " + "; ".join(bad_daily[:10]))
    raise SystemExit("Refusing cleanup: " + " | ".join(details))

if not monthly_files:
    raise SystemExit(f"Refusing cleanup: no monthly NetCDF found in {monthly_dir} for {month}.")

valid_monthly = []
bad_monthly = []
expected_index = pd.date_range(f"{year:04d}-{month_num:02d}-01T00:00:00", periods=expected_times, freq="h")
expected_np = expected_index.to_numpy(dtype="datetime64[ns]")

for path in monthly_files:
    with xr.open_dataset(path) as ds:
        ntime = int(ds.sizes.get("time", 0))
        if ntime != expected_times:
            bad_monthly.append(f"{path.name}: time={ntime}, expected {expected_times}")
            continue
        if "time" not in ds.coords:
            bad_monthly.append(f"{path.name}: missing time coordinate")
            continue
        actual = ds["time"].values.astype("datetime64[ns]")
        if actual.shape[0] != expected_np.shape[0] or (actual != expected_np).any():
            bad_monthly.append(f"{path.name}: time coordinate is not contiguous hourly for {month}")
            continue
        valid_monthly.append(path)

if not valid_monthly:
    raise SystemExit(
        "Refusing cleanup: no monthly file passed completeness checks. "
        + "; ".join(bad_monthly[:10])
    )

if delete_run_dirs:
    run_dirs = [root / f"{domain}_{site}_{day}" for day in sorted(expected_daily_dates)]
    missing_run_dirs = [path for path in run_dirs if not path.is_dir()]
    if missing_run_dirs:
        print(
            "WARNING: some run directories are already absent: "
            + ", ".join(path.name for path in missing_run_dirs[:10])
        )

    for path in run_dirs:
        if path.is_dir():
            shutil.rmtree(path)
            print(f"DELETED run directory {path}")

if delete_daily_files:
    for path in daily_files:
        path.unlink()
        print(f"DELETED daily postprocess file {path}")

print(f"Verified {len(daily_files)} daily files and {len(valid_monthly)} monthly file(s) for {domain}_{site} {month}.")
PY
