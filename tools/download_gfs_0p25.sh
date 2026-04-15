#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Download GFS 0.25 degree files from NOMADS and generate FLEXPART AVAILABLE.

Usage:
  ./tools/download_gfs_0p25.sh --start YYYYMMDDHH --end YYYYMMDDHH [options]

Required:
  --start YYYYMMDDHH   Start valid time (UTC)
  --end YYYYMMDDHH     End valid time (UTC), inclusive

Options:
  --outdir DIR         Output directory for meteorological files (default: ./inputs)
  --available FILE     AVAILABLE file path (default: ./AVAILABLE)
  --step-hours N       Time step in hours (default: 3)
  --force              Re-download files even if they already exist
  --dry-run            Print actions only, do not download/write
  -h, --help           Show this help

Notes:
  - Files are renamed to FLEXPART-style names: GFyymmddhh
  - The script picks the nearest previous 6-hour GFS cycle for each valid time.
  - Source URL pattern:
    https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.YYYYMMDD/HH/atmos/gfs.tHHz.pgrb2.0p25.fFFF
EOF
}

START=""
END=""
OUTDIR="/net/fs01/data/AGAGE/meteorology/gfs_grib/"
AVAILABLE_FILE="./AVAILABLE"
STEP_HOURS=3
FORCE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start)
      START="$2"
      shift 2
      ;;
    --end)
      END="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --available)
      AVAILABLE_FILE="$2"
      shift 2
      ;;
    --step-hours)
      STEP_HOURS="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$START" || -z "$END" ]]; then
  echo "Error: --start and --end are required." >&2
  usage
  exit 1
fi

if ! [[ "$START" =~ ^[0-9]{10}$ && "$END" =~ ^[0-9]{10}$ ]]; then
  echo "Error: start/end must be YYYYMMDDHH." >&2
  exit 1
fi

if ! [[ "$STEP_HOURS" =~ ^[0-9]+$ ]] || [[ "$STEP_HOURS" -le 0 ]]; then
  echo "Error: --step-hours must be a positive integer." >&2
  exit 1
fi

to_epoch() {
  local ymdh="$1"
  date -u -d "${ymdh:0:8} ${ymdh:8:2}:00:00" +%s
}

from_epoch_ymdh() {
  local epoch="$1"
  date -u -d "@${epoch}" +%Y%m%d%H
}

from_epoch_yymmddhh() {
  local epoch="$1"
  date -u -d "@${epoch}" +%y%m%d%H
}

START_EPOCH="$(to_epoch "$START")"
END_EPOCH="$(to_epoch "$END")"
STEP_SEC=$((STEP_HOURS * 3600))

if [[ "$START_EPOCH" -gt "$END_EPOCH" ]]; then
  echo "Error: --start must be <= --end." >&2
  exit 1
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
  mkdir -p "$OUTDIR"
fi

if command -v curl >/dev/null 2>&1; then
  DOWNLOADER="curl"
elif command -v wget >/dev/null 2>&1; then
  DOWNLOADER="wget"
else
  echo "Error: neither curl nor wget is available." >&2
  exit 1
fi

AVAILABLE_TMP="${AVAILABLE_FILE}.tmp"
if [[ "$DRY_RUN" -eq 0 ]]; then
  {
    echo "XXXXXX EMPTY LINES XXXXXXXXX"
    echo "XXXXXX EMPTY LINES XXXXXXXX"
    echo "YYYYMMDD HHMMSS   name of the file(up to 80 characters)"
  } > "$AVAILABLE_TMP"
fi

cur="$START_EPOCH"
count=0
while [[ "$cur" -le "$END_EPOCH" ]]; do
  valid_ymdh="$(from_epoch_ymdh "$cur")"
  valid_short="$(from_epoch_yymmddhh "$cur")"
  out_name="GF${valid_short}"

  # Nearest previous 6-hour cycle for this valid time.
  cycle_epoch=$((cur - (cur % 21600)))
  cycle_ymdh="$(from_epoch_ymdh "$cycle_epoch")"
  cycle_date="${cycle_ymdh:0:8}"
  cycle_hour="${cycle_ymdh:8:2}"

  fh=$(((cur - cycle_epoch) / 3600))
  printf -v fff "%03d" "$fh"

  src_url="https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${cycle_date}/${cycle_hour}/atmos/gfs.t${cycle_hour}z.pgrb2.0p25.f${fff}"
  dst_path="${OUTDIR}/${out_name}"

  if [[ -f "$dst_path" && "$FORCE" -eq 0 ]]; then
    echo "Skip existing: ${dst_path}"
  else
    echo "Fetch ${valid_ymdh} -> ${out_name}"
    echo "  ${src_url}"
    if [[ "$DRY_RUN" -eq 0 ]]; then
      if [[ "$DOWNLOADER" == "curl" ]]; then
        curl -fL --retry 4 --retry-delay 3 --connect-timeout 20 --max-time 300 -o "$dst_path" "$src_url"
      else
        wget -O "$dst_path" "$src_url"
      fi
    fi
  fi

  if [[ "$DRY_RUN" -eq 0 ]]; then
    printf "%s %s      %-12s ON DISK\n" "${valid_ymdh:0:8}" "${valid_ymdh:8:2}0000" "$out_name" >> "$AVAILABLE_TMP"
  fi

  count=$((count + 1))
  cur=$((cur + STEP_SEC))
done

if [[ "$DRY_RUN" -eq 0 ]]; then
  mv "$AVAILABLE_TMP" "$AVAILABLE_FILE"
  echo "Wrote ${AVAILABLE_FILE} with ${count} entries."
else
  echo "Dry run complete (${count} timestamps)."
fi
