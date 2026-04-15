#!/usr/bin/env python3
"""
Download historical GFS data from NOAA archives and generate FLEXPART AVAILABLE file.

This script fetches GFS forecast GRIB2 files from NOAA's archive servers and organizes
them in FLEXPART-compatible format (naming: GFyymmddhh, AVAILABLE metadata).

Supports three archive sources (auto-selected by default based on date):
  1. AWS S3 Registry (for >= 2021-01-01): NOAA GFS on AWS OpenData
     - 0.25 degree resolution (gfs_0p25)
     - Fast, reliable, no authentication required
     - Note: path format changed on 2021-03-23 (added /atmos/ subdirectory)
  2. NOAA NOMADS: Direct NOAA operational archives
     - 0.25 degree, recent data only (~30-40 days)
     - Alternative if AWS is unavailable
  3. NCEI archive (for < 2021-01-01): NOAA NCEI historical archive
     - 0.5 degree resolution (gfsanl_4 / gfs_4 / grid-004)
     - THREDDS fileServer for 2004-04 to 2020-05
     - Object Store for 2020-06 to 2020-12
     - THREDDS URL: https://www.ncei.noaa.gov/thredds/fileServer/model-gfs-g4-anl-files-old/
     - Object Store URL: https://www.ncei.noaa.gov/oa/prod-model/global-forecast-system/access/grid-004-0.5-degree/analysis/

By default (--source auto) the script automatically selects the appropriate source
and resolution for each timestamp based on the date.
"""

import argparse
import os
import sys
from datetime import datetime, timedelta
import urllib.request
import urllib.error
import re

try:
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    HAS_BOTO3 = True
except ImportError:
    HAS_BOTO3 = False

# Date from which 0.25° GFS is available on AWS S3 (noaa-gfs-bdp-pds)
GFS_025_START = datetime(2021, 1, 1)

# AWS S3 introduced /atmos/ subdirectory on this date; earlier cycles live directly under cycle/
GFS_S3_ATMOS_START = datetime(2021, 3, 23)

# NCEI THREDDS continuous archive end; Object Store takes over from this date
# THREDDS (model-gfs-g4-anl-files-old) covers up to and including 2020-05-15
# Object Store (grid-004-0.5-degree)    covers from 2020-05-16 through end of 2020
# THREDDS: gfsanl_4_YYYYMMDD_HHMM_FFF.grb2
# Object Store: gfs_4_YYYYMMDD_HHMM_FFF.grb2
NCEI_OBJSTORE_START = datetime(2020, 5, 16)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download historical GFS data from NOAA archives (auto-selects source by date).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-select source by date (default):
  #   >= 2021-01-01: AWS S3 0.25 deg
  #   2020-06-01 to 2020-12-31: NCEI Object Store 0.5 deg
  #   2004-04 to 2020-05-31: NCEI THREDDS 0.5 deg
  %(prog)s --start 2019010100 --end 2019010200 --outdir ./gfs_data --available ./AVAILABLE

  # Force AWS S3 for recent data
  %(prog)s --start 2026030100 --end 2026040200 --source aws --outdir ./gfs_data

  # Force NCEI for historical 0.5 degree data (routes internally to THREDDS or Object Store)
  %(prog)s --start 2020080100 --end 2020080200 --source ncei --outdir ./gfs_data

  # Dry-run to preview
  %(prog)s --start 2019010100 --end 2019020100 --dry-run

Data sources and coverage:
  auto   (default) >= 2021-01-01: AWS S3 0.25 deg
                   2020-06 to 2020-12: NCEI Object Store 0.5 deg
                   2004-04 to 2020-05: NCEI THREDDS 0.5 deg
  aws    AWS OpenData S3 bucket (noaa-gfs-bdp-pds), 0.25 deg, from 2021-01-01
  nomads NOAA NOMADS direct HTTP, 0.25 deg, recent data only (~30-40 days)
  ncei   NCEI archive (THREDDS <= 2020-05, Object Store >= 2020-06), 0.5 deg
""",
    )

    parser.add_argument(
        "--start",
        required=True,
        metavar="YYYYMMDDHH",
        help="Start valid time (UTC)",
    )
    parser.add_argument(
        "--end", required=True, metavar="YYYYMMDDHH", help="End valid time (UTC), inclusive"
    )
    parser.add_argument(
        "--outdir",
        default="./inputs",
        metavar="DIR",
        help="Output directory for meteorological files (default: %(default)s)",
    )
    parser.add_argument(
        "--available",
        default="./AVAILABLE",
        metavar="FILE",
        help="AVAILABLE file path (default: %(default)s)",
    )
    parser.add_argument(
        "--step-hours",
        type=int,
        default=3,
        metavar="N",
        help="Time step in hours (default: %(default)s)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions only, do not download/write",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        metavar="SEC",
        help="Download timeout in seconds (default: %(default)s)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        metavar="N",
        help="Number of retry attempts (default: %(default)s)",
    )
    parser.add_argument(
        "--source",
        choices=["auto", "aws", "nomads", "ncei"],
        default="auto",
        metavar="SOURCE",
        help="Data source: auto (default), aws, nomads, or ncei",
    )

    return parser.parse_args()


def validate_ymdh(ymdh_str):
    """Validate YYYYMMDDHH format."""
    if not re.match(r"^\d{10}$", ymdh_str):
        raise ValueError(f"Invalid YYYYMMDDHH format: {ymdh_str}")
    try:
        datetime.strptime(ymdh_str, "%Y%m%d%H")
        return True
    except ValueError:
        raise ValueError(f"Invalid date/time: {ymdh_str}")


def ymdh_to_datetime(ymdh_str):
    """Convert YYYYMMDDHH to datetime."""
    return datetime.strptime(ymdh_str, "%Y%m%d%H")


def datetime_to_ymdh(dt):
    """Convert datetime to YYYYMMDDHH."""
    return dt.strftime("%Y%m%d%H")


def datetime_to_yymmddhh(dt):
    """Convert datetime to yymmddhh."""
    return dt.strftime("%y%m%d%H")


def nearest_gfs_cycle(valid_dt):
    """
    Find nearest previous GFS cycle (00, 06, 12, 18 UTC) for a valid time.
    Returns (cycle_datetime, forecast_hour).
    """
    # Find the nearest previous 6-hour cycle
    cycle_hour = (valid_dt.hour // 6) * 6
    cycle_dt = valid_dt.replace(hour=cycle_hour, minute=0, second=0, microsecond=0)

    if cycle_dt > valid_dt:
        # Move back to previous cycle if forecast would be negative
        cycle_dt -= timedelta(hours=6)

    forecast_hour = int((valid_dt - cycle_dt).total_seconds() // 3600)
    return cycle_dt, forecast_hour


def construct_noaa_url(cycle_dt, forecast_hour, abbrev="gfs_0p25"):
    """
    Construct NOAA operational archive URL for GFS GRIB2 file.

    Pattern:
      https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/
      gfs.YYYYMMDD/HH/atmos/gfs.tHHz.pgrb2.0p25.fFFF

    Args:
      cycle_dt: datetime object for the GFS cycle
      forecast_hour: forecast hour (int)
      abbrev: abbreviation for file type (e.g., 'gfs_0p25', 'gfs_0p50')

    Returns:
      URL string
    """
    cycle_date = cycle_dt.strftime("%Y%m%d")
    cycle_hour = cycle_dt.strftime("%H")
    forecast_str = f"{forecast_hour:03d}"

    # Map abbreviation to resolution and file suffix
    if abbrev == "gfs_0p25":
        resolution = "0p25"
    elif abbrev == "gfs_0p50":
        resolution = "0p50"
    elif abbrev == "gfs_1p00":
        resolution = "1p00"
    else:
        # Default to 0.25°
        resolution = "0p25"

    url = (
        f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/"
        f"gfs.{cycle_date}/{cycle_hour}/atmos/"
        f"gfs.t{cycle_hour}z.pgrb2.{resolution}.f{forecast_str}"
    )
    return url


def construct_aws_s3_path(cycle_dt, forecast_hour, abbrev="gfs_0p25"):
    """
    Construct AWS S3 path for GFS GRIB2 file.

    Pattern:
      noaa-gfs-bdp-pds/gfs.YYYYMMDD/HH/atmos/gfs.tHHz.pgrb2.0p25.fFFF

    Args:
      cycle_dt: datetime object for the GFS cycle
      forecast_hour: forecast hour (int)
      abbrev: abbreviation for file type (e.g., 'gfs_0p25', 'gfs_0p50')

    Returns:
      S3 object key string (without bucket name)
    """
    cycle_date = cycle_dt.strftime("%Y%m%d")
    cycle_hour = cycle_dt.strftime("%H")
    forecast_str = f"{forecast_hour:03d}"

    # Map abbreviation to resolution
    if abbrev == "gfs_0p25":
        resolution = "0p25"
    elif abbrev == "gfs_0p50":
        resolution = "0p50"
    elif abbrev == "gfs_1p00":
        resolution = "1p00"
    else:
        resolution = "0p25"

    # AWS S3 introduced /atmos/ subdirectory on 2021-03-23; earlier data is under cycle/ directly
    if cycle_dt >= GFS_S3_ATMOS_START:
        path = (
            f"gfs.{cycle_date}/{cycle_hour}/atmos/"
            f"gfs.t{cycle_hour}z.pgrb2.{resolution}.f{forecast_str}"
        )
    else:
        path = (
            f"gfs.{cycle_date}/{cycle_hour}/"
            f"gfs.t{cycle_hour}z.pgrb2.{resolution}.f{forecast_str}"
        )
    return path


def construct_ncei_url(cycle_dt, forecast_hour):
    """
    Construct NCEI URL for GFS 0.5 degree GRIB2 file.

    Routes internally based on date:
      - <= 2020-05-31: NCEI THREDDS fileServer (model-gfs-g4-anl-files-old)
        Pattern: gfsanl_4_YYYYMMDD_HHMM_FFF.grb2
      - >= 2020-06-01: NCEI Object Store (grid-004-0.5-degree/analysis)
        Pattern: gfs_4_YYYYMMDD_HHMM_FFF.grb2

    Available cycles: 0000, 0600, 1200, 1800 UTC
    Available forecast hours: 000, 003, 006

    Args:
      cycle_dt: datetime object for the GFS cycle
      forecast_hour: forecast hour (int, 0-6)

    Returns:
      URL string
    """
    yyyymm = cycle_dt.strftime("%Y%m")
    yyyymmdd = cycle_dt.strftime("%Y%m%d")
    hhmm = cycle_dt.strftime("%H") + "00"  # e.g. "0000", "0600", "1200", "1800"
    forecast_str = f"{forecast_hour:03d}"

    if cycle_dt >= NCEI_OBJSTORE_START:
        # NCEI Object Store: gfs_4 naming (no 'anl'), covers 2020-05-16 to end of 2020
        url = (
            f"https://www.ncei.noaa.gov/oa/prod-model/global-forecast-system/access/"
            f"grid-004-0.5-degree/analysis/"
            f"{yyyymm}/{yyyymmdd}/gfs_4_{yyyymmdd}_{hhmm}_{forecast_str}.grb2"
        )
    else:
        # NCEI THREDDS: gfsanl_4 naming, covers 2004-04 to 2020-05-15
        url = (
            f"https://www.ncei.noaa.gov/thredds/fileServer/model-gfs-g4-anl-files-old/"
            f"{yyyymm}/{yyyymmdd}/gfsanl_4_{yyyymmdd}_{hhmm}_{forecast_str}.grb2"
        )
    return url


def select_source(valid_dt, requested_source):
    """
    Determine the actual data source to use for a given valid time.

    In 'auto' mode, selects based on date:
      - >= GFS_025_START (2021-01-01): AWS S3 (0.25 degree)
      - < GFS_025_START: NCEI (0.5 degree; internally routes to THREDDS or Object Store)

    Args:
      valid_dt: datetime for the requested valid time
      requested_source: source string from CLI ('auto', 'aws', 'nomads', 'ncei')

    Returns:
      Effective source string: 'aws', 'nomads', or 'ncei'
    """
    if requested_source != "auto":
        return requested_source
    return "aws" if valid_dt >= GFS_025_START else "ncei"


def download_file(url, dst_path, timeout=60, retries=3):
    """
    Download a file from URL to destination path with retries.

    Args:
      url: Source URL
      dst_path: Destination file path
      timeout: Timeout in seconds
      retries: Number of retry attempts

    Returns:
      True if successful, False otherwise
    """
    for attempt in range(retries):
        try:
            print(f"  Downloading (attempt {attempt + 1}/{retries})...", file=sys.stderr)
            req = urllib.request.urlopen(url, timeout=timeout)
            with open(dst_path, "wb") as f:
                f.write(req.read())
            return True
        except urllib.error.URLError as e:
            if attempt == retries - 1:
                print(f"  Error: {e}", file=sys.stderr)
                return False
            print(f"  Retry after error: {e}", file=sys.stderr)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"  Error: File not found (404)", file=sys.stderr)
                return False
            if attempt == retries - 1:
                print(f"  Error (HTTP {e.code}): {e}", file=sys.stderr)
                return False
            print(f"  Retry after HTTP {e.code}", file=sys.stderr)
        except Exception as e:
            print(f"  Error: {e}", file=sys.stderr)
            return False

    return False


def download_from_s3(s3_key, dst_path, bucket="noaa-gfs-bdp-pds", retries=3):
    """
    Download a file from AWS S3 bucket to destination path with retries.

    Args:
      s3_key: S3 object key (path within bucket)
      dst_path: Destination file path
      bucket: S3 bucket name
      retries: Number of retry attempts

    Returns:
      True if successful, False otherwise
    """
    if not HAS_BOTO3:
        print("Error: boto3 not available. Install with: pip install boto3", file=sys.stderr)
        return False

    for attempt in range(retries):
        try:
            print(f"  Downloading from S3 (attempt {attempt + 1}/{retries})...", file=sys.stderr)
            s3client = boto3.client(
                "s3",
                config=Config(signature_version=UNSIGNED),
                region_name="us-east-1"
            )
            s3client.download_file(bucket, s3_key, dst_path)
            return True
        except Exception as e:
            if attempt == retries - 1:
                print(f"  Error: {e}", file=sys.stderr)
                return False
            print(f"  Retry after error: {e}", file=sys.stderr)

    return False


def generate_available_header():
    """Generate the header lines for AVAILABLE file."""
    return ["XXXXXX EMPTY LINES XXXXXXXXX\n", "XXXXXX EMPTY LINES XXXXXXXX\n"]


def write_available_line(ymdh_str, short_yymmddhh, filename):
    """
    Format a line for AVAILABLE file.

    Format: YYYYMMDD HHMMSS   filename
    """
    yyyy_mm_dd = ymdh_str[:8]
    hh_mmss = ymdh_str[8:10] + "0000"
    return f"{yyyy_mm_dd} {hh_mmss}      {filename:<12s} ON DISK\n"


def main():
    args = parse_args()

    # Validate inputs
    try:
        validate_ymdh(args.start)
        validate_ymdh(args.end)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.step_hours <= 0:
        print("Error: --step-hours must be positive", file=sys.stderr)
        sys.exit(1)

    # Check boto3 availability when AWS could be used
    if args.source in ("aws", "auto") and not HAS_BOTO3:
        if args.source == "aws":
            print("Warning: boto3 not available. Install with: pip install boto3", file=sys.stderr)
            print("Falling back to NOMADS source.", file=sys.stderr)
            args.source = "nomads"
        else:
            print("Warning: boto3 not available; 'auto' will use NOMADS instead of AWS S3 for recent dates.", file=sys.stderr)
            print("Install boto3 with: pip install boto3", file=sys.stderr)

    start_dt = ymdh_to_datetime(args.start)
    end_dt = ymdh_to_datetime(args.end)

    if start_dt > end_dt:
        print("Error: --start must be <= --end", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    if not args.dry_run:
        os.makedirs(args.outdir, exist_ok=True)

    # Prepare AVAILABLE file
    available_lines = generate_available_header()
    available_lines.append("YYYYMMDD HHMMSS   name of the file(up to 80 characters)\n")

    source_label = {"auto": "auto (date-based)", "aws": "AWS S3", "nomads": "NOMADS", "ncei": "NCEI THREDDS"}
    print(f"Using {source_label.get(args.source, args.source)} as data source")
    if args.source == "auto":
        print(f"  >= {GFS_025_START.strftime('%Y-%m-%d')}: AWS S3 (0.25 deg)")
        print(f"  {NCEI_OBJSTORE_START.strftime('%Y-%m-%d')} to {GFS_025_START.strftime('%Y-%m-%d')}: NCEI Object Store (0.5 deg)")
        print(f"  < {NCEI_OBJSTORE_START.strftime('%Y-%m-%d')}: NCEI THREDDS (0.5 deg)")

    count = 0
    current_dt = start_dt

    while current_dt <= end_dt:
        ymdh = datetime_to_ymdh(current_dt)
        yymmddhh = datetime_to_yymmddhh(current_dt)
        out_name = f"GF{yymmddhh}"

        # Find nearest GFS cycle and forecast hour
        cycle_dt, forecast_hour = nearest_gfs_cycle(current_dt)

        # Resolve effective source (may differ per timestamp in auto mode)
        effective_source = select_source(current_dt, args.source)

        # NCEI only has forecast hours 0-6; skip if step puts us beyond that
        if effective_source == "ncei" and forecast_hour > 6:
            print(f"Warning: {ymdh} requires forecast +{forecast_hour:03d}h which is not available in NCEI (max +006h). Skipping.", file=sys.stderr)
            current_dt += timedelta(hours=args.step_hours)
            continue

        # Destination path
        dst_path = os.path.join(args.outdir, out_name)

        # Check if file exists and --force not set
        if os.path.isfile(dst_path) and not args.force:
            print(f"Skip existing: {dst_path}")
        else:
            res_label = "0.5 deg" if effective_source == "ncei" else "0.25 deg"
            print(f"Fetch {ymdh} -> {out_name} [{res_label} via {effective_source}]")
            print(f"  Cycle: {datetime_to_ymdh(cycle_dt)}, Forecast: +{forecast_hour:03d}h")

            if effective_source == "aws":
                s3_key = construct_aws_s3_path(cycle_dt, forecast_hour, abbrev="gfs_0p25")
                print(f"  S3: noaa-gfs-bdp-pds/{s3_key}")
            elif effective_source == "ncei":
                url = construct_ncei_url(cycle_dt, forecast_hour)
                print(f"  URL: {url}")
            else:  # nomads
                url = construct_noaa_url(cycle_dt, forecast_hour, abbrev="gfs_0p25")
                print(f"  URL: {url}")

            if not args.dry_run:
                success = False
                if effective_source == "aws":
                    s3_key = construct_aws_s3_path(cycle_dt, forecast_hour, abbrev="gfs_0p25")
                    success = download_from_s3(s3_key, dst_path, retries=args.retries)
                elif effective_source == "ncei":
                    url = construct_ncei_url(cycle_dt, forecast_hour)
                    success = download_file(url, dst_path, timeout=args.timeout, retries=args.retries)
                else:  # nomads
                    url = construct_noaa_url(cycle_dt, forecast_hour, abbrev="gfs_0p25")
                    success = download_file(url, dst_path, timeout=args.timeout, retries=args.retries)

                if success:
                    print(f"  SUCCESS: {out_name}")
                else:
                    print(f"  FAILED: {out_name} (skipping AVAILABLE entry)")
                    current_dt += timedelta(hours=args.step_hours)
                    continue

        # Add to AVAILABLE
        available_lines.append(write_available_line(ymdh, yymmddhh, out_name))

        count += 1
        current_dt += timedelta(hours=args.step_hours)

    # Write AVAILABLE file
    if not args.dry_run:
        try:
            with open(args.available, "w") as f:
                f.writelines(available_lines)
            print(f"\nWrote {args.available} with {count} entries.")
        except IOError as e:
            print(f"Error writing {args.available}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"\nDry run complete ({count} timestamps).")


if __name__ == "__main__":
    main()
