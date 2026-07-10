#!/usr/bin/env python3
"""
Postprocess daily FLEXPART outputs and write monthly footprint files.

This is the daily-run analogue of the hourly postprocess-all workflow. It scans
daily run directories, writes one daily footprint file containing all 24 hourly
footprints, then aggregates those daily files into monthly files. By default the
final products are written into a folder named after the site/receptor, e.g. THD.

Example usage:
python run_scripts_daily/postprocess_daily_to_monthly.py \
  --root-dir /net/fs06/d2/${USER}/flexpart_outs_daily \
  --site THD --domain WESTUSA --month 201801 \
  --workers 8 --include-exits --no-exit-csv

"""

from __future__ import annotations

import argparse
import calendar
import concurrent.futures
import datetime as dt
import os
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "run_scripts"))

from postprocess_footprint import DEFAULT_CFSV2_GRIB_DIR  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
POSTPROCESS_SCRIPT = SCRIPT_DIR / "postprocess_daily_footprints.py"
MONTHLY_SCRIPT = SCRIPT_DIR / "aggregate_monthly_daily_footprints.py"
RE_RUN_DIR = re.compile(r"^(.+)_([^_]+)_(\d{8})$")


def _run(cmd: list[str], dry_run: bool) -> int:
    print("RUN " + " ".join(str(part) for part in cmd))
    if dry_run:
        return 0

    env = os.environ.copy()
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    result = subprocess.run(cmd, env=env)
    return int(result.returncode)


def _parse_run_dir(path: Path) -> tuple[str, str, dt.date] | None:
    match = RE_RUN_DIR.match(path.name)
    if match is None:
        return None
    domain, site, yyyymmdd = match.groups()
    try:
        run_date = dt.datetime.strptime(yyyymmdd, "%Y%m%d").date()
    except ValueError:
        return None
    return domain, site, run_date


def _discover_grid_files(
    root_dir: Path,
    site: str,
    domain: str | None,
    month: str | None,
    complete_months_only: bool,
) -> list[Path]:
    candidates: list[tuple[str, str, dt.date, Path]] = []
    for grid_file in sorted(root_dir.glob("**/output/grid_time_*.nc")):
        parsed = _parse_run_dir(grid_file.parent.parent)
        if parsed is None:
            continue
        run_domain, run_site, run_date = parsed
        if run_site != site:
            continue
        if domain is not None and run_domain != domain:
            continue
        if month is not None and run_date.strftime("%Y%m") != month:
            continue
        candidates.append((run_domain, run_site, run_date, grid_file))

    if complete_months_only:
        month_days: dict[tuple[str, str, str], set[dt.date]] = {}
        for run_domain, run_site, run_date, _ in candidates:
            key = (run_domain, run_site, run_date.strftime("%Y%m"))
            month_days.setdefault(key, set()).add(run_date)

        complete_keys = set()
        for key, dates in month_days.items():
            yyyymm = key[2]
            expected = calendar.monthrange(int(yyyymm[:4]), int(yyyymm[4:6]))[1]
            if len(dates) == expected:
                complete_keys.add(key)
            else:
                print(
                    f"SKIP incomplete month {key[0]}_{key[1]} {yyyymm}: "
                    f"{len(dates)} days found, expected {expected}"
                )

        candidates = [
            item
            for item in candidates
            if (item[0], item[1], item[2].strftime("%Y%m")) in complete_keys
        ]

    return [grid_file for _, _, _, grid_file in candidates]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Postprocess daily grid_time outputs and aggregate monthly footprints."
    )
    parser.add_argument(
        "--root-dir",
        required=True,
        type=Path,
        help="Root directory containing daily run directories like DOMAIN_SITE_YYYYMMDD",
    )
    parser.add_argument("--site", required=True, help="Site/receptor code, e.g. THD")
    parser.add_argument(
        "--domain",
        default=None,
        help="Optional domain filter/name; inferred by postprocess script when omitted",
    )
    parser.add_argument(
        "--final-dir",
        type=Path,
        default=None,
        help="Directory for daily and monthly products (default: --root-dir/--site)",
    )
    parser.add_argument(
        "--daily-dir",
        type=Path,
        default=None,
        help="Directory for daily footprint files (default: --final-dir)",
    )
    parser.add_argument(
        "--monthly-dir",
        type=Path,
        default=None,
        help="Directory for monthly footprint files (default: --final-dir)",
    )
    parser.add_argument(
        "--month",
        type=str,
        default=None,
        help="Optional month filter in YYYYMM format",
    )
    parser.add_argument("--python", default=sys.executable, help="Python executable for subprocesses")
    parser.add_argument("--postprocess-script", type=Path, default=POSTPROCESS_SCRIPT)
    parser.add_argument("--monthly-script", type=Path, default=MONTHLY_SCRIPT)
    parser.add_argument("--species", default="inert")
    parser.add_argument("--model", default="FLEXPART")
    parser.add_argument("--met-model", default="CFSv2")
    parser.add_argument(
        "--meteo-grib-dir",
        type=Path,
        default=Path(DEFAULT_CFSV2_GRIB_DIR),
        help="Directory containing GFyymmddhh GRIB files for release-time meteorology extraction.",
    )
    parser.add_argument(
        "--meteo-grib-file",
        type=Path,
        default=None,
        help="Optional explicit GRIB file path to use for all release-time meteorology extraction.",
    )
    parser.add_argument(
        "--disable-grib-meteo",
        action="store_true",
        help="Disable GRIB-based extraction of release-time meteorology.",
    )
    parser.add_argument("--source-layer-thickness-m", type=float, default=100.0)
    parser.add_argument("--footprint-outheight-m", type=float, default=100.0)
    parser.add_argument("--compression-level", type=int, default=4)
    parser.add_argument(
        "--include-exits",
        action="store_true",
        help="Add domain-exit diagnostics from partoutput_*.nc beside each daily grid file.",
    )
    parser.add_argument(
        "--no-exit-csv",
        action="store_true",
        help="When --include-exits is set, skip writing individual daily exit-point CSV files.",
    )
    parser.add_argument(
        "--monthly-compression-level",
        type=int,
        default=4,
        help="NetCDF compression level for monthly files",
    )
    parser.add_argument(
        "--complete-months-only",
        action="store_true",
        help="Process and aggregate only months with every calendar day present",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int, default=0, help="Optional max daily grid files to process")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers for daily postprocessing (default: 1)",
    )
    args = parser.parse_args()

    root_dir = args.root_dir.resolve()
    if not root_dir.is_dir():
        raise FileNotFoundError(f"root directory not found: {root_dir}")
    if not args.postprocess_script.exists():
        raise FileNotFoundError(f"postprocess script not found: {args.postprocess_script}")
    if not args.monthly_script.exists():
        raise FileNotFoundError(f"monthly script not found: {args.monthly_script}")

    final_dir = (args.final_dir if args.final_dir is not None else root_dir / args.site).resolve()
    if args.final_dir is not None:
        # If --final-dir is explicitly set, use it as default for both daily and monthly (backward compatible)
        daily_dir = (args.daily_dir if args.daily_dir is not None else final_dir).resolve()
        monthly_dir = (args.monthly_dir if args.monthly_dir is not None else final_dir).resolve()
    else:
        # If --final-dir is not set, use sensible defaults: {site}_daily for daily, {site} for monthly
        daily_dir = (args.daily_dir if args.daily_dir is not None else root_dir / f"{args.site}_daily").resolve()
        monthly_dir = (args.monthly_dir if args.monthly_dir is not None else final_dir).resolve()

    if not args.dry_run:
        daily_dir.mkdir(parents=True, exist_ok=True)
        monthly_dir.mkdir(parents=True, exist_ok=True)

    grid_files = _discover_grid_files(
        root_dir=root_dir,
        site=args.site,
        domain=args.domain,
        month=args.month,
        complete_months_only=args.complete_months_only,
    )
    if args.limit > 0:
        grid_files = grid_files[: args.limit]

    if not grid_files:
        print(f"No matching daily grid_time files found under {root_dir}")
    else:
        print(f"Discovered {len(grid_files)} matching daily grid_time file(s)")

    def _build_postprocess_cmd(grid_file: Path) -> list[str]:
        cmd = [
            str(args.python),
            str(args.postprocess_script),
            "--grid-file",
            str(grid_file),
            "--output-dir",
            str(daily_dir),
            "--site",
            args.site,
            "--species",
            args.species,
            "--model",
            args.model,
            "--met-model",
            args.met_model,
            "--meteo-grib-dir",
            str(args.meteo_grib_dir),
            "--source-layer-thickness-m",
            str(args.source_layer_thickness_m),
            "--footprint-outheight-m",
            str(args.footprint_outheight_m),
            "--compression-level",
            str(args.compression_level),
        ]
        if args.domain is not None:
            cmd.extend(["--domain", args.domain])
        if args.meteo_grib_file is not None:
            cmd.extend(["--meteo-grib-file", str(args.meteo_grib_file)])
        if args.disable_grib_meteo:
            cmd.append("--disable-grib-meteo")
        if args.overwrite:
            cmd.append("--overwrite")
        if args.include_exits:
            cmd.append("--auto-partoutput")
        if args.no_exit_csv:
            cmd.append("--no-exit-csv")
        if args.dry_run:
            cmd.append("--dry-run")
        return cmd

    fail = 0
    workers = max(1, args.workers)
    if workers == 1:
        for grid_file in grid_files:
            rc = _run(_build_postprocess_cmd(grid_file), dry_run=False)
            if rc != 0:
                fail += 1
    else:
        print(f"Running daily postprocess with {workers} parallel workers")
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_run, _build_postprocess_cmd(gf), False): gf
                for gf in grid_files
            }
            for future in concurrent.futures.as_completed(futures):
                rc = future.result()
                if rc != 0:
                    fail += 1

    if fail:
        print(f"Daily postprocess failures: {fail}")
        return 1

    monthly_cmd = [
        str(args.python),
        str(args.monthly_script),
        "--input-dir",
        str(daily_dir),
        "--output-dir",
        str(monthly_dir),
        "--compression-level",
        str(args.monthly_compression_level),
    ]
    if args.month is not None:
        monthly_cmd.extend(["--month", args.month])
    if args.complete_months_only:
        monthly_cmd.append("--complete-months-only")
    if args.overwrite:
        monthly_cmd.append("--overwrite")

    return _run(monthly_cmd, dry_run=args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())
