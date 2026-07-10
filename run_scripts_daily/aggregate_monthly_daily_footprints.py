#!/usr/bin/env python3
"""
Aggregate daily FLEXPART footprint files into monthly NetCDF files.

Daily footprint products are expected to contain all hourly releases for one
calendar day on their internal ``time`` dimension, with file names like:

  <site>-<height>magl_FLEXPART_CFSv2_<domain>_inert_<YYYYMMDD>.nc

The monthly output keeps the same prefix and replaces the day stamp with
``YYYYMM``.
"""

from __future__ import annotations

import argparse
import calendar
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr


RE_DAILY = re.compile(r"^(?P<prefix>.+)_(?P<date>\d{8})\.nc$")
RE_MONTH = re.compile(r"^\d{6}$")


def _discover_daily_files(input_dir: Path) -> list[Path]:
    files = set(input_dir.glob("*_FLEXPART_GFS_*_inert_????????.nc"))
    files.update(input_dir.glob("*_FLEXPART_CFSv2_*_inert_????????.nc"))
    return sorted(files)


def _group_by_month(paths: list[Path]) -> dict[tuple[str, str], list[tuple[str, Path]]]:
    groups: dict[tuple[str, str], list[tuple[str, Path]]] = defaultdict(list)
    for path in paths:
        match = RE_DAILY.match(path.name)
        if match is None:
            continue
        prefix = match.group("prefix")
        yyyymmdd = match.group("date")
        groups[(prefix, yyyymmdd[:6])].append((yyyymmdd, path))

    for key in groups:
        groups[key].sort(key=lambda item: item[0])
    return groups


def _normalize_daily_time(ds: xr.Dataset, path: Path) -> xr.Dataset:
    match = RE_DAILY.match(path.name)
    if match is None:
        raise ValueError(f"cannot infer date from daily file name: {path.name}")
    yyyymmdd = match.group("date")
    base_date = np.datetime64(
        f"{yyyymmdd[0:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}T00:00:00",
        "ns",
    )

    if "time" in ds and "time" not in ds.coords:
        ds = ds.set_coords("time")

    if "time" not in ds.dims:
        ds = ds.expand_dims(time=[base_date])
    else:
        # Daily FLEXPART release times can be relative to each run. Use the
        # filename date so monthly products have absolute hourly coordinates.
        ntime = ds.sizes["time"]
        time_values = base_date + np.arange(ntime).astype("timedelta64[h]")
        ds = ds.assign_coords(time=("time", time_values))

    return ds


def _high_compression_encoding(ds: xr.Dataset, compression_level: int) -> dict[str, dict]:
    encoding: dict[str, dict] = {}
    for name in ds.data_vars:
        encoding[name] = {
            "zlib": True,
            "complevel": int(compression_level),
            "shuffle": True,
        }
    return encoding


def _is_complete_month(yyyymm: str, entries: list[tuple[str, Path]]) -> bool:
    year = int(yyyymm[:4])
    month = int(yyyymm[4:6])
    expected_days = calendar.monthrange(year, month)[1]
    found_days = {yyyymmdd for yyyymmdd, _ in entries}
    return len(found_days) == expected_days


def _write_monthly(paths: list[Path], out_file: Path, compression_level: int) -> None:
    datasets: list[xr.Dataset] = []
    out: xr.Dataset | None = None
    try:
        for path in paths:
            ds = xr.open_dataset(path, decode_times=False)
            ds = _normalize_daily_time(ds, path)
            datasets.append(ds)

        if not datasets:
            return

        out = xr.concat(datasets, dim="time", join="override", coords="minimal")
        if "time" in out.coords:
            order = np.argsort(np.asarray(out["time"].values))
            out = out.isel(time=order)

        out.attrs["aggregation"] = "monthly_from_daily"
        out.attrs["source_daily_files"] = str(len(paths))
        encoding = _high_compression_encoding(out, compression_level=compression_level)
        out.to_netcdf(out_file, encoding=encoding)
    finally:
        if out is not None:
            out.close()
        for ds in datasets:
            ds.close()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate daily footprint NetCDF files into monthly files"
    )
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing daily footprint files")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where monthly files are written (default: --input-dir)",
    )
    parser.add_argument(
        "--month",
        type=str,
        default=None,
        help="Optional month filter in YYYYMM format; only this month is written",
    )
    parser.add_argument(
        "--complete-months-only",
        action="store_true",
        help="Write only months with every calendar day present",
    )
    parser.add_argument("--compression-level", type=int, default=9)
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing monthly files")
    args = parser.parse_args()

    target_month: str | None = None
    if args.month is not None:
        if RE_MONTH.match(args.month) is None:
            raise ValueError(f"--month must match YYYYMM, got: {args.month}")
        target_month = args.month

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {input_dir}")

    output_dir = (args.output_dir if args.output_dir is not None else input_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    daily_files = _discover_daily_files(input_dir)
    if not daily_files:
        print(f"No daily footprint files found in {input_dir}")
        return 0

    groups = _group_by_month(daily_files)
    if target_month is not None:
        groups = {key: entries for key, entries in groups.items() if key[1] == target_month}

    if args.complete_months_only:
        groups = {
            key: entries
            for key, entries in groups.items()
            if _is_complete_month(key[1], entries)
        }

    if not groups:
        if target_month is None:
            print(f"No monthly groups built from files in {input_dir}")
        else:
            print(f"No monthly groups for {target_month} built from files in {input_dir}")
        return 0

    print(f"Discovered {len(daily_files)} daily files in {input_dir}")
    print(f"Building {len(groups)} monthly files in {output_dir}")
    if target_month is not None:
        print(f"Applying month filter: {target_month}")

    wrote = 0
    skipped = 0
    for (prefix, yyyymm), entries in sorted(groups.items()):
        monthly_name = f"{prefix}_{yyyymm}.nc"
        out_file = output_dir / monthly_name

        if out_file.exists() and not args.overwrite:
            skipped += 1
            print(f"SKIP exists: {out_file}")
            continue

        paths = [path for _, path in entries]
        _write_monthly(paths, out_file, compression_level=args.compression_level)
        wrote += 1
        print(f"WROTE {out_file} from {len(paths)} daily files")

    print(f"Monthly summary: wrote={wrote}, skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
