#!/usr/bin/env python3
"""
Aggregate hourly FLEXPART footprint files into monthly NetCDF files.

Input file naming convention expected:
  <prefix>_<YYYYMMDDHH>.nc
where <prefix> is the full descriptor before the timestamp.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import xarray as xr


RE_HOURLY = re.compile(r"^(?P<prefix>.+)_(?P<ts>\d{10})\.nc$")


def _discover_hourly_files(input_dir: Path) -> list[Path]:
    return sorted(input_dir.glob("*_FLEXPART_GFS_*_inert_??????????.nc"))


def _group_by_month(paths: list[Path]) -> dict[tuple[str, str], list[tuple[str, Path]]]:
    groups: dict[tuple[str, str], list[tuple[str, Path]]] = defaultdict(list)
    for path in paths:
        m = RE_HOURLY.match(path.name)
        if not m:
            continue
        prefix = m.group("prefix")
        ts10 = m.group("ts")
        yyyymm = ts10[:6]
        groups[(prefix, yyyymm)].append((ts10, path))

    for key in groups:
        groups[key].sort(key=lambda item: item[0])
    return groups


def _write_monthly(paths: list[Path], out_file: Path) -> None:
    datasets: list[xr.Dataset] = []
    try:
        for path in paths:
            ds = xr.open_dataset(path)
            datasets.append(ds)

        if not datasets:
            return

        out = xr.concat(datasets, dim="time")
        if "time" in out.coords:
            order = np.argsort(np.asarray(out["time"].values))
            out = out.isel(time=order)

        out.attrs["aggregation"] = "monthly"
        out.attrs["source_files"] = str(len(paths))
        out.to_netcdf(out_file)
    finally:
        for ds in datasets:
            ds.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate hourly footprint NetCDF files into monthly files")
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing hourly footprint NetCDF files")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where monthly files are written (default: --input-dir)",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing monthly files")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {input_dir}")

    output_dir = (args.output_dir if args.output_dir is not None else input_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    hourly_files = _discover_hourly_files(input_dir)
    if not hourly_files:
        print(f"No hourly footprint files found in {input_dir}")
        return 0

    groups = _group_by_month(hourly_files)
    if not groups:
        print(f"No monthly groups built from files in {input_dir}")
        return 0

    print(f"Discovered {len(hourly_files)} hourly files in {input_dir}")
    print(f"Building {len(groups)} monthly files in {output_dir}")

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
        _write_monthly(paths, out_file)
        wrote += 1
        print(f"WROTE {out_file} from {len(paths)} hourly files")

    print(f"Monthly summary: wrote={wrote}, skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
