#!/usr/bin/env python3
"""
Postprocess all FLEXPART grid_time files discovered under a directory tree.

This script is intentionally lightweight (stdlib only) and delegates actual
NetCDF/xarray work to postprocess_footprint.py via subprocess, so you can run
it with any Python that has xarray/netCDF dependencies working.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


RE_GRID_TIME = re.compile(r"^grid_time_(\d{14})\.nc$")


def _run_postprocess_cmd(cmd: list[str]) -> int:
    env = os.environ.copy()
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    result = subprocess.run(cmd, env=env)
    return int(result.returncode)


def _load_release_height_map(site_info_path: Path) -> dict[str, float]:
    if not site_info_path.exists():
        return {}
    with open(site_info_path) as f:
        data = json.load(f)

    out: dict[str, float] = {}
    for code, location_data in data.items():
        if not isinstance(location_data, dict):
            continue
        station_info = None
        for _, network_data in location_data.items():
            if isinstance(network_data, dict) and "height" in network_data:
                station_info = network_data
                break
        if station_info is None:
            continue

        release_height = None
        for val in station_info.get("height", []):
            m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*m\s*$", str(val), flags=re.IGNORECASE)
            if m:
                release_height = float(m.group(1))
                break
        if release_height is None:
            release_height = 10.0
        out[str(code)] = float(release_height)

    return out


def _parse_run_name(run_dir_name: str) -> tuple[str, str, str] | None:
    # Expect: <domain>_<receptor>_<YYYYMMDDHH>; domain may include underscores.
    parts = run_dir_name.rsplit("_", 2)
    if len(parts) != 3:
        return None
    domain, receptor, yyyymmddhh = parts
    if not re.match(r"^\d{10}$", yyyymmddhh):
        return None
    return domain, receptor, yyyymmddhh


def _format_magl_label(height: float) -> str:
    if float(height).is_integer():
        return str(int(height))
    return f"{float(height):.3f}".rstrip("0").rstrip(".")


def _discover_grid_files(root_dir: Path) -> list[Path]:
    return sorted(root_dir.glob("**/output/grid_time_*.nc"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Postprocess all FLEXPART grid_time outputs under a root directory")
    parser.add_argument("--root-dir", required=True, type=Path, help="Root directory containing FLEXPART run folders")
    parser.add_argument("--python", default=sys.executable, help="Python executable used to run postprocess_footprint.py")
    parser.add_argument("--postprocess-script", type=Path, default=Path(__file__).parent / "postprocess_footprint.py")
    parser.add_argument("--site-info", type=Path, default=Path(__file__).parent.parent / "site_domains" / "site_info.json")
    parser.add_argument(
        "--postprocess-footprint-outheight-m",
        type=float,
        default=100.0,
        help="Outheight layer (m agl) used to derive footprints (default: 100).",
    )
    parser.add_argument(
        "--postprocess-lowest-magl",
        type=float,
        default=None,
        help="Deprecated alias for --postprocess-footprint-outheight-m.",
    )
    parser.add_argument("--postprocess-source-layer-thickness-m", type=float, default=100.0)
    parser.add_argument(
        "--final-dir",
        type=Path,
        default=None,
        help="Directory where final footprint NetCDF files are stored (default: --root-dir)",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing postprocessed footprint files")
    parser.add_argument(
        "--keep-run-dirs",
        action="store_true",
        help="Keep original FLEXPART run directories after successful postprocessing",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands but do not execute")
    parser.add_argument("--limit", type=int, default=0, help="Optional max number of grid files to process (0 = all)")
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of concurrent postprocess workers for NetCDF runs (default: 8)",
    )
    parser.add_argument(
        "--write-monthly",
        action="store_true",
        help="After per-hour postprocessing, write monthly aggregated NetCDF files",
    )
    parser.add_argument(
        "--monthly-script",
        type=Path,
        default=Path(__file__).parent / "aggregate_monthly_footprints.py",
        help="Script used to aggregate hourly files into monthly files",
    )
    parser.add_argument(
        "--monthly-dir",
        type=Path,
        default=None,
        help="Directory where monthly aggregated files are written (default: --final-dir)",
    )
    args = parser.parse_args()

    workers = max(1, int(args.workers))

    footprint_outheight_m = args.postprocess_footprint_outheight_m
    if args.postprocess_lowest_magl is not None:
        footprint_outheight_m = args.postprocess_lowest_magl
        print(
            "WARNING: --postprocess-lowest-magl is deprecated; use "
            "--postprocess-footprint-outheight-m instead."
        )

    root_dir = args.root_dir
    if not root_dir.is_dir():
        raise FileNotFoundError(f"root directory not found: {root_dir}")
    if not args.postprocess_script.exists():
        raise FileNotFoundError(f"postprocess script not found: {args.postprocess_script}")

    final_dir = args.final_dir if args.final_dir is not None else root_dir
    final_dir = final_dir.resolve()
    if not final_dir.exists() and not args.dry_run:
        final_dir.mkdir(parents=True, exist_ok=True)

    release_heights = _load_release_height_map(args.site_info)
    grid_files = _discover_grid_files(root_dir)
    if args.limit > 0:
        grid_files = grid_files[: args.limit]

    if not grid_files:
        print(f"No grid_time_*.nc files found under {root_dir}")
        return 0

    print(f"Discovered {len(grid_files)} grid_time files under {root_dir}")
    if workers > 1 and args.dry_run:
        print(f"DRY-RUN: requested workers={workers} (parallel execution disabled in dry-run)")
    elif workers > 1:
        print(f"Using {workers} concurrent workers for NetCDF postprocessing")

    ok = 0
    fail = 0
    tasks: list[dict] = []

    for grid_file in grid_files:
        m = RE_GRID_TIME.match(grid_file.name)
        if not m:
            continue
        ts14 = m.group(1)
        try:
            ts = datetime.strptime(ts14, "%Y%m%d%H%M%S")
        except ValueError:
            print(f"SKIP invalid timestamp in filename: {grid_file}")
            continue

        output_dir = grid_file.parent
        run_dir = output_dir.parent
        parsed = _parse_run_name(run_dir.name)
        if parsed is None:
            print(f"SKIP unrecognized run dir naming pattern: {run_dir.name}")
            continue

        domain, receptor, _ = parsed
        rel_h = release_heights.get(receptor, 10.0)
        magl = _format_magl_label(rel_h)
        out_name = f"{receptor}-{magl}magl_FLEXPART_CFSv2_{domain}_inert_{ts.strftime('%Y%m%d%H')}.nc"
        out_file = output_dir / out_name
        final_file = final_dir / out_name
        exit_csv = output_dir / out_name.replace(".nc", "_domain_exit_points.csv")

        if final_file.exists() and not args.overwrite:
            print(f"SKIP exists: {final_file}")
            continue

        task = {
            "grid_file": grid_file,
            "run_dir": run_dir,
            "out_file": out_file,
            "final_file": final_file,
            "needs_run": False,
            "cmd": None,
            "produced_ok": False,
            "returncode": 0,
        }

        if out_file.exists() and not args.overwrite:
            print(f"USE existing local postprocessed file: {out_file}")
            task["produced_ok"] = True
        else:
            cmd = [
                str(args.python),
                "-X",
                "faulthandler",
                str(args.postprocess_script),
                "--grid-file",
                str(grid_file),
                "--out-file",
                str(out_file),
                "--source-layer-thickness-m",
                str(args.postprocess_source_layer_thickness_m),
                "--footprint-outheight-m",
                str(footprint_outheight_m),
                "--partoutput",
                str(output_dir),
                "--exit-csv",
                str(exit_csv),
                "--site",
                receptor,
                "--domain",
                domain,
                "--species",
                "inert",
                "--model",
                "FLEXPART",
                "--met-model",
                "CFSv2",
            ]
            task["needs_run"] = True
            task["cmd"] = cmd

            print(f"RUN {grid_file} -> {out_file}")
            if args.dry_run:
                print(f"MOVE {out_file} -> {final_file}")
                if not args.keep_run_dirs:
                    print(f"DELETE run dir {run_dir}")

        tasks.append(task)

    if not args.dry_run:
        run_indices = [i for i, t in enumerate(tasks) if t["needs_run"] and t["cmd"] is not None]

        if workers == 1 or len(run_indices) <= 1:
            for i in run_indices:
                rc = _run_postprocess_cmd(tasks[i]["cmd"])
                tasks[i]["returncode"] = rc
                tasks[i]["produced_ok"] = (rc == 0)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                future_map = {
                    executor.submit(_run_postprocess_cmd, tasks[i]["cmd"]): i
                    for i in run_indices
                }
                for fut in as_completed(future_map):
                    i = future_map[fut]
                    rc = int(fut.result())
                    tasks[i]["returncode"] = rc
                    tasks[i]["produced_ok"] = (rc == 0)

    for t in tasks:
        grid_file = t["grid_file"]
        run_dir = t["run_dir"]
        out_file = t["out_file"]
        final_file = t["final_file"]

        if args.dry_run:
            continue

        if t["needs_run"] and not t["produced_ok"]:
            fail += 1
            print(f"FAIL ({t['returncode']}): {grid_file}")
            continue

        if final_file.exists() and args.overwrite:
            final_file.unlink()

        shutil.move(str(out_file), str(final_file))
        print(f"MOVED {out_file} -> {final_file}")

        if not args.keep_run_dirs:
            run_dir_resolved = run_dir.resolve()
            final_resolved = final_file.resolve()
            if final_resolved.is_relative_to(run_dir_resolved):
                print(
                    f"WARNING: final file is inside run directory; skipping deletion of {run_dir_resolved}"
                )
            else:
                shutil.rmtree(run_dir_resolved)
                print(f"DELETED run dir {run_dir_resolved}")

        ok += 1

    print(f"Summary: success={ok}, fail={fail}")

    if args.write_monthly:
        if not args.monthly_script.exists():
            print(f"ERROR: monthly aggregation script not found: {args.monthly_script}")
            fail += 1
        else:
            monthly_dir = args.monthly_dir if args.monthly_dir is not None else final_dir
            monthly_cmd = [
                str(args.python),
                str(args.monthly_script),
                "--input-dir",
                str(final_dir),
                "--output-dir",
                str(monthly_dir),
            ]
            if args.overwrite:
                monthly_cmd.append("--overwrite")

            if args.dry_run:
                print(f"MONTHLY {' '.join(monthly_cmd)}")
            else:
                print(f"RUN monthly aggregation in {monthly_dir}")
                env = os.environ.copy()
                env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
                monthly_result = subprocess.run(monthly_cmd, env=env)
                if monthly_result.returncode != 0:
                    fail += 1
                    print(f"FAIL monthly aggregation ({monthly_result.returncode})")

    return 1 if fail > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())