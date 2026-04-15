#!/usr/bin/env python3
"""Extract first particle domain-exit locations from FLEXPART partoutput NetCDF files.

Domain limits are derived from OUTGRID:
  lon in [OUTLON0, OUTLON0 + (NUMXGRID-1)*DXOUT]
  lat in [OUTLAT0, OUTLAT0 + (NUMYGRID-1)*DYOUT]
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import xarray as xr


@dataclass
class Domain:
    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float


def parse_outgrid(path: str) -> Domain:
    text = open(path, "r", encoding="utf-8").read()

    def get_float(key: str) -> float:
        m = re.search(rf"^\s*{key}\s*=\s*([-+0-9.eEdD]+)", text, re.MULTILINE)
        if not m:
            raise ValueError(f"Missing {key} in {path}")
        return float(m.group(1).replace("D", "E").replace("d", "e"))

    def get_int(key: str) -> int:
        return int(round(get_float(key)))

    outlon0 = get_float("OUTLON0")
    outlat0 = get_float("OUTLAT0")
    numx = get_int("NUMXGRID")
    numy = get_int("NUMYGRID")
    dx = get_float("DXOUT")
    dy = get_float("DYOUT")

    lon_max = outlon0 + (numx - 1) * dx
    lat_max = outlat0 + (numy - 1) * dy
    return Domain(lon_min=min(outlon0, lon_max), lon_max=max(outlon0, lon_max), lat_min=min(outlat0, lat_max), lat_max=max(outlat0, lat_max))


def find_var_name(ds: xr.Dataset, candidates: List[str]) -> Optional[str]:
    lower_map = {name.lower(): name for name in ds.variables}

    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]

    # fallback: variable contains candidate substring
    for cand in candidates:
        for name in ds.variables:
            if cand in name.lower():
                return name

    return None


def to_time_particle(da: xr.DataArray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dims = list(da.dims)
    if len(dims) < 2:
        raise ValueError(f"Variable {da.name} must have at least 2 dims, got {dims}")

    time_dim = next((d for d in dims if "time" in d.lower()), dims[0])
    part_dim = next((d for d in dims if "part" in d.lower()), None)
    if part_dim is None:
        part_dim = next((d for d in dims if d != time_dim), dims[1])

    indexers = {d: 0 for d in dims if d not in (time_dim, part_dim)}
    arr = da.isel(indexers).transpose(time_dim, part_dim)

    times = np.asarray(arr[time_dim].values)
    particles = np.asarray(arr[part_dim].values)
    values = np.asarray(arr.values)
    return values, times, particles


def format_time(val: object) -> str:
    if isinstance(val, np.datetime64):
        return str(np.datetime_as_string(val, unit="s"))
    return str(val)


def main() -> None:
    ap = argparse.ArgumentParser(description="Extract first exit location per particle from FLEXPART partoutput files.")
    ap.add_argument("--outgrid", required=True, help="Path to OUTGRID used for the run")
    ap.add_argument("--partglob", required=True, help="Glob for partoutput NetCDF files, e.g. output/partoutput*.nc")
    ap.add_argument("--output", required=True, help="Output CSV path")
    args = ap.parse_args()

    domain = parse_outgrid(args.outgrid)
    files = sorted(glob.glob(args.partglob))
    if not files:
        raise SystemExit(f"No partoutput files found for glob: {args.partglob}")

    exits: Dict[object, Dict[str, object]] = {}

    for nc in files:
        ds = xr.open_dataset(nc)
        try:
            lon_name = find_var_name(ds, ["longitude", "lon"])
            lat_name = find_var_name(ds, ["latitude", "lat"])
            hgt_name = find_var_name(ds, ["height", "z", "altitude"])

            if lon_name is None or lat_name is None:
                raise ValueError(f"Could not find lon/lat variables in {nc}")

            lon_vals, times, particles = to_time_particle(ds[lon_name])
            lat_vals, _, _ = to_time_particle(ds[lat_name])
            hgt_vals = None
            if hgt_name is not None:
                hgt_vals, _, _ = to_time_particle(ds[hgt_name])

            finite = np.isfinite(lon_vals) & np.isfinite(lat_vals)
            outside = finite & (
                (lon_vals < domain.lon_min)
                | (lon_vals > domain.lon_max)
                | (lat_vals < domain.lat_min)
                | (lat_vals > domain.lat_max)
            )

            nt, npart = outside.shape
            for j in range(npart):
                idx = np.where(outside[:, j])[0]
                if idx.size == 0:
                    continue
                i = int(idx[0])
                pid = particles[j].item() if hasattr(particles[j], "item") else particles[j]
                rec = {
                    "particle": pid,
                    "time": format_time(times[i]),
                    "lon": float(lon_vals[i, j]),
                    "lat": float(lat_vals[i, j]),
                    "height": float(hgt_vals[i, j]) if hgt_vals is not None and math.isfinite(float(hgt_vals[i, j])) else "",
                    "source_file": os.path.basename(nc),
                }

                # Keep earliest by lexical compare on time string as a practical default.
                old = exits.get(pid)
                if old is None or str(rec["time"]) < str(old["time"]):
                    exits[pid] = rec
        finally:
            ds.close()

    rows = sorted(exits.values(), key=lambda r: str(r["particle"]))
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["particle", "time", "lon", "lat", "height", "source_file"])
        w.writeheader()
        w.writerows(rows)

    print(
        f"Domain: lon [{domain.lon_min}, {domain.lon_max}], lat [{domain.lat_min}, {domain.lat_max}]\n"
        f"Scanned {len(files)} file(s), wrote {len(rows)} first-exit records to {args.output}."
    )


if __name__ == "__main__":
    main()
