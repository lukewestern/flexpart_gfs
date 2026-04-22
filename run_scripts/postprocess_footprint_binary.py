#!/usr/bin/env python3
"""
Postprocess FLEXPART binary grid_time files into a single footprint NetCDF.

This is intended for backward runs with SFC_ONLY=1, where FLEXPART writes
binary gridded output (grid_time_YYYYMMDDHHMMSS_NNN) instead of grid_time_*.nc.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import re
import struct
import sys
from pathlib import Path

import numpy as np
import xarray as xr


RE_GRID_TIME_BIN = re.compile(r"^grid_time_(\d{14})_(\d{3})$")
MOLAR_MASS_AIR_KG_PER_MOL = 0.02897


def _read_fortran_record(fh):
    """Read one little-endian unformatted Fortran record payload."""
    head = fh.read(4)
    if not head:
        return None
    if len(head) != 4:
        raise ValueError("truncated Fortran record header")

    n = struct.unpack("<i", head)[0]
    if n < 0:
        raise ValueError(f"invalid record size: {n}")

    payload = fh.read(n)
    if len(payload) != n:
        raise ValueError("truncated Fortran record payload")

    tail = fh.read(4)
    if len(tail) != 4:
        raise ValueError("truncated Fortran record trailer")

    n_tail = struct.unpack("<i", tail)[0]
    if n_tail != n:
        raise ValueError(f"record size mismatch: head={n}, tail={n_tail}")

    return payload


def _read_int_record(fh):
    payload = _read_fortran_record(fh)
    if payload is None:
        raise EOFError("unexpected EOF while reading int record")
    if len(payload) != 4:
        raise ValueError(f"expected int record of 4 bytes, got {len(payload)}")
    return struct.unpack("<i", payload)[0]


def _read_int_array_record(fh, n):
    payload = _read_fortran_record(fh)
    if payload is None:
        raise EOFError("unexpected EOF while reading int-array record")
    expected = 4 * int(n)
    if len(payload) != expected:
        raise ValueError(f"expected {expected} bytes for int-array record, got {len(payload)}")
    if n == 0:
        return np.empty((0,), dtype=np.int32)
    return np.frombuffer(payload, dtype="<i4")


def _read_float_array_record(fh, n):
    payload = _read_fortran_record(fh)
    if payload is None:
        raise EOFError("unexpected EOF while reading float-array record")
    expected = 4 * int(n)
    if len(payload) != expected:
        raise ValueError(f"expected {expected} bytes for float-array record, got {len(payload)}")
    if n == 0:
        return np.empty((0,), dtype=np.float32)
    return np.frombuffer(payload, dtype="<f4")


def _parse_header_txt(header_txt: Path):
    rows = []
    with open(header_txt) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line.split())

    if len(rows) < 4:
        raise ValueError(f"header_txt appears incomplete: {header_txt}")

    grid_row = rows[2]
    if len(grid_row) < 6:
        raise ValueError(f"could not parse grid setup row in {header_txt}")

    outlon0 = float(grid_row[0])
    outlat0 = float(grid_row[1])
    nx = int(float(grid_row[2]))
    ny = int(float(grid_row[3]))
    dx = float(grid_row[4])
    dy = float(grid_row[5])

    species_row = rows[5] if len(rows) > 5 else ["3", "1"]
    nspec3 = int(float(species_row[0]))
    nspec = max(1, nspec3 // 3)

    return {
        "outlon0": outlon0,
        "outlat0": outlat0,
        "nx": nx,
        "ny": ny,
        "dx": dx,
        "dy": dy,
        "nspec": nspec,
    }


def _parse_release_info(header_txt_releases: Path, nspec: int):
    if not header_txt_releases.exists():
        return np.nan, np.nan, np.nan

    rows = []
    with open(header_txt_releases) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line.split())

    if len(rows) < 3:
        return np.nan, np.nan, np.nan

    try:
        numpoint = int(float(rows[0][0]))
    except Exception:
        return np.nan, np.nan, np.nan

    if numpoint < 1:
        return np.nan, np.nan, np.nan

    # First release block: row1 timing/kind, row2 coords/heights, row3 npart,
    # row4 name, then 3*nspec mass rows.
    if len(rows) < 3:
        return np.nan, np.nan, np.nan

    try:
        coord_row = rows[2]
        xp1 = float(coord_row[0])
        yp1 = float(coord_row[1])
        xp2 = float(coord_row[2])
        yp2 = float(coord_row[3])
        z1 = float(coord_row[4])
        z2 = float(coord_row[5])
        rel_lon = 0.5 * (xp1 + xp2)
        rel_lat = 0.5 * (yp1 + yp2)
        rel_h = 0.5 * (z1 + z2)
        return rel_lon, rel_lat, rel_h
    except Exception:
        return np.nan, np.nan, np.nan


def _split_runs_by_sign(vals: np.ndarray):
    runs = []
    if vals.size == 0:
        return runs

    current = [0]
    current_sign = 1 if vals[0] >= 0 else -1

    for i in range(1, vals.size):
        s = 1 if vals[i] >= 0 else -1
        if s == current_sign:
            current.append(i)
        else:
            runs.append(np.array(current, dtype=np.int32))
            current = [i]
            current_sign = s
    runs.append(np.array(current, dtype=np.int32))
    return runs


def _decode_sparse_conc(nx: int, ny: int, run_starts: np.ndarray, run_vals: np.ndarray):
    # FLEXPART stores kz with 1-based offset in flattened index:
    # idx = ix + jy*nx + kz*nx*ny. For SFC_ONLY, kz=1 only.
    layer_offset = nx * ny
    flat_len = 2 * layer_offset

    field_flat = np.zeros(flat_len, dtype=np.float32)

    if run_starts.size == 0 or run_vals.size == 0:
        return field_flat[layer_offset : 2 * layer_offset].reshape((ny, nx))

    value_runs = _split_runs_by_sign(run_vals)

    if len(value_runs) < int(run_starts.size):
        raise ValueError(
            f"sparse decode mismatch: {run_starts.size} run starts but only {len(value_runs)} sign-runs"
        )

    for j in range(int(run_starts.size)):
        start = int(run_starts[j])
        vals = np.abs(run_vals[value_runs[j]])
        end = start + vals.size
        if start < 0:
            raise ValueError(f"negative sparse start index: {start}")
        if end > flat_len:
            # Expand defensively if needed for unusual indexing layouts.
            grow = end - flat_len
            field_flat = np.pad(field_flat, (0, grow), mode="constant")
            flat_len = field_flat.size
        field_flat[start:end] = vals

    return field_flat[layer_offset : 2 * layer_offset].reshape((ny, nx))


def _read_single_grid_time_file(path: Path, nx: int, ny: int):
    with open(path, "rb") as fh:
        _itime = _read_int_record(fh)

        # Block 1: wet deposition sparse payload.
        n_i_wet = _read_int_record(fh)
        _ = _read_int_array_record(fh, n_i_wet)
        n_r_wet = _read_int_record(fh)
        _ = _read_float_array_record(fh, n_r_wet)

        # Block 2: dry deposition sparse payload.
        n_i_dry = _read_int_record(fh)
        _ = _read_int_array_record(fh, n_i_dry)
        n_r_dry = _read_int_record(fh)
        _ = _read_float_array_record(fh, n_r_dry)

        # Block 3: concentration/time sparse payload.
        n_i = _read_int_record(fh)
        run_starts = _read_int_array_record(fh, n_i)
        n_r = _read_int_record(fh)
        run_vals = _read_float_array_record(fh, n_r)

    return _decode_sparse_conc(nx, ny, run_starts, run_vals)


def _set_netcdf_compression(ds, compression_level=4):
    encoding = {}
    for var in ds.data_vars:
        encoding[var] = {"zlib": True, "complevel": int(compression_level), "shuffle": True}
    return encoding


def _convert_to_m2s_per_mol(da_2d: np.ndarray, source_layer_thickness_m: float):
    """Convert backward SRR from native s m3 kg-1 to m2 s mol-1."""
    h = float(source_layer_thickness_m)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("source layer thickness must be a positive finite number (m)")

    factor = MOLAR_MASS_AIR_KG_PER_MOL / h
    return da_2d * factor, factor


def _infer_reference_time(output_dir: Path) -> dt.datetime:
    """Infer release/reference datetime from run directory name.

    Expected run directory naming: <domain>_<receptor>_<YYYYMMDDHH>/output
    """
    run_name = output_dir.parent.name
    parts = run_name.rsplit("_", 2)
    if len(parts) == 3 and re.match(r"^\d{10}$", parts[2]):
        try:
            return dt.datetime.strptime(parts[2], "%Y%m%d%H")
        except ValueError:
            pass
    # Conservative fallback if run naming is unexpected.
    return dt.datetime(1970, 1, 1, 0, 0, 0)


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert FLEXPART binary grid_time files to one footprint NetCDF")
    ap.add_argument("--output-dir", required=True, type=Path, help="FLEXPART output directory containing binary grid_time files")
    ap.add_argument("--out-file", required=True, type=Path, help="Output NetCDF footprint path")
    ap.add_argument("--site", default="UNKNOWN")
    ap.add_argument("--domain", default="UNKNOWN")
    ap.add_argument("--species", default="inert")
    ap.add_argument("--model", default="FLEXPART")
    ap.add_argument("--met-model", default="GFS")
    ap.add_argument(
        "--source-layer-thickness-m",
        type=float,
        default=100.0,
        help="Accepted for interface compatibility; binary path writes native seconds units.",
    )
    args = ap.parse_args()

    output_dir = args.output_dir
    if not output_dir.is_dir():
        raise FileNotFoundError(f"output directory not found: {output_dir}")

    header_txt = output_dir / "header_txt"
    if not header_txt.exists():
        raise FileNotFoundError(f"header_txt not found in {output_dir}")

    info = _parse_header_txt(header_txt)
    nx = info["nx"]
    ny = info["ny"]
    outlon0 = info["outlon0"]
    outlat0 = info["outlat0"]
    dx = info["dx"]
    dy = info["dy"]

    binary_files = sorted(
        p for p in output_dir.iterdir() if p.is_file() and RE_GRID_TIME_BIN.match(p.name)
    )
    if not binary_files:
        raise FileNotFoundError(f"no binary grid_time_YYYYMMDDHHMMSS_NNN files found in {output_dir}")

    summed = np.zeros((ny, nx), dtype=np.float64)
    for p in binary_files:
        summed += _read_single_grid_time_file(p, nx=nx, ny=ny)

    lons = outlon0 + np.arange(nx, dtype=np.float64) * dx
    lats = outlat0 + np.arange(ny, dtype=np.float64) * dy
    ref_time = _infer_reference_time(output_dir)
    time_units = f"seconds since {ref_time.strftime('%Y-%m-%d %H:%M:%S')}"

    rel_lon, rel_lat, rel_h = _parse_release_info(output_dir / "header_txt_releases", nspec=info["nspec"])

    summed_conv, conv_factor = _convert_to_m2s_per_mol(
        summed.astype(np.float32),
        args.source_layer_thickness_m,
    )

    out = xr.Dataset()
    out["time"] = xr.DataArray(np.array([0.0], dtype=np.float64), dims=("time",))
    out["longitude"] = xr.DataArray(lons, dims=("longitude",), attrs={"units": "degrees_east", "long_name": "longitude"})
    out["latitude"] = xr.DataArray(lats, dims=("latitude",), attrs={"units": "degrees_north", "long_name": "latitude"})

    out["srr"] = xr.DataArray(
        summed_conv[None, :, :],
        dims=("time", "latitude", "longitude"),
        coords={"time": out["time"], "latitude": out["latitude"], "longitude": out["longitude"]},
        attrs={
            "long_name": "source_receptor_relationship",
            "units": "m2 s mol-1",
            "source_format": "FLEXPART binary grid_time",
            "description": "sum over all binary grid_time slices in run output directory",
            "loss_lifetime_hrs": -9.0,
            "loss_lifetime_comment": "lifetime in hours; -9 corresponds to inert",
            "conversion_from_native_units": "s m3 kg-1",
            "conversion_factor_applied": float(conv_factor),
            "source_layer_thickness_m": float(args.source_layer_thickness_m),
            "molar_mass_air_kg_per_mol": float(MOLAR_MASS_AIR_KG_PER_MOL),
        },
    )

    ntime = int(out.sizes.get("time", 1))
    out["release_lon"] = xr.DataArray(np.full(ntime, rel_lon, dtype=np.float32), dims=("time",), coords={"time": out["time"]}, attrs={"units": "degrees_east", "long_name": "Release longitude"})
    out["release_lat"] = xr.DataArray(np.full(ntime, rel_lat, dtype=np.float32), dims=("time",), coords={"time": out["time"]}, attrs={"units": "degrees_north", "long_name": "Release latitude"})
    out["release_height"] = xr.DataArray(np.full(ntime, rel_h, dtype=np.float32), dims=("time",), coords={"time": out["time"]}, attrs={"units": "m", "long_name": "Release height above model ground"})

    out["time"].attrs.update(
        {
            "long_name": "time",
            "standard_name": "time",
            "label": "left",
            "comment": "single release-time-integrated footprint",
            "units": time_units,
            "calendar": "gregorian",
            "period": "release-integrated",
        }
    )

    out.attrs.update(
        {
            "title": "Derived FLEXPART footprint products from binary grid_time",
            "input_output_dir": str(output_dir.resolve()),
            "species": args.species,
            "model": args.model,
            "met_model": args.met_model,
            "domain": args.domain,
            "site": args.site,
            "lpdm_native_output_unit": "s m3 kg-1",
            "note": "Binary fallback path for SFC_ONLY=1 where NetCDF gridded output is unavailable.",
            "created": dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z"),
            "author": "postprocess_footprint_binary.py",
        }
    )

    out_file = args.out_file
    out_file.parent.mkdir(parents=True, exist_ok=True)
    encoding = _set_netcdf_compression(out, compression_level=4)
    out.to_netcdf(out_file, encoding=encoding)

    print(f"Read {len(binary_files)} binary grid_time files from {output_dir}")
    print(f"Wrote binary-derived footprint NetCDF: {out_file}")
    print(
        "Applied unit conversion: native s m3 kg-1 -> m2 s mol-1 "
        f"with factor {conv_factor:.8g} (source layer {args.source_layer_thickness_m:g} m)"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
