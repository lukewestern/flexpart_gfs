#!/usr/bin/env python3
"""
Postprocess daily FLEXPART grid_time files into AGAGE-style footprints.

Daily FLEXPART runs can write one grid_time_*.nc file containing 24 hourly
release groups. In those files the gridded sensitivity variable has a
``pointspec`` dimension with one slice per release. This script preserves that
dimension and writes it as the output ``time`` dimension, giving one daily
NetCDF file with 24 hourly footprints.
"""

from __future__ import annotations

import argparse
import calendar
import csv
import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_ROOT / "run_scripts"))

from postprocess_footprint import (  # noqa: E402
    DEFAULT_CFSV2_GRIB_DIR,
    HEIGHT_BINS_M_AGL,
    MOLAR_MASS_AIR_KG_PER_MOL,
    _build_boundary_exit_fractions,
    _build_time_attrs,
    _classify_exit_side,
    _convert_to_m2s_per_mol,
    _derive_exit_points,
    _extract_release_meteo_from_grib,
    _find_spatial_dims,
    _find_gf_file_for_time,
    _find_particle_vars,
    _open_dataset_auto,
    _open_partoutput,
    _pick_sensitivity_var,
    _set_netcdf_compression,
)


RE_GRID_TIME = re.compile(r"^grid_time_(\d{14})\.nc$")
RE_RUN_DIR = re.compile(r"^(.+)_([^_]+)_(\d{8})$")

HELP_EPILOG = r"""
Examples:
  # Postprocess one daily run directory in place.
  python run_scripts_daily/postprocess_daily_footprints.py \
      --run-dir /net/fs06/d2/lwestern/flexpart_outs_daily/WESTUSA_THD_20180101

  # Postprocess one explicit grid file and choose the output file name.
  python run_scripts_daily/postprocess_daily_footprints.py \
      --grid-file /net/fs06/d2/lwestern/flexpart_outs_daily/WESTUSA_THD_20180101/output/grid_time_20180101230000.nc \
      --out-file /net/fs06/d2/lwestern/flexpart_outs_daily/WESTUSA_THD_20180101/output/THD-10magl_FLEXPART_CFSv2_WESTUSA_inert_20180101.nc

  # Scan a root directory and postprocess all daily outputs found below it.
  python run_scripts_daily/postprocess_daily_footprints.py \
      --root-dir /net/fs06/d2/lwestern/flexpart_outs_daily

  # Preview work without writing files.
  python run_scripts_daily/postprocess_daily_footprints.py \
      --root-dir /net/fs06/d2/lwestern/flexpart_outs_daily \
      --dry-run --limit 5

Notes:
  The default output is written beside each grid_time_*.nc file using a name like
  SITE-HEIGHTmagl_FLEXPART_CFSv2_DOMAIN_inert_YYYYMMDD.nc. Use --output-dir to
  redirect generated files, or --out-file for a single --grid-file run.
"""


def _parse_run_dir(path: Path) -> tuple[str, str, dt.date] | None:
    """Parse DOMAIN_SITE_YYYYMMDD daily run directory names."""
    match = RE_RUN_DIR.match(path.name)
    if match is None:
        return None
    domain, site, yyyymmdd = match.groups()
    try:
        return domain, site, dt.datetime.strptime(yyyymmdd, "%Y%m%d").date()
    except ValueError:
        return None


def _parse_grid_datetime(path: Path) -> dt.datetime | None:
    match = RE_GRID_TIME.match(path.name)
    if match is None:
        return None
    try:
        return dt.datetime.strptime(match.group(1), "%Y%m%d%H%M%S")
    except ValueError:
        return None


def _release_datetimes_for_grid(grid_file: Path, run_date: dt.date | None, ntime: int) -> list[dt.datetime]:
    grid_dt = _parse_grid_datetime(grid_file)
    if ntime <= 1:
        if grid_dt is not None:
            return [grid_dt]
        if run_date is not None:
            return [dt.datetime.combine(run_date, dt.time())]
        return []

    if run_date is not None:
        base = dt.datetime.combine(run_date, dt.time())
    elif grid_dt is not None:
        base = dt.datetime.combine(grid_dt.date(), dt.time())
    else:
        return []
    return [base + dt.timedelta(hours=hour) for hour in range(ntime)]


def _extract_daily_release_meteo(
    grid_file: Path,
    ds: xr.Dataset,
    run_date: dt.date | None,
    ntime: int,
    meteo_grib_dir: Path | None,
    meteo_grib_file: Path | None,
) -> dict[str, np.ndarray]:
    release_times = _release_datetimes_for_grid(grid_file, run_date, ntime)
    if len(release_times) != ntime:
        print("WARNING: could not infer daily release datetimes; keeping release meteorology as NaN")
        return {}

    rel_lons = _release_metadata_1d(ds, "RELLNG1", ntime, dtype=np.float64)
    rel_lats = _release_metadata_1d(ds, "RELLAT1", ntime, dtype=np.float64)
    fields = [
        "air_temperature",
        "air_pressure",
        "wind_speed",
        "wind_from_direction",
        "atmosphere_boundary_layer_thickness",
    ]
    release_meteo = {
        name: np.full(ntime, np.nan, dtype=np.float32)
        for name in fields
    }
    extracted = 0
    last_grib_file: Path | str | None = None

    for i, release_dt in enumerate(release_times):
        rel_lon = float(rel_lons[i]) if i < rel_lons.size else np.nan
        rel_lat = float(rel_lats[i]) if i < rel_lats.size else np.nan
        grib_file = meteo_grib_file
        if grib_file is None:
            grib_file = _find_gf_file_for_time(
                str(meteo_grib_dir) if meteo_grib_dir is not None else None,
                release_dt,
            )
        if grib_file is None or not Path(grib_file).is_file():
            continue

        values = _extract_release_meteo_from_grib(str(grib_file), release_dt, rel_lon, rel_lat)
        if values:
            extracted += 1
            last_grib_file = grib_file
        for name, value in values.items():
            if name in release_meteo and np.isfinite(value):
                release_meteo[name][i] = np.float32(value)

    if extracted:
        print(f"Extracted release meteorology for {extracted}/{ntime} release time(s)")
        if last_grib_file is not None:
            print(f"Last GRIB meteorology source: {last_grib_file}")
    else:
        print("WARNING: could not extract release meteorology from GRIB; keeping NaNs")
    return release_meteo


def _format_magl_label(height: float) -> str:
    height = float(height)
    if height.is_integer():
        return str(int(height))
    return f"{height:.3f}".rstrip("0").rstrip(".")


def _as_float_seconds(values) -> np.ndarray:
    arr = np.asarray(values)
    if np.issubdtype(arr.dtype, np.timedelta64):
        return arr.astype("timedelta64[s]").astype(np.float64)
    if np.issubdtype(arr.dtype, np.datetime64):
        return arr.astype("datetime64[s]").astype(np.float64)
    return arr.astype(np.float64)


def _release_time_values(ds: xr.Dataset, release_dim: str | None, nrelease: int) -> np.ndarray:
    """Return one native numeric time value per release."""
    for name in ["RELSTART", "RELEND"]:
        if name not in ds:
            continue
        da = ds[name]
        if release_dim in da.dims or "numpoint" in da.dims or da.size == nrelease:
            return _as_float_seconds(da.values).ravel()[:nrelease]
    return np.arange(nrelease, dtype=np.float64)


def _coord_1d(ds: xr.Dataset, dim: str, out_name: str) -> xr.DataArray:
    if dim in ds.coords:
        return ds[dim].rename(out_name)
    if dim in ds:
        return ds[dim].rename(out_name)
    return xr.DataArray(np.arange(ds.sizes[dim], dtype=np.float32), dims=(out_name,))


def _release_metadata_1d(
    ds: xr.Dataset,
    name: str,
    nrelease: int,
    dtype,
    default=np.nan,
) -> np.ndarray:
    if name not in ds:
        return np.full(nrelease, default, dtype=dtype)

    da = ds[name]
    if "numpoint" in da.dims:
        values = np.asarray(da.values)[:nrelease]
    else:
        values = np.asarray(da.values).ravel()
        if values.size != nrelease:
            values = np.resize(values, nrelease)
    return values.astype(dtype)


def _infer_release_height(ds: xr.Dataset) -> float:
    for name in ["RELZZ1", "RELZZ2"]:
        if name in ds:
            values = np.asarray(ds[name].values, dtype=float).ravel()
            finite = values[np.isfinite(values)]
            if finite.size:
                return float(finite[0])
    return float("nan")


def _input_release_count(ds: xr.Dataset) -> int:
    try:
        var_name = _pick_sensitivity_var(ds)
    except ValueError:
        return 1

    da = ds[var_name]
    for dim in ["pointspec", "numpoint"]:
        if dim in da.dims:
            return int(da.sizes[dim])
    return 1


def compute_daily_srr(
    ds: xr.Dataset,
    footprint_outheight_m: float,
    source_layer_thickness_m: float,
) -> tuple[xr.DataArray, str, str | None, float | None, int | None, float]:
    """Integrate sensitivity while preserving one slice per release."""
    var_name = _pick_sensitivity_var(ds)
    da = ds[var_name]

    if "pointspec" in da.dims:
        release_dim = "pointspec"
    elif "numpoint" in da.dims:
        release_dim = "numpoint"
    else:
        release_dim = None

    selected_height = None
    integrated_height_levels = None
    if "height" in da.dims:
        if "height" not in ds:
            raise ValueError("Sensitivity has a height dimension but no height coordinate")
        hvals = np.asarray(ds["height"].values, dtype=float)
        if hvals.size == 0:
            raise ValueError("Height coordinate is empty")
        hidx = int(np.argmin(np.abs(hvals - float(footprint_outheight_m))))
        selected_height = float(hvals[hidx])
        integrated_height_levels = int(hidx + 1)
        da = da.isel(height=slice(0, hidx + 1))
        if not np.isclose(selected_height, float(footprint_outheight_m), atol=1e-6):
            print(
                "WARNING: requested footprint outheight "
                f"{float(footprint_outheight_m):.3f} m not found; "
                f"using nearest {selected_height:.3f} m"
            )

    reduce_dims = ["time", "height", "nageclass", "numspec"]
    active_reduce_dims = [dim for dim in reduce_dims if dim in da.dims]
    if active_reduce_dims:
        da = da.sum(dim=active_reduce_dims, skipna=True)

    native_units = ds[var_name].attrs.get("units", "")
    da, conversion_factor = _convert_to_m2s_per_mol(
        da,
        native_units,
        source_layer_thickness_m,
    )

    if release_dim is None:
        nrelease = 1
        da = da.expand_dims(time=np.array([0.0], dtype=np.float64))
    else:
        nrelease = int(da.sizes[release_dim])
        da = da.rename({release_dim: "time"})
        da = da.assign_coords(time=_release_time_values(ds, release_dim, nrelease))

    lat_dim, lon_dim = _find_spatial_dims(ds)
    rename_dims = {}
    if lat_dim != "latitude":
        rename_dims[lat_dim] = "latitude"
    if lon_dim != "longitude":
        rename_dims[lon_dim] = "longitude"
    if rename_dims:
        da = da.rename(rename_dims)

    da = da.transpose("time", "latitude", "longitude").astype(np.float32)
    da.attrs.update(
        {
            "long_name": "source_receptor_relationship",
            "loss_lifetime_hrs": -9.0,
            "loss_lifetime_comment": "lifetime in hours; -9 corresponds to inert",
            "units": "m2 s mol-1",
            "source_variable": var_name,
            "description": "time-integrated hourly-release footprint from daily FLEXPART grid_time output",
            "conversion_from_native_units": str(native_units),
            "conversion_factor_applied": float(conversion_factor),
            "source_layer_thickness_m": float(source_layer_thickness_m),
            "molar_mass_air_kg_per_mol": float(MOLAR_MASS_AIR_KG_PER_MOL),
        }
    )
    if selected_height is not None:
        da.attrs.update(
            {
                "footprint_outheight_m": float(selected_height),
                "requested_footprint_outheight_m": float(footprint_outheight_m),
                "footprint_height_integration": "surface_to_outheight_inclusive",
                "footprint_height_levels_integrated": int(integrated_height_levels),
            }
        )

    return da, var_name, release_dim, selected_height, integrated_height_levels, conversion_factor


def build_daily_footprint_dataset(
    grid_file: Path,
    domain: str,
    site: str,
    run_date: dt.date | None,
    species: str,
    model: str,
    met_model: str,
    footprint_outheight_m: float,
    source_layer_thickness_m: float,
    meteo_grib_dir: Path | None = None,
    meteo_grib_file: Path | None = None,
    disable_grib_meteo: bool = False,
) -> xr.Dataset:
    ds = _open_dataset_auto(str(grid_file))
    try:
        lat_dim, lon_dim = _find_spatial_dims(ds)
        srr, var_name, release_dim, _, _, _ = compute_daily_srr(
            ds,
            footprint_outheight_m=footprint_outheight_m,
            source_layer_thickness_m=source_layer_thickness_m,
        )
        ntime = int(srr.sizes["time"])

        out = xr.Dataset(
            coords={
                "time": srr["time"],
                "latitude": _coord_1d(ds, lat_dim, "latitude"),
                "longitude": _coord_1d(ds, lon_dim, "longitude"),
            }
        )
        out["srr"] = srr

        out["time"].attrs.update(_build_time_attrs(ds))
        out["time"].attrs["comment"] = (
            "hourly release-time stamps from RELSTART/RELEND in daily FLEXPART output"
        )
        out["latitude"].attrs.update({"units": "degrees_north", "long_name": "latitude"})
        out["longitude"].attrs.update({"units": "degrees_east", "long_name": "longitude"})

        release_meteo: dict[str, np.ndarray] = {}
        if not disable_grib_meteo:
            release_meteo = _extract_daily_release_meteo(
                grid_file=grid_file,
                ds=ds,
                run_date=run_date,
                ntime=ntime,
                meteo_grib_dir=meteo_grib_dir,
                meteo_grib_file=meteo_grib_file,
            )

        release_vars = {
            "release_lon": ("RELLNG1", np.float32, "degrees_east", "Release longitude"),
            "release_lat": ("RELLAT1", np.float32, "degrees_north", "Release latitude"),
            "release_height": ("RELZZ1", np.float32, "m", "Release height above model ground"),
            "release_height_top": ("RELZZ2", np.float32, "m", "Release height top above model ground"),
            "release_particles": ("RELPART", np.int32, "1", "Number of release particles"),
        }
        for out_name, (src_name, dtype, units, long_name) in release_vars.items():
            out[out_name] = xr.DataArray(
                _release_metadata_1d(ds, src_name, ntime, dtype=dtype),
                dims=("time",),
                coords={"time": out["time"]},
                attrs={"units": units, "long_name": long_name},
            )

        for name, units, long_name in [
            ("air_temperature", "K", "air temperature at release"),
            ("air_pressure", "hPa", "air pressure at release"),
            ("wind_speed", "m s-1", "wind speed at release"),
            ("wind_from_direction", "degree", "wind direction at release"),
            ("atmosphere_boundary_layer_thickness", "m", "atmospheric boundary layer thickness at release"),
        ]:
            values = release_meteo.get(name)
            if values is None:
                values = np.full(ntime, np.nan, dtype=np.float32)
            else:
                values = np.asarray(values, dtype=np.float32).ravel()
                if values.size != ntime:
                    values = np.resize(values, ntime).astype(np.float32)
            out[name] = xr.DataArray(
                values,
                dims=("time",),
                coords={"time": out["time"]},
                attrs={"units": units, "long_name": long_name},
            )

        if release_meteo:
            out.attrs["release_meteo_source"] = str(meteo_grib_file or meteo_grib_dir or "unknown")
            out.attrs["release_meteo_extraction"] = "nearest GRIB grid point at each hourly release"

        out.attrs.update(
            {
                "title": "Derived FLEXPART daily hourly-release footprint products",
                "input_grid_file": str(grid_file.resolve()),
                "note": "Backward grid_time is source-receptor sensitivity; pointspec is written as hourly release time.",
                "lpdm_native_output_unit": "s",
                "species": species,
                "model": model,
                "met_model": met_model,
                "output_folder": "output",
                "model_version": str(ds.attrs.get("source", "FLEXPART")),
                "domain": domain,
                "site": site,
                "author": "postprocess_daily_footprints.py",
                "created": dt.datetime.utcnow().isoformat() + "Z",
                "source_sensitivity_variable": var_name,
                "release_dimension_in_input": str(release_dim or "none"),
            }
        )
        return out
    finally:
        ds.close()


def _daily_release_particle_slices(
    out: xr.Dataset,
    particle_offset: int = 0,
) -> list[tuple[int, int, int]]:
    """Return (release_index, start, stop) particle slices for daily releases."""
    if "release_particles" not in out:
        return []

    counts = np.asarray(out["release_particles"].values, dtype=float).ravel()
    if counts.size != int(out.sizes.get("time", counts.size)):
        return []
    if not np.all(np.isfinite(counts)):
        return []

    counts = counts.astype(np.int64)
    if np.any(counts < 0):
        return []

    starts = np.concatenate([[0], np.cumsum(counts[:-1])]) + int(particle_offset)
    stops = starts + counts
    return [(int(i), int(s), int(e)) for i, (s, e) in enumerate(zip(starts, stops))]


def _daily_release_particle_count(out: xr.Dataset) -> int:
    if "release_particles" not in out:
        return 0
    counts = np.asarray(out["release_particles"].values, dtype=float).ravel()
    finite = counts[np.isfinite(counts)]
    if finite.size != counts.size:
        return 0
    return int(np.sum(finite.astype(np.int64)))


def _subset_exits_by_particle_range(
    exits: list[tuple[int, float, float, float, float]],
    start: int,
    stop: int,
) -> list[tuple[int, float, float, float, float]]:
    return [exit_rec for exit_rec in exits if start <= exit_rec[0] < stop]


def _prepare_exit_lookup(
    exits: list[tuple[int, float, float, float, float]],
) -> tuple[np.ndarray, list[tuple[int, float, float, float, float]]]:
    if not exits:
        return np.array([], dtype=np.int64), []

    order = np.argsort(np.fromiter((rec[0] for rec in exits), dtype=np.int64, count=len(exits)))
    sorted_exits = [exits[int(i)] for i in order]
    particle_ids = np.fromiter((rec[0] for rec in sorted_exits), dtype=np.int64, count=len(sorted_exits))
    return particle_ids, sorted_exits


def _subset_exit_lookup_by_particle_range(
    particle_ids: np.ndarray,
    sorted_exits: list[tuple[int, float, float, float, float]],
    start: int,
    stop: int,
) -> list[tuple[int, float, float, float, float]]:
    if particle_ids.size == 0:
        return []
    i0 = int(np.searchsorted(particle_ids, int(start), side="left"))
    i1 = int(np.searchsorted(particle_ids, int(stop), side="left"))
    return sorted_exits[i0:i1]


def _write_daily_exit_csv(
    exit_csv: Path,
    exits_by_release: list[list[tuple[int, float, float, float, float]]],
    release_times: np.ndarray,
    lon_centers: np.ndarray,
    lat_centers: np.ndarray,
) -> None:
    exit_csv.parent.mkdir(parents=True, exist_ok=True)
    lon_min = float(np.min(lon_centers))
    lon_max = float(np.max(lon_centers))
    lat_min = float(np.min(lat_centers))
    lat_max = float(np.max(lat_centers))

    with exit_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "release_index",
            "release_time",
            "particle_index",
            "exit_time",
            "exit_longitude",
            "exit_latitude",
            "exit_height_magl",
            "exit_side",
        ])
        for release_index, exits in enumerate(exits_by_release):
            release_time = release_times[release_index] if release_index < len(release_times) else np.nan
            for pidx, tval, lon, lat, z in exits:
                side = _classify_exit_side(lon, lat, lon_min, lon_max, lat_min, lat_max)
                writer.writerow([release_index, release_time, pidx, tval, lon, lat, z, side])


def load_daily_exit_records(partoutput: str | Path) -> tuple[list[tuple[int, float, float, float, float]], int] | None:
    """Load partoutput once and derive all daily particle exits."""
    pds = _open_partoutput(str(partoutput))
    if pds is None:
        print(f"No partoutput files supplied/found at {partoutput}; skipping domain-exit diagnostics.")
        return None

    try:
        lon_var, lat_var, z_var = _find_particle_vars(pds)
        if lon_var is None or lat_var is None:
            print("Could not find particle longitude/latitude vars; skipping domain exits.")
            return None
        return _derive_exit_points(pds, lon_var, lat_var, z_var=z_var)
    finally:
        pds.close()


def add_daily_exit_diagnostics(
    out: xr.Dataset,
    partoutput: str | Path,
    exit_csv: Path | None = None,
    exit_records: tuple[list[tuple[int, float, float, float, float]], int] | None = None,
    particle_offset: int = 0,
) -> bool:
    """Add boundary exit-fraction variables after deriving exits once per day."""
    if exit_records is None:
        exit_records = load_daily_exit_records(partoutput)
    if exit_records is None:
        return False

    exits, npart = exit_records
    particle_ids, sorted_exits = _prepare_exit_lookup(exits)

    release_slices = _daily_release_particle_slices(out, particle_offset=particle_offset)
    if not release_slices:
        print("Could not infer release particle ranges; skipping daily domain-exit diagnostics.")
        return False

    expected_particles = release_slices[-1][2] if release_slices else 0
    if expected_particles > npart:
        print(
            "WARNING: release particle ranges expect "
            f"{expected_particles} particles, but partoutput contains {npart}; clipping to available particles."
        )

    lon_centers = np.asarray(out["longitude"].values, dtype=float)
    lat_centers = np.asarray(out["latitude"].values, dtype=float)
    release_times = np.asarray(out["time"].values, dtype=float)
    nt = int(out.sizes["time"])
    nh = len(HEIGHT_BINS_M_AGL)

    frac_n = np.zeros((nt, nh, len(lon_centers)), dtype=np.float32)
    frac_s = np.zeros_like(frac_n)
    frac_e = np.zeros((nt, nh, len(lat_centers)), dtype=np.float32)
    frac_w = np.zeros_like(frac_e)

    exits_by_release: list[list[tuple[int, float, float, float, float]]] = []
    for release_index, start, stop in release_slices:
        release_exits = _subset_exit_lookup_by_particle_range(
            particle_ids,
            sorted_exits,
            start,
            min(stop, npart),
        )
        exits_by_release.append(release_exits)
        n_slice, e_slice, s_slice, w_slice = _build_boundary_exit_fractions(
            release_exits,
            time_vals=np.asarray([release_times[release_index]], dtype=float),
            lon_centers=lon_centers,
            lat_centers=lat_centers,
            height_centers=HEIGHT_BINS_M_AGL,
        )
        frac_n[release_index] = n_slice[0]
        frac_e[release_index] = e_slice[0]
        frac_s[release_index] = s_slice[0]
        frac_w[release_index] = w_slice[0]

    out["height"] = xr.DataArray(
        HEIGHT_BINS_M_AGL.astype(np.float32),
        dims=("height",),
        attrs={"long_name": "height at layer midpoints", "units": "m", "positive": "up"},
    )
    out["particle_locations_n"] = xr.DataArray(
        frac_n,
        dims=("time", "height", "longitude"),
        coords={"time": out["time"], "height": out["height"], "longitude": out["longitude"]},
        attrs={"long_name": "Fraction of exiting particles leaving domain (N side)", "units": "1"},
    )
    out["particle_locations_e"] = xr.DataArray(
        frac_e,
        dims=("time", "height", "latitude"),
        coords={"time": out["time"], "height": out["height"], "latitude": out["latitude"]},
        attrs={"long_name": "Fraction of exiting particles leaving domain (E side)", "units": "1"},
    )
    out["particle_locations_s"] = xr.DataArray(
        frac_s,
        dims=("time", "height", "longitude"),
        coords={"time": out["time"], "height": out["height"], "longitude": out["longitude"]},
        attrs={"long_name": "Fraction of exiting particles leaving domain (S side)", "units": "1"},
    )
    out["particle_locations_w"] = xr.DataArray(
        frac_w,
        dims=("time", "height", "latitude"),
        coords={"time": out["time"], "height": out["height"], "latitude": out["latitude"]},
        attrs={"long_name": "Fraction of exiting particles leaving domain (W side)", "units": "1"},
    )

    out.attrs["domain_exit_particle_count"] = int(len(exits))
    out.attrs["domain_exit_partoutput"] = str(partoutput)
    out.attrs["domain_exit_attribution"] = "particle index ranges from daily RELEASES order"
    out.attrs["domain_exit_particle_offset"] = int(particle_offset)

    if exit_csv is not None:
        _write_daily_exit_csv(exit_csv, exits_by_release, release_times, lon_centers, lat_centers)
        print(f"Wrote daily domain exit points CSV: {exit_csv}")

    print(f"Derived {len(exits)} particle exit points once for {len(release_slices)} release group(s)")
    return True


def _infer_context(grid_file: Path, domain: str | None, site: str | None) -> tuple[str, str, dt.date | None]:
    run_dir = grid_file.parent.parent if grid_file.parent.name == "output" else grid_file.parent
    parsed = _parse_run_dir(run_dir)
    run_date = None
    if parsed is not None:
        inferred_domain, inferred_site, run_date = parsed
        domain = domain or inferred_domain
        site = site or inferred_site

    grid_dt = _parse_grid_datetime(grid_file)
    if run_date is None and grid_dt is not None:
        run_date = grid_dt.date()

    if domain is None or site is None:
        raise ValueError(
            "Could not infer --domain/--site from the grid file path; pass both explicitly."
        )
    return domain, site, run_date


def _default_output_file(
    grid_file: Path,
    domain: str,
    site: str,
    species: str,
    output_dir: Path | None,
) -> Path:
    ds = _open_dataset_auto(str(grid_file))
    try:
        release_height = _infer_release_height(ds)
        release_count = _input_release_count(ds)
    finally:
        ds.close()

    height_label = "unknown" if not np.isfinite(release_height) else _format_magl_label(release_height)
    grid_dt = _parse_grid_datetime(grid_file)
    if grid_dt is None:
        date_label = grid_file.stem
    elif release_count > 1:
        date_label = grid_dt.strftime("%Y%m%d")
    else:
        date_label = grid_dt.strftime("%Y%m%d%H")
    out_name = f"{site}-{height_label}magl_FLEXPART_CFSv2_{domain}_{species}_{date_label}.nc"
    return (output_dir or grid_file.parent) / out_name


def process_grid_file(
    grid_file: Path,
    args: argparse.Namespace,
    output_dir: Path | None = None,
    exit_records_cache: dict[Path, tuple[list[tuple[int, float, float, float, float]], int] | None] | None = None,
    particle_offset: int = 0,
) -> Path | None:
    domain, site, run_date = _infer_context(grid_file, args.domain, args.site)
    out_file = args.out_file
    if out_file is None:
        out_file = _default_output_file(grid_file, domain, site, args.species, output_dir)
    else:
        out_file = Path(out_file)

    if out_file.exists() and not args.overwrite:
        print(f"SKIP existing: {out_file}")
        return out_file

    print(f"POSTPROCESS {grid_file} -> {out_file}")
    if args.dry_run:
        return out_file

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out = build_daily_footprint_dataset(
        grid_file=grid_file,
        domain=domain,
        site=site,
        run_date=run_date,
        species=args.species,
        model=args.model,
        met_model=args.met_model,
        footprint_outheight_m=args.footprint_outheight_m,
        source_layer_thickness_m=args.source_layer_thickness_m,
        meteo_grib_dir=args.meteo_grib_dir,
        meteo_grib_file=args.meteo_grib_file,
        disable_grib_meteo=args.disable_grib_meteo,
    )
    try:
        if args.partoutput is not None:
            partoutput = args.partoutput
        elif args.auto_partoutput:
            partoutput = grid_file.parent
        else:
            partoutput = None

        if partoutput is not None:
            partoutput_path = Path(partoutput)
            exit_records = None
            if exit_records_cache is not None:
                cache_key = partoutput_path.resolve()
                if cache_key not in exit_records_cache:
                    exit_records_cache[cache_key] = load_daily_exit_records(partoutput_path)
                exit_records = exit_records_cache[cache_key]

            if args.exit_csv is not None:
                exit_csv = args.exit_csv
            elif args.no_exit_csv:
                exit_csv = None
            else:
                exit_csv = out_file.with_name(out_file.stem + "_domain_exit_points.csv")

            add_daily_exit_diagnostics(
                out,
                partoutput=partoutput_path,
                exit_csv=exit_csv,
                exit_records=exit_records,
                particle_offset=particle_offset,
            )

        if args.no_compression:
            encoding = None
        else:
            encoding = _set_netcdf_compression(out, compression_level=args.compression_level)
        out.to_netcdf(out_file, encoding=encoding)
    finally:
        out.close()

    print(f"WROTE {out_file}")
    return out_file


def _discover_from_run_dir(run_dir: Path) -> list[Path]:
    output_dir = run_dir / "output"
    if not output_dir.is_dir():
        output_dir = run_dir
    return sorted(output_dir.glob("grid_time_*.nc"))


def _discover_from_root(root_dir: Path) -> list[Path]:
    return sorted(root_dir.glob("**/output/grid_time_*.nc"))


def _complete_month_filter(grid_files: list[Path]) -> list[Path]:
    by_month: dict[tuple[str, str, int, int], set[dt.date]] = {}
    keep_keys: set[tuple[str, str, int, int]] = set()

    for grid_file in grid_files:
        parsed = _parse_run_dir(grid_file.parent.parent)
        if parsed is None:
            continue
        domain, site, run_date = parsed
        key = (domain, site, run_date.year, run_date.month)
        by_month.setdefault(key, set()).add(run_date)

    for key, dates in by_month.items():
        _, _, year, month = key
        expected = calendar.monthrange(year, month)[1]
        if len(dates) == expected:
            keep_keys.add(key)
        else:
            domain, site, _, _ = key
            print(
                f"SKIP incomplete month {domain}_{site} {year}-{month:02d}: "
                f"{len(dates)} days found, expected {expected}"
            )

    kept = []
    for grid_file in grid_files:
        parsed = _parse_run_dir(grid_file.parent.parent)
        if parsed is None:
            continue
        domain, site, run_date = parsed
        key = (domain, site, run_date.year, run_date.month)
        if key in keep_keys:
            kept.append(grid_file)
    return kept


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Postprocess daily FLEXPART grid_time files containing 24 hourly releases.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=HELP_EPILOG,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--grid-file", type=Path, help="Single grid_time_*.nc file to postprocess")
    source.add_argument("--run-dir", type=Path, help="Daily run directory, usually DOMAIN_SITE_YYYYMMDD")
    source.add_argument("--root-dir", type=Path, help="Root directory containing daily run directories")

    parser.add_argument("--out-file", type=Path, default=None, help="Output file for --grid-file mode")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated files (default: alongside each grid file)",
    )
    parser.add_argument("--domain", default=None, help="Domain name if it cannot be inferred")
    parser.add_argument("--site", default=None, help="Site/receptor code if it cannot be inferred")
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
    parser.add_argument(
        "--partoutput",
        type=Path,
        default=None,
        help="Partoutput file, directory, or glob for domain-exit analysis.",
    )
    parser.add_argument(
        "--auto-partoutput",
        action="store_true",
        help="Look for partoutput_*.nc beside each grid file and add domain-exit diagnostics.",
    )
    parser.add_argument(
        "--exit-csv",
        type=Path,
        default=None,
        help="Optional CSV path for individual daily exit points; only valid with one grid file.",
    )
    parser.add_argument(
        "--no-exit-csv",
        action="store_true",
        help="Add particle location variables without writing the individual exit-point CSV.",
    )
    parser.add_argument("--compression-level", type=int, default=4)
    parser.add_argument("--no-compression", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int, default=0, help="Maximum number of files to process; 0 means all")
    parser.add_argument(
        "--complete-months-only",
        action="store_true",
        help="In --root-dir mode, process only domain/site/month groups with all days present",
    )
    args = parser.parse_args()

    if args.out_file is not None and args.grid_file is None:
        raise ValueError("--out-file can only be used with --grid-file")

    if args.grid_file is not None:
        grid_files = [args.grid_file]
    elif args.run_dir is not None:
        grid_files = _discover_from_run_dir(args.run_dir)
    else:
        grid_files = _discover_from_root(args.root_dir)
        if args.complete_months_only:
            grid_files = _complete_month_filter(grid_files)

    grid_files = [path for path in grid_files if path.name.startswith("grid_time_")]
    if args.limit > 0:
        grid_files = grid_files[: args.limit]

    if args.exit_csv is not None and len(grid_files) != 1:
        raise ValueError("--exit-csv can only be used when processing exactly one grid file")

    if not grid_files:
        print("No grid_time_*.nc files found")
        return 0

    print(f"Discovered {len(grid_files)} grid_time file(s)")
    exit_records_cache: dict[Path, tuple[list[tuple[int, float, float, float, float]], int] | None] = {}
    particle_offsets: dict[Path, int] = {}
    ok = 0
    fail = 0
    for grid_file in grid_files:
        try:
            output_dir = args.output_dir
            partoutput_key = (args.partoutput if args.partoutput is not None else grid_file.parent).resolve()
            particle_offset = particle_offsets.get(partoutput_key, 0)
            out_path = process_grid_file(
                grid_file,
                args,
                output_dir=output_dir,
                exit_records_cache=exit_records_cache,
                particle_offset=particle_offset,
            )
            if out_path is not None and (args.partoutput is not None or args.auto_partoutput):
                processed = xr.open_dataset(out_path, decode_times=False)
                try:
                    particle_offsets[partoutput_key] = particle_offset + _daily_release_particle_count(processed)
                finally:
                    processed.close()
            ok += 1
        except Exception as exc:
            fail += 1
            print(f"FAIL {grid_file}: {exc}")

    print(f"Summary: success={ok}, fail={fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
