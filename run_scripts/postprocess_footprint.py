#!/usr/bin/env python3
"""
Postprocess FLEXPART backward output into practical footprint products.

Features:
- Reads gridded backward sensitivity from grid_time_*.nc
- Creates time-integrated 2D column footprint (all heights)
- Optionally creates time-integrated 2D low-level footprint up to N m agl
- Optionally derives domain-exit locations from partoutput_*.nc particle data

Notes:
- Backward gridded output (grid_time_*.nc) is already the source-receptor
  sensitivity / footprint quantity in FLEXPART.
- Domain-exit locations require particle output from FLEXPART (IPOUT=1 or 2).
"""

import argparse
import csv
import datetime as dt
import glob
import os
import sys

import numpy as np
import xarray as xr


HEIGHT_BINS_M_AGL = np.array([
    500.0, 1500.0, 2500.0, 3500.0, 4500.0,
    5500.0, 6500.0, 7500.0, 8500.0, 9500.0,
    10500.0, 11500.0, 12500.0, 13500.0, 14500.0,
    15500.0, 16500.0, 17500.0, 18500.0, 19500.0,
], dtype=float)

MOLAR_MASS_AIR_KG_PER_MOL = 0.02897


def _open_dataset_auto(path):
    """Open NetCDF with whatever xarray backend is available in the environment."""
    try:
        return xr.open_dataset(path, engine="netcdf4", decode_times=False)
    except Exception:
        return xr.open_dataset(path, decode_times=False)


def _pick_sensitivity_var(ds):
    """Pick the main backward sensitivity variable from grid_time file."""
    preferred_dims = set(["time", "latitude", "longitude"])
    candidates = []

    for name, da in ds.data_vars.items():
        dims = set(da.dims)
        if preferred_dims.issubset(dims):
            candidates.append(name)

    if not candidates:
        raise ValueError("No variable found with dims including time/latitude/longitude")

    for name in candidates:
        units = str(ds[name].attrs.get("units", "")).strip().lower()
        if units == "s":
            return name

    for name in candidates:
        if name.lower().endswith("_mr"):
            return name

    return candidates[0]


def _sum_dims(da, dims):
    active = [d for d in dims if d in da.dims]
    if not active:
        return da
    return da.sum(dim=active, skipna=True)


def _compute_footprints(ds, var_name, lowest_magl=None):
    """
    Compute 2D footprints from sensitivity variable.

    Returns:
      column_2d: sum over time and vertical dimensions
      low_2d: same, but with height <= lowest_magl (or None)
    """
    da = ds[var_name]
    reduce_dims = ["time", "height", "nageclass", "pointspec", "numspec"]

    column_2d = _sum_dims(da, reduce_dims)

    low_2d = None
    if lowest_magl is not None and "height" in da.dims and "height" in ds.coords:
        mask = ds["height"] <= float(lowest_magl)
        if bool(mask.any()):
            da_low = da.where(mask, other=0.0)
            low_2d = _sum_dims(da_low, reduce_dims)
        else:
            # If no layer center falls below threshold, use the lowest layer.
            hmin = float(ds["height"].values.min())
            da_low = da.isel(height=0)
            low_2d = _sum_dims(da_low, ["time", "nageclass", "pointspec", "numspec"])
            print(
                "WARNING: no height level <= {:.1f} m. Using lowest model layer center ({:.1f} m).".format(
                    float(lowest_magl), hmin
                )
            )

    return column_2d, low_2d


def _compute_srr_timeint_2d(ds, var_name):
    """Compute time-integrated 2D SRR by summing time and non-horizontal dims."""
    da = ds[var_name]
    reduce_dims = ["time", "height", "nageclass", "pointspec", "numspec"]
    return _sum_dims(da, reduce_dims)


def _infer_release_time_value(ds):
    """Infer instantaneous release-time coordinate value in native time units."""
    if "RELSTART" in ds:
        relstart = np.asarray(np.ravel(ds["RELSTART"].values), dtype=float)
        finite = np.isfinite(relstart)
        if np.any(finite):
            return float(relstart[finite][0])

    if "time" not in ds or ds["time"].size == 0:
        return 0.0
    tvals = np.asarray(ds["time"].values, dtype=float)
    finite = np.isfinite(tvals)
    if not np.any(finite):
        return 0.0
    tvals = tvals[finite]
    return float(tvals[np.argmin(np.abs(tvals))])


def _convert_to_m2s_per_mol(da_2d, native_units, source_layer_thickness_m):
    """Convert backward SRR from s m3 kg-1 to m2 s mol-1 using source-layer thickness."""
    units = str(native_units or "").strip().lower()
    if units != "s m3 kg-1":
        raise ValueError(
            "Expected source SRR units 's m3 kg-1' for conversion, got '{}'".format(native_units)
        )

    h = float(source_layer_thickness_m)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("source layer thickness must be a positive finite number (m)")

    factor = MOLAR_MASS_AIR_KG_PER_MOL / h
    return da_2d * factor, factor


def _build_time_attrs(ds):
    """Build AGAGE-like time metadata from FLEXPART time coordinate."""
    attrs = {
        "long_name": "time",
        "standard_name": "time",
        "label": "left",
        "comment": "time stamp corresponds to the beginning of each averaging period",
    }

    time_units = str(ds["time"].attrs.get("units", "")) if "time" in ds else ""
    if time_units:
        attrs["units"] = time_units
    calendar = str(ds["time"].attrs.get("calendar", "gregorian")) if "time" in ds else "gregorian"
    attrs["calendar"] = calendar

    lout = ds.attrs.get("loutaver", ds.attrs.get("loutstep", None))
    try:
        period_hours = abs(float(lout)) / 3600.0
    except Exception:
        period_hours = 1.0
    attrs["period"] = f"{period_hours:.1f} hours"
    return attrs


def _append_agage_style_variables(out, ds, site, domain, species, model, met_model):
    """Append AGAGE-like variables and metadata for compatibility with existing tooling."""
    ntime = int(out.sizes.get("time", 0))

    rel_lon = np.nan
    rel_lat = np.nan
    rel_h = np.nan
    if "RELLNG1" in ds:
        rel_lon = float(np.ravel(ds["RELLNG1"].values)[0])
    if "RELLAT1" in ds:
        rel_lat = float(np.ravel(ds["RELLAT1"].values)[0])
    if "RELZZ1" in ds:
        rel_h = float(np.ravel(ds["RELZZ1"].values)[0])

    out["release_lon"] = xr.DataArray(
        np.full(ntime, rel_lon, dtype=np.float32),
        dims=("time",),
        coords={"time": out["time"]},
        attrs={"units": "degrees_east", "long_name": "Release longitude"},
    )
    out["release_lat"] = xr.DataArray(
        np.full(ntime, rel_lat, dtype=np.float32),
        dims=("time",),
        coords={"time": out["time"]},
        attrs={"units": "degrees_north", "long_name": "Release latitude"},
    )
    out["release_height"] = xr.DataArray(
        np.full(ntime, rel_h, dtype=np.float32),
        dims=("time",),
        coords={"time": out["time"]},
        attrs={"units": "m", "long_name": "Release height above model ground"},
    )

    # Present in AGAGE files; fill with NaN if unavailable from FLEXPART grid output.
    missing_series = [
        ("air_temperature", "K", "air temperature at release"),
        ("air_pressure", "hPa", "air pressure at release"),
        ("wind_speed", "m s-1", "wind speed at release"),
        ("wind_from_direction", "degree", "wind direction at release"),
        ("atmosphere_boundary_layer_thickness", "m", "atmospheric boundary layer thickness at release"),
    ]
    for name, units, long_name in missing_series:
        out[name] = xr.DataArray(
            np.full(ntime, np.nan, dtype=np.float32),
            dims=("time",),
            coords={"time": out["time"]},
            attrs={"units": units, "long_name": long_name},
        )

    out.attrs.update({
        "lpdm_native_output_unit": "s",
        "species": species,
        "model": model,
        "met_model": met_model,
        "output_folder": "output",
        "model_version": str(ds.attrs.get("source", "FLEXPART")),
        "domain": domain,
        "site": site,
        "author": "FLEXPART postprocess_footprint.py",
        "created": dt.datetime.utcnow().isoformat() + "Z",
    })


def _find_particle_vars(ds):
    """Find longitude, latitude, and height particle variables in partoutput dataset."""
    lon_name = None
    lat_name = None
    z_name = None

    for name, da in ds.data_vars.items():
        attrs = {k.lower(): str(v).lower() for k, v in da.attrs.items()}
        txt = " ".join([
            attrs.get("standard_name", ""),
            attrs.get("long_name", ""),
            attrs.get("description", ""),
            attrs.get("axis", ""),
            attrs.get("units", ""),
            name.lower(),
        ])
        dims = [d.lower() for d in da.dims]
        looks_particle_field = any("part" in d for d in dims) and any("time" in d for d in dims)

        if not looks_particle_field:
            continue

        if lon_name is None and ("longitude" in txt or "degrees_east" in txt):
            lon_name = name
        if lat_name is None and ("latitude" in txt or "degrees_north" in txt):
            lat_name = name
        if z_name is None and ("height" in txt or "altitude" in txt or " m" in txt):
            z_name = name

    return lon_name, lat_name, z_name


def _time_particle_array(da):
    """Return array with shape (time, particle) and numeric values."""
    dims = list(da.dims)
    tdim = None
    pdim = None
    for d in dims:
        dl = d.lower()
        if "time" in dl:
            tdim = d
        if "part" in dl:
            pdim = d

    if tdim is None or pdim is None:
        raise ValueError("Could not identify time/particle dims for variable {}".format(da.name))

    arr = da.transpose(tdim, pdim).values
    return np.asarray(arr, dtype=float), tdim, pdim


def _derive_exit_points(part_ds, lon_var, lat_var, z_var=None):
    """
    Derive particle exit points from first NaN transition in particle tracks.

    Assumes terminated particles are represented by NaN fields (default FLEXPART
    behavior when IPOUT>0).
    """
    lon, tdim, _ = _time_particle_array(part_ds[lon_var])
    lat, _, _ = _time_particle_array(part_ds[lat_var])
    z = None
    if z_var is not None and z_var in part_ds:
        z, _, _ = _time_particle_array(part_ds[z_var])

    if lon.shape != lat.shape:
        raise ValueError("Longitude and latitude arrays have different shapes")
    if z is not None and z.shape != lon.shape:
        raise ValueError("Height array shape does not match lon/lat arrays")

    if tdim in part_ds:
        tvals = np.asarray(part_ds[tdim].values)
    else:
        tvals = np.arange(lon.shape[0], dtype=float)

    exits = []
    ntime, npart = lon.shape
    for p in range(npart):
        valid = np.isfinite(lon[:, p]) & np.isfinite(lat[:, p])
        if not np.any(valid):
            continue

        invalid_idx = np.where(~valid)[0]
        if invalid_idx.size == 0:
            # For IPOUT=2 (end-only dumps), use the last valid point as exit location.
            i_prev = int(np.where(valid)[0][-1])
        else:
            i_prev = int(invalid_idx[0]) - 1
            if i_prev < 0:
                continue
        if not np.isfinite(lon[i_prev, p]) or not np.isfinite(lat[i_prev, p]):
            continue

        z_val = np.nan
        if z is not None and np.isfinite(z[i_prev, p]):
            z_val = float(z[i_prev, p])

        exits.append((
            p,
            float(tvals[i_prev]) if i_prev < len(tvals) else float(i_prev),
            float(lon[i_prev, p]),
            float(lat[i_prev, p]),
            z_val,
        ))

    return exits, int(npart)


def _nearest_index(values, value):
    arr = np.asarray(values, dtype=float)
    return int(np.argmin(np.abs(arr - float(value))))


def _classify_exit_side(lon, lat, lon_min, lon_max, lat_min, lat_max):
    dist = {
        "w": abs(float(lon) - float(lon_min)),
        "e": abs(float(lon) - float(lon_max)),
        "s": abs(float(lat) - float(lat_min)),
        "n": abs(float(lat) - float(lat_max)),
    }
    return min(dist, key=dist.get)


def _build_boundary_exit_fractions(exits, time_vals, lon_centers, lat_centers, height_centers):
    nt = len(time_vals)
    nh = len(height_centers)
    nlon = len(lon_centers)
    nlat = len(lat_centers)

    n_counts = np.zeros((nt, nh, nlon), dtype=np.float32)
    s_counts = np.zeros((nt, nh, nlon), dtype=np.float32)
    e_counts = np.zeros((nt, nh, nlat), dtype=np.float32)
    w_counts = np.zeros((nt, nh, nlat), dtype=np.float32)

    lon_min = float(np.min(lon_centers))
    lon_max = float(np.max(lon_centers))
    lat_min = float(np.min(lat_centers))
    lat_max = float(np.max(lat_centers))

    for _, tval, lon, lat, z in exits:
        tidx = _nearest_index(time_vals, tval)
        hval = z if np.isfinite(z) else float(height_centers[0])
        hidx = _nearest_index(height_centers, hval)

        side = _classify_exit_side(lon, lat, lon_min, lon_max, lat_min, lat_max)
        if side == "n":
            xidx = _nearest_index(lon_centers, lon)
            n_counts[tidx, hidx, xidx] += 1.0
        elif side == "s":
            xidx = _nearest_index(lon_centers, lon)
            s_counts[tidx, hidx, xidx] += 1.0
        elif side == "e":
            yidx = _nearest_index(lat_centers, lat)
            e_counts[tidx, hidx, yidx] += 1.0
        else:
            yidx = _nearest_index(lat_centers, lat)
            w_counts[tidx, hidx, yidx] += 1.0

    # Normalize by exiting particles so the sum over N/E/S/W boundary fractions is 1.
    denom = float(len(exits)) if len(exits) > 0 else 1.0
    return n_counts / denom, e_counts / denom, s_counts / denom, w_counts / denom


def _centers_to_edges(centers):
    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 1 or centers.size < 2:
        raise ValueError("Need at least 2 coordinate centers to infer edges")
    d = np.diff(centers)
    first = centers[0] - d[0] / 2.0
    last = centers[-1] + d[-1] / 2.0
    mid = centers[:-1] + d / 2.0
    return np.concatenate([[first], mid, [last]])


def _build_exit_histogram(exits, lon_centers, lat_centers):
    if not exits:
        return np.zeros((len(lat_centers), len(lon_centers)), dtype=np.int32)

    lons = np.array([e[2] for e in exits], dtype=float)
    lats = np.array([e[3] for e in exits], dtype=float)

    lon_edges = _centers_to_edges(lon_centers)
    lat_edges = _centers_to_edges(lat_centers)

    hist, _, _ = np.histogram2d(lats, lons, bins=[lat_edges, lon_edges])
    return hist.astype(np.int32)


def _open_partoutput(path_arg):
    """Open one or more partoutput NetCDF files and concatenate on time."""
    if path_arg is None:
        return None

    paths = []
    if os.path.isdir(path_arg):
        paths = sorted(glob.glob(os.path.join(path_arg, "partoutput_*.nc")))
    elif any(ch in path_arg for ch in ["*", "?", "["]):
        paths = sorted(glob.glob(path_arg))
    elif os.path.isfile(path_arg):
        paths = [path_arg]

    if not paths:
        return None

    if len(paths) == 1:
        return _open_dataset_auto(paths[0])

    datasets = [_open_dataset_auto(p) for p in paths]
    try:
        merged = xr.concat(datasets, dim="time")
    except Exception:
        for ds in datasets:
            ds.close()
        raise
    for ds in datasets:
        ds.close()
    return merged


def main():
    parser = argparse.ArgumentParser(description="Postprocess FLEXPART backward footprints")
    parser.add_argument("--grid-file", required=True, help="Path to grid_time_*.nc")
    parser.add_argument(
        "--out-file",
        default=None,
        help="Output NetCDF for derived footprints (default: <grid-file stem>_footprints.nc)",
    )
    parser.add_argument(
        "--lowest-magl",
        type=float,
        default=100.0,
        help="Optional low-level footprint top in m agl (default: 100). Set <0 to disable.",
    )
    parser.add_argument(
        "--partoutput",
        default=None,
        help="Optional partoutput file, directory, or glob for domain-exit analysis",
    )
    parser.add_argument(
        "--exit-csv",
        default=None,
        help="Optional CSV path for individual exit points (default: alongside out-file)",
    )
    parser.add_argument("--site", default="UNKNOWN", help="Site/receptor code for AGAGE-like metadata")
    parser.add_argument("--domain", default="UNKNOWN", help="Domain name for AGAGE-like metadata")
    parser.add_argument("--species", default="inert", help="Species label for AGAGE-like metadata")
    parser.add_argument("--model", default="FLEXPART", help="LPDM model label for AGAGE-like metadata")
    parser.add_argument("--met-model", default="GFS", help="Meteorological model label for AGAGE-like metadata")
    parser.add_argument(
        "--source-layer-thickness-m",
        type=float,
        default=100.0,
        help="Source-layer thickness (m) for converting SRR to m2 s mol-1 (default: 100).",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.grid_file):
        raise FileNotFoundError("grid file not found: {}".format(args.grid_file))

    out_file = args.out_file
    if out_file is None:
        base, _ = os.path.splitext(args.grid_file)
        out_file = base + "_footprints.nc"

    print("Reading grid file: {}".format(args.grid_file))
    ds = _open_dataset_auto(args.grid_file)
    try:
        var_name = _pick_sensitivity_var(ds)
        print("Using sensitivity variable: {}".format(var_name))

        lowest = args.lowest_magl if args.lowest_magl is not None and args.lowest_magl >= 0 else None
        col_2d, low_2d = _compute_footprints(ds, var_name, lowest_magl=lowest)
        srr_timeint_2d = _compute_srr_timeint_2d(ds, var_name)
        native_units = ds[var_name].attrs.get("units", "")
        srr_timeint_2d, conv_factor = _convert_to_m2s_per_mol(
            srr_timeint_2d,
            native_units,
            args.source_layer_thickness_m,
        )
        col_2d, _ = _convert_to_m2s_per_mol(
            col_2d,
            native_units,
            args.source_layer_thickness_m,
        )
        if low_2d is not None:
            low_2d, _ = _convert_to_m2s_per_mol(
                low_2d,
                native_units,
                args.source_layer_thickness_m,
            )
        release_time_value = _infer_release_time_value(ds)

        out = xr.Dataset()
        out["time"] = xr.DataArray(
            np.array([release_time_value], dtype=np.float64),
            dims=("time",),
        )
        out["latitude"] = ds["latitude"]
        out["longitude"] = ds["longitude"]

        out["srr"] = srr_timeint_2d.astype(np.float32).expand_dims({"time": out["time"]}).transpose("time", "latitude", "longitude")
        out["srr"].attrs.update({
            "long_name": "source_receptor_relationship",
            "loss_lifetime_hrs": -9.0,
            "loss_lifetime_comment": "lifetime in hours; -9 corresponds to inert",
            "units": "m2 s mol-1",
            "source_variable": var_name,
            "description": "time-integrated footprint at receptor release time (converted to molar flux sensitivity)",
            "conversion_from_native_units": str(native_units),
            "conversion_factor_applied": float(conv_factor),
            "source_layer_thickness_m": float(args.source_layer_thickness_m),
            "molar_mass_air_kg_per_mol": float(MOLAR_MASS_AIR_KG_PER_MOL),
        })

        # Intentionally keep only `srr` as the gridded footprint field in output.

        # Optional domain-exit diagnostics from particle output.
        pds = _open_partoutput(args.partoutput)
        try:
            if pds is None:
                print("No partoutput files supplied/found; skipping domain-exit diagnostics.")
            else:
                lon_var, lat_var, z_var = _find_particle_vars(pds)
                if lon_var is None or lat_var is None:
                    print("Could not find particle longitude/latitude vars; skipping domain exits.")
                else:
                    exits, npart = _derive_exit_points(pds, lon_var, lat_var, z_var=z_var)
                    print("Derived {} particle exit points".format(len(exits)))

                    lon_centers = np.asarray(ds["longitude"].values, dtype=float)
                    lat_centers = np.asarray(ds["latitude"].values, dtype=float)
                    exit_count = _build_exit_histogram(exits, lon_centers, lat_centers)

                    out["height"] = xr.DataArray(
                        HEIGHT_BINS_M_AGL.astype(np.float32),
                        dims=("height",),
                        attrs={
                            "long_name": "height at layer midpoints",
                            "units": "m",
                            "positive": "up",
                        },
                    )

                    time_vals = np.asarray(out["time"].values, dtype=float)
                    frac_n, frac_e, frac_s, frac_w = _build_boundary_exit_fractions(
                        exits,
                        time_vals=time_vals,
                        lon_centers=lon_centers,
                        lat_centers=lat_centers,
                        height_centers=HEIGHT_BINS_M_AGL,
                    )

                    out["particle_locations_n"] = xr.DataArray(
                        frac_n,
                        dims=("time", "height", "longitude"),
                        coords={"time": out["time"], "height": out["height"], "longitude": out["longitude"]},
                        attrs={
                            "long_name": "Fraction of exiting particles leaving domain (N side)",
                            "units": "1",
                        },
                    )
                    out["particle_locations_e"] = xr.DataArray(
                        frac_e,
                        dims=("time", "height", "latitude"),
                        coords={"time": out["time"], "height": out["height"], "latitude": out["latitude"]},
                        attrs={
                            "long_name": "Fraction of exiting particles leaving domain (E side)",
                            "units": "1",
                        },
                    )
                    out["particle_locations_s"] = xr.DataArray(
                        frac_s,
                        dims=("time", "height", "longitude"),
                        coords={"time": out["time"], "height": out["height"], "longitude": out["longitude"]},
                        attrs={
                            "long_name": "Fraction of exiting particles leaving domain (S side)",
                            "units": "1",
                        },
                    )
                    out["particle_locations_w"] = xr.DataArray(
                        frac_w,
                        dims=("time", "height", "latitude"),
                        coords={"time": out["time"], "height": out["height"], "latitude": out["latitude"]},
                        attrs={
                            "long_name": "Fraction of exiting particles leaving domain (W side)",
                            "units": "1",
                        },
                    )

                    exit_csv = args.exit_csv
                    if exit_csv is None:
                        base, _ = os.path.splitext(out_file)
                        exit_csv = base + "_domain_exit_points.csv"

                    with open(exit_csv, "w", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow([
                            "particle_index",
                            "exit_time",
                            "exit_longitude",
                            "exit_latitude",
                            "exit_height_magl",
                            "exit_side",
                        ])
                        lon_min = float(np.min(lon_centers))
                        lon_max = float(np.max(lon_centers))
                        lat_min = float(np.min(lat_centers))
                        lat_max = float(np.max(lat_centers))
                        for pidx, tval, lon, lat, z in exits:
                            side = _classify_exit_side(lon, lat, lon_min, lon_max, lat_min, lat_max)
                            writer.writerow([pidx, tval, lon, lat, z, side])
                    print("Wrote domain exit points CSV: {}".format(exit_csv))
        finally:
            if pds is not None:
                pds.close()

        out.attrs.update({
            "title": "Derived FLEXPART footprint products",
            "input_grid_file": os.path.abspath(args.grid_file),
            "note": "Backward grid_time is already source-receptor sensitivity; this file is postprocessed summaries.",
        })

        out["time"].attrs.update(_build_time_attrs(ds))
        out["time"].attrs["comment"] = "single instantaneous release-time stamp for time-integrated footprint"
        out["latitude"].attrs.update({"units": "degrees_north", "long_name": "latitude"})
        out["longitude"].attrs.update({"units": "degrees_east", "long_name": "longitude"})

        _append_agage_style_variables(
            out,
            ds,
            site=args.site,
            domain=args.domain,
            species=args.species,
            model=args.model,
            met_model=args.met_model,
        )

        out.to_netcdf(out_file)
        print("Wrote footprint products: {}".format(out_file))
    finally:
        ds.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
