#!/usr/bin/env python3
"""
Postprocess FLEXPART backward output into practical footprint products.

Features:
- Reads gridded backward sensitivity from grid_time_*.nc
- Creates time-integrated 2D footprint from a selected output height level
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
import importlib
import os
import re
import sys

import numpy as np
import xarray as xr


HEIGHT_BINS_M_AGL = np.array([
    500.0, 1500.0, 2500.0, 3500.0, 4500.0,
    5500.0, 6500.0, 7500.0, 8500.0, 9500.0,
    10500.0, 11500.0, 12500.0, 13500.0, 14500.0,
    15500.0, 16500.0, 17500.0, 18500.0, 19500.0,
], dtype=float)

MOLAR_MASS_AIR_KG_PER_MOL = 29./1000.
DEFAULT_CFSV2_GRIB_DIR = "/net/fs01/data/AGAGE/meteorology/cfsv2/flexpart_inputs"


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


def _find_spatial_dims(ds):
    """Robustly find latitude and longitude dimension names."""
    lat_names = ["latitude", "lat", "y"]
    lon_names = ["longitude", "lon", "x"]
    
    lat_dim = next((n for n in lat_names if n in ds.dims), None)
    lon_dim = next((n for n in lon_names if n in ds.dims), None)
    
    if lat_dim is None or lon_dim is None:
        raise ValueError(
            f"Could not identify lat/lon dimensions in {list(ds.dims)}. "
            f"Expected one of {lat_names} and {lon_names}"
        )
    
    return lat_dim, lon_dim


def _sum_dims(da, dims):
    active = [d for d in dims if d in da.dims]
    if not active:
        return da
    return da.sum(dim=active, skipna=True)



def _compute_srr_timeint_2d(ds, var_name, footprint_outheight_m=100.0):
    """
    Compute time-integrated 2D SRR cumulatively up to a requested outheight.

    If a height dimension exists, the nearest level to footprint_outheight_m is
    selected as the top integration level and all lower levels are summed.
    """
    da = ds[var_name]
    selected_height = None
    integrated_height_levels = None

    if "height" in da.dims:
        if "height" in da.coords:
            hvals = np.asarray(da["height"].values, dtype=float)
        elif "height" in ds:
            hvals = np.asarray(ds["height"].values, dtype=float)
        else:
            raise ValueError("Variable has 'height' dimension but no usable height coordinate")

        if hvals.size == 0:
            raise ValueError("Height coordinate is empty; cannot select footprint outheight")

        hidx = int(np.argmin(np.abs(hvals - float(footprint_outheight_m))))
        selected_height = float(hvals[hidx])
        # Integrate all layers from the surface up to the selected top outheight.
        da = da.isel(height=slice(0, hidx + 1))
        integrated_height_levels = int(hidx + 1)

        if not np.isclose(selected_height, float(footprint_outheight_m), atol=1e-6):
            print(
                "WARNING: requested footprint outheight {:.3f} m not found; using nearest {:.3f} m".format(
                    float(footprint_outheight_m), selected_height
                )
            )

    reduce_dims = ["time", "height", "nageclass", "pointspec", "numspec"]
    return _sum_dims(da, reduce_dims), selected_height, integrated_height_levels


def _infer_release_time_value(ds):
    """Infer instantaneous release-time coordinate value in native time units."""
    if "RELSTART" in ds:
        relstart = np.asarray(np.ravel(ds["RELSTART"].values), dtype=float)
        finite = np.isfinite(relstart)
        if np.any(finite):
            return float(relstart[finite][0])

    if "time" not in ds or ds["time"].size == 0:
        return None
    tvals = np.asarray(ds["time"].values, dtype=float)
    finite = np.isfinite(tvals)
    if not np.any(finite):
        return None
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


def _append_agage_style_variables(out, ds, site, domain, species, model, met_model, release_meteo=None):
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

    # Present in AGAGE files; fill with extracted release-time meteo where available.
    missing_series = [
        ("air_temperature", "K", "air temperature at release"),
        ("air_pressure", "hPa", "air pressure at release"),
        ("wind_speed", "m s-1", "wind speed at release"),
        ("wind_from_direction", "degree", "wind direction at release"),
        ("atmosphere_boundary_layer_thickness", "m", "atmospheric boundary layer thickness at release"),
    ]
    release_meteo = release_meteo or {}
    for name, units, long_name in missing_series:
        value = release_meteo.get(name, np.nan)
        out[name] = xr.DataArray(
            np.full(ntime, value, dtype=np.float32),
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


def _parse_release_datetime_from_grid_filename(path):
    """Parse release datetime from grid_time_YYYYMMDDHHMMSS.nc filename."""
    m = re.match(r"^grid_time_(\d{14})\.nc$", os.path.basename(path))
    if not m:
        return None
    try:
        return dt.datetime.strptime(m.group(1), "%Y%m%d%H%M%S")
    except ValueError:
        return None


def _infer_release_location(ds):
    """Infer release longitude/latitude from grid metadata."""
    rel_lon = np.nan
    rel_lat = np.nan
    if "RELLNG1" in ds:
        rel_lon = float(np.ravel(ds["RELLNG1"].values)[0])
    if "RELLAT1" in ds:
        rel_lat = float(np.ravel(ds["RELLAT1"].values)[0])
    return rel_lon, rel_lat


def _find_gf_file_for_time(grib_dir, release_dt):
    """Find nearest GFyymmddhh file to the release timestamp."""
    if release_dt is None or grib_dir is None or not os.path.isdir(grib_dir):
        return None

    best_path = None
    best_delta = None
    for name in os.listdir(grib_dir):
        m = re.match(r"^GF(\d{8})$", name)
        if not m:
            continue
        try:
            ts = dt.datetime.strptime(m.group(1), "%y%m%d%H")
        except ValueError:
            continue
        delta = abs((ts - release_dt).total_seconds())
        if best_delta is None or delta < best_delta:
            best_delta = delta
            best_path = os.path.join(grib_dir, name)

    return best_path


def _normalize_lon_to_grid(lon_value, lon_grid):
    lon = float(lon_value)
    lon_vals = np.asarray(lon_grid, dtype=float)
    if lon_vals.size == 0:
        return lon

    lon_min = float(np.nanmin(lon_vals))
    lon_max = float(np.nanmax(lon_vals))

    if lon_min >= 0.0 and lon_max > 180.0 and lon < 0.0:
        return lon % 360.0
    if lon_max <= 180.0 and lon > 180.0:
        return ((lon + 180.0) % 360.0) - 180.0
    return lon


def _select_nearest_valid_time(da, target_dt64):
    """Select nearest valid-time slice when valid_time/time coordinates exist."""
    if target_dt64 is None:
        return da

    if "valid_time" in da.coords:
        vt = da["valid_time"]
        vt_vals = np.asarray(vt.values)
        if vt_vals.size > 1:
            vt_s = vt_vals.astype("datetime64[s]")
            idx_flat = int(np.argmin(np.abs(vt_s - target_dt64)))
            idx_multi = np.unravel_index(idx_flat, vt_vals.shape)
            indexers = {dim: int(i) for dim, i in zip(vt.dims, idx_multi)}
            return da.isel(indexers)

    if "time" in da.dims and "time" in da.coords:
        tvals = np.asarray(da["time"].values)
        if tvals.size > 1 and np.issubdtype(tvals.dtype, np.datetime64):
            t_s = tvals.astype("datetime64[s]")
            tidx = int(np.argmin(np.abs(t_s - target_dt64)))
            da = da.isel(time=tidx)

    if "step" in da.dims and da.sizes.get("step", 1) > 1:
        da = da.isel(step=0)

    return da


def _extract_point_value(da, release_dt64, rel_lon, rel_lat):
    """Extract nearest scalar value from a gridded field at release time/location."""
    da = _select_nearest_valid_time(da, release_dt64)

    lat_name = next((n for n in ["latitude", "lat", "y"] if n in da.coords or n in da.dims), None)
    lon_name = next((n for n in ["longitude", "lon", "x"] if n in da.coords or n in da.dims), None)
    if lat_name is None or lon_name is None:
        return np.nan

    lon_target = _normalize_lon_to_grid(rel_lon, da[lon_name].values)
    try:
        da = da.sel({lat_name: float(rel_lat), lon_name: float(lon_target)}, method="nearest")
    except Exception:
        lat_idx = _nearest_index(da[lat_name].values, rel_lat)
        lon_idx = _nearest_index(da[lon_name].values, lon_target)
        da = da.isel({lat_name: lat_idx, lon_name: lon_idx})

    remaining = [d for d in da.dims if d not in [lat_name, lon_name]]
    for d in remaining:
        if da.sizes.get(d, 0) > 0:
            da = da.isel({d: 0})

    arr = np.asarray(da.values, dtype=float).ravel()
    if arr.size == 0:
        return np.nan
    return float(arr[0])


def _find_var_ci(datasets, preferred_names):
    """Find first data variable by case-insensitive preferred name order."""
    lname_to_da = {}
    for ds in datasets:
        for name in ds.data_vars:
            lname_to_da.setdefault(name.lower(), ds[name])

    for name in preferred_names:
        da = lname_to_da.get(name.lower())
        if da is not None:
            return da
    return None


def _extract_release_meteo_from_grib(grib_file, release_dt, rel_lon, rel_lat):
    """Extract release-time meteorology from a CFSv2 GRIB file using cfgrib."""
    release_meteo = {}
    if grib_file is None:
        return release_meteo
    if not np.isfinite(rel_lon) or not np.isfinite(rel_lat):
        return release_meteo

    try:
        cfgrib = importlib.import_module("cfgrib")
    except Exception:
        print("WARNING: cfgrib not available; leaving release meteorology as NaN")
        return release_meteo

    try:
        datasets = cfgrib.open_datasets(grib_file)
    except Exception as e:
        print("WARNING: failed to open GRIB file '{}': {}".format(grib_file, e))
        return release_meteo

    if not datasets:
        print("WARNING: no datasets found in GRIB file '{}'".format(grib_file))
        return release_meteo

    release_dt64 = np.datetime64(release_dt, "s") if release_dt is not None else None

    t_da = _find_var_ci(datasets, ["t2m", "2t", "tmp", "t"])
    p_da = _find_var_ci(datasets, ["sp", "pres", "surface_pressure", "msl", "prmsl"])
    u_da = _find_var_ci(datasets, ["u10", "10u", "u"])
    v_da = _find_var_ci(datasets, ["v10", "10v", "v"])
    blh_da = _find_var_ci(datasets, ["blh", "hpbl", "boundary_layer_height"])

    if t_da is not None:
        t_val = _extract_point_value(t_da, release_dt64, rel_lon, rel_lat)
        if np.isfinite(t_val):
            release_meteo["air_temperature"] = t_val

    if p_da is not None:
        p_val = _extract_point_value(p_da, release_dt64, rel_lon, rel_lat)
        if np.isfinite(p_val):
            p_units = str(p_da.attrs.get("units", "")).strip().lower()
            if p_units in ["pa", "pascal", "pascals"] or p_val > 20000.0:
                p_val = p_val / 100.0
            release_meteo["air_pressure"] = p_val

    u_val = np.nan
    v_val = np.nan
    if u_da is not None:
        u_val = _extract_point_value(u_da, release_dt64, rel_lon, rel_lat)
    if v_da is not None:
        v_val = _extract_point_value(v_da, release_dt64, rel_lon, rel_lat)
    if np.isfinite(u_val) and np.isfinite(v_val):
        speed = float(np.hypot(u_val, v_val))
        direction = float((270.0 - np.degrees(np.arctan2(v_val, u_val))) % 360.0)
        release_meteo["wind_speed"] = speed
        release_meteo["wind_from_direction"] = direction

    if blh_da is not None:
        blh_val = _extract_point_value(blh_da, release_dt64, rel_lon, rel_lat)
        if np.isfinite(blh_val):
            release_meteo["atmosphere_boundary_layer_thickness"] = blh_val

    return release_meteo


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
    
    Vectorized for performance: processes all 20k+ particles at once instead of
    looping, achieving ~50-100x speedup over naive per-particle iteration.
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
        tvals = np.asarray(part_ds[tdim].values, dtype=float)
    else:
        tvals = np.arange(lon.shape[0], dtype=float)

    ntime, npart = lon.shape

    valid = np.isfinite(lon) & np.isfinite(lat)
    invalid = ~valid
    has_valid = valid.any(axis=0)
    has_invalid = invalid.any(axis=0)

    first_invalid = np.argmax(invalid, axis=0)
    exit_idx = np.where(has_invalid, first_invalid - 1, ntime - 1)
    usable = has_valid & (exit_idx >= 0)

    particles = np.nonzero(usable)[0]
    exit_idx = exit_idx[usable].astype(np.intp, copy=False)
    finite_exit = valid[exit_idx, particles]
    particles = particles[finite_exit]
    exit_idx = exit_idx[finite_exit]

    exit_times = tvals[exit_idx] if len(tvals) >= ntime else exit_idx.astype(float)
    z_vals = np.full(len(particles), np.nan, dtype=float)
    if z is not None:
        z_at_exit = z[exit_idx, particles]
        finite_z = np.isfinite(z_at_exit)
        z_vals[finite_z] = z_at_exit[finite_z]

    exits = [
        (
            int(p),
            float(t),
            float(x),
            float(y),
            float(zz) if np.isfinite(zz) else np.nan,
        )
        for p, t, x, y, zz in zip(
            particles,
            exit_times,
            lon[exit_idx, particles],
            lat[exit_idx, particles],
            z_vals,
        )
    ]

    return exits, int(npart)


def _load_exit_points_csv(path):
    """Load previously written domain-exit CSV rows into exit tuples."""
    exits = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                pidx = int(row.get("particle_index", ""))
                tval = float(row.get("exit_time", "nan"))
                lon = float(row.get("exit_longitude", "nan"))
                lat = float(row.get("exit_latitude", "nan"))
                z_raw = row.get("exit_height_magl", "nan")
                z = float(z_raw) if z_raw not in [None, "", "nan", "NaN"] else np.nan
            except Exception:
                continue
            if np.isfinite(lon) and np.isfinite(lat):
                exits.append((pidx, tval, lon, lat, z))
    return exits


def _nearest_index(values, value):
    arr = np.asarray(values, dtype=float)
    return int(np.argmin(np.abs(arr - float(value))))


def _nearest_indices(values, targets):
    arr = np.asarray(values, dtype=float)
    targets = np.asarray(targets, dtype=float)
    if arr.size == 0:
        raise ValueError("Cannot locate nearest indices in an empty coordinate")
    if arr.size == 1:
        return np.zeros(targets.shape, dtype=np.intp)

    diffs = np.diff(arr)
    if np.all(diffs > 0):
        idx = np.searchsorted(arr, targets)
        right = np.clip(idx, 0, arr.size - 1)
        left = np.clip(right - 1, 0, arr.size - 1)
    elif np.all(diffs < 0):
        idx = np.searchsorted(-arr, -targets)
        right = np.clip(idx, 0, arr.size - 1)
        left = np.clip(right - 1, 0, arr.size - 1)
    else:
        return np.argmin(np.abs(targets[:, None] - arr[None, :]), axis=1)

    choose_right = np.abs(arr[right] - targets) < np.abs(arr[left] - targets)
    return np.where(choose_right, right, left).astype(np.intp, copy=False)


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

    if exits:
        exit_arr = np.asarray(exits, dtype=float)
        tvals = exit_arr[:, 1]
        lons = exit_arr[:, 2]
        lats = exit_arr[:, 3]
        zvals = exit_arr[:, 4]
        zvals = np.where(np.isfinite(zvals), zvals, float(height_centers[0]))

        tidx = _nearest_indices(time_vals, tvals)
        hidx = _nearest_indices(height_centers, zvals)
        xidx = _nearest_indices(lon_centers, lons)
        yidx = _nearest_indices(lat_centers, lats)

        distances = np.column_stack([
            np.abs(lons - lon_min),
            np.abs(lons - lon_max),
            np.abs(lats - lat_min),
            np.abs(lats - lat_max),
        ])
        sides = np.argmin(distances, axis=1)

        west = sides == 0
        east = sides == 1
        south = sides == 2
        north = sides == 3

        np.add.at(w_counts, (tidx[west], hidx[west], yidx[west]), 1.0)
        np.add.at(e_counts, (tidx[east], hidx[east], yidx[east]), 1.0)
        np.add.at(s_counts, (tidx[south], hidx[south], xidx[south]), 1.0)
        np.add.at(n_counts, (tidx[north], hidx[north], xidx[north]), 1.0)

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


def _open_partoutput(path_arg, particle_start=None, particle_end=None):
    """Open one or more partoutput NetCDF files and concatenate on time.

    particle_start / particle_end: optional 0-based inclusive index range to
    select a slice of particles (used in daily-batch mode to isolate the
    particles belonging to one specific release group).
    """
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
        ds = _open_dataset_auto(paths[0])
    else:
        datasets = [_open_dataset_auto(p) for p in paths]
        try:
            ds = xr.concat(datasets, dim="time")
        except Exception:
            for d in datasets:
                d.close()
            raise
        for d in datasets:
            d.close()

    if particle_start is not None or particle_end is not None:
        part_dim = next((d for d in ds.dims if "part" in d.lower()), None)
        if part_dim is not None:
            s = particle_start if particle_start is not None else 0
            e = (particle_end + 1) if particle_end is not None else ds.dims[part_dim]
            ds = ds.isel({part_dim: slice(s, e)})
            print("Particle filter: indices {}–{} (dim '{}', {} particles)".format(
                s, e - 1, part_dim, e - s))
        else:
            print("WARNING: --particle-start/--particle-end specified but no "
                  "'part*' dimension found in partoutput; using all particles.")

    return ds


def _set_netcdf_compression(ds, compression_level=4):
    """
    Configure zlib compression for all data variables in a dataset.
    
    Args:
        ds: xarray Dataset
        compression_level: zlib compression level (1-9, default 4)
    
    Returns:
        dict of encoding specifications for to_netcdf()
    """
    encoding = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {
            "zlib": True,
            "complevel": compression_level,
        }
    return encoding


def main():
    parser = argparse.ArgumentParser(description="Postprocess FLEXPART backward footprints")
    parser.add_argument("--grid-file", required=True, help="Path to grid_time_*.nc")
    parser.add_argument(
        "--out-file",
        default=None,
        help="Output NetCDF for derived footprints (default: <grid-file stem>_footprints.nc)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file if present.",
    )
    parser.add_argument(
        "--skip-existing",
        dest="skip_existing",
        action="store_true",
        default=True,
        help="Skip processing when output file already exists (default).",
    )
    parser.add_argument(
        "--no-skip-existing",
        dest="skip_existing",
        action="store_false",
        help="Do not skip existing output file; process and replace it.",
    )

    parser.add_argument(
        "--partoutput",
        default=None,
        help="Partoutput file, directory, or glob for domain-exit analysis (default: grid-file directory)",
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
    parser.add_argument("--met-model", default="CFSv2", help="Meteorological model label for AGAGE-like metadata")
    parser.add_argument(
        "--source-layer-thickness-m",
        type=float,
        default=100.0,
        help="Source-layer thickness (m) for converting SRR to m2 s mol-1 (default: 100).",
    )
    parser.add_argument(
        "--footprint-outheight-m",
        type=float,
        default=100.0,
        help="Outheight layer (m agl) used to derive the footprint (default: 100).",
    )
    parser.add_argument(
        "--meteo-grib-dir",
        default=DEFAULT_CFSV2_GRIB_DIR,
        help="Directory containing GFyymmddhh GRIB files for release-time meteorology extraction.",
    )
    parser.add_argument(
        "--meteo-grib-file",
        default=None,
        help="Optional explicit GRIB file path to use for release-time meteorology extraction.",
    )
    parser.add_argument(
        "--disable-grib-meteo",
        action="store_true",
        help="Disable GRIB-based extraction of release-time meteorology.",
    )
    parser.add_argument(
        "--use-existing-exit-csv",
        action="store_true",
        help="Reuse --exit-csv if it exists instead of recomputing domain-exit points from particle output.",
    )
    parser.add_argument(
        "--particle-start",
        type=int,
        default=None,
        metavar="N",
        help="0-based start index (inclusive) of particles to use for exit-point analysis. "
             "Use with --particle-end to isolate one release group in daily-batch partoutput files.",
    )
    parser.add_argument(
        "--particle-end",
        type=int,
        default=None,
        metavar="N",
        help="0-based end index (inclusive) of particles to use for exit-point analysis.",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.grid_file):
        raise FileNotFoundError("grid file not found: {}".format(args.grid_file))

    out_file = args.out_file
    if out_file is None:
        base, _ = os.path.splitext(args.grid_file)
        out_file = base + "_footprints.nc"

    exit_csv = args.exit_csv
    if exit_csv is None:
        base, _ = os.path.splitext(out_file)
        exit_csv = base + "_domain_exit_points.csv"

    if os.path.exists(out_file) and not args.overwrite and args.skip_existing:
        print("SKIP output exists: {}".format(out_file))
        return 0

    partoutput_arg = args.partoutput
    if partoutput_arg is None:
        partoutput_arg = os.path.dirname(os.path.abspath(args.grid_file))
        print("No --partoutput supplied; defaulting to grid directory: {}".format(partoutput_arg))

    print("Reading grid file: {}".format(args.grid_file))
    ds = _open_dataset_auto(args.grid_file)
    try:
        var_name = _pick_sensitivity_var(ds)
        print("Using sensitivity variable: {}".format(var_name))
        
        # Validate that we can find spatial dimensions
        try:
            lat_dim, lon_dim = _find_spatial_dims(ds)
            print("Using lat/lon dimensions: {}/{}".format(lat_dim, lon_dim))
        except ValueError as e:
            print("WARNING: {}; output may have dimension naming issues".format(e))

        srr_timeint_2d, selected_height, integrated_height_levels = _compute_srr_timeint_2d(
            ds,
            var_name,
            footprint_outheight_m=args.footprint_outheight_m,
        )
        native_units = ds[var_name].attrs.get("units", "")
        srr_timeint_2d, conv_factor = _convert_to_m2s_per_mol(
            srr_timeint_2d,
            native_units,
            args.source_layer_thickness_m,
        )
        if selected_height is not None:
            print(
                "Using cumulative footprint integration up to {:.3f} m agl ({} levels)".format(
                    selected_height,
                    integrated_height_levels,
                )
            )
        else:
            print("No height dimension found in sensitivity variable; footprint computed from native field")
        release_time_value = _infer_release_time_value(ds)
        if release_time_value is None:
            print("WARNING: Could not infer release time from RELSTART or time coordinate; using 0.0")
            release_time_value = 0.0

        release_dt = _parse_release_datetime_from_grid_filename(args.grid_file)
        rel_lon, rel_lat = _infer_release_location(ds)
        release_meteo = {}
        if not args.disable_grib_meteo:
            grib_file = args.meteo_grib_file
            if grib_file is None:
                grib_file = _find_gf_file_for_time(args.meteo_grib_dir, release_dt)
            if grib_file is not None and os.path.isfile(grib_file):
                release_meteo = _extract_release_meteo_from_grib(grib_file, release_dt, rel_lon, rel_lat)
                if release_meteo:
                    print("Extracted release meteorology from GRIB: {}".format(grib_file))
                else:
                    print("WARNING: could not extract release meteorology from GRIB {}; keeping NaNs".format(grib_file))
            else:
                print("WARNING: no suitable GRIB file found for release meteorology; keeping NaNs")

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
            "description": "time-integrated footprint at receptor release time from selected outheight (converted to molar flux sensitivity)",
            "conversion_from_native_units": str(native_units),
            "conversion_factor_applied": float(conv_factor),
            "source_layer_thickness_m": float(args.source_layer_thickness_m),
            "molar_mass_air_kg_per_mol": float(MOLAR_MASS_AIR_KG_PER_MOL),
        })
        if selected_height is not None:
            out["srr"].attrs["footprint_outheight_m"] = float(selected_height)
            out["srr"].attrs["requested_footprint_outheight_m"] = float(args.footprint_outheight_m)
            out["srr"].attrs["footprint_height_integration"] = "surface_to_outheight_inclusive"
            out["srr"].attrs["footprint_height_levels_integrated"] = int(integrated_height_levels)

        # Optional domain-exit diagnostics from particle output.
        exits = None
        reused_exit_csv = False
        if args.use_existing_exit_csv and os.path.isfile(exit_csv):
            exits = _load_exit_points_csv(exit_csv)
            reused_exit_csv = True
            print("Using existing domain exit points CSV: {} ({} rows)".format(exit_csv, len(exits)))
        else:
            pds = _open_partoutput(partoutput_arg,
                                   particle_start=args.particle_start,
                                   particle_end=args.particle_end)
            try:
                if pds is None:
                    print("No partoutput files supplied/found; skipping domain-exit diagnostics.")
                else:
                    lon_var, lat_var, z_var = _find_particle_vars(pds)
                    if lon_var is None or lat_var is None:
                        print("Could not find particle longitude/latitude vars; skipping domain exits.")
                    else:
                        exits, _ = _derive_exit_points(pds, lon_var, lat_var, z_var=z_var)
                        print("Derived {} particle exit points".format(len(exits)))
            finally:
                if pds is not None:
                    pds.close()

        if exits is not None:
            lon_centers = np.asarray(ds["longitude"].values, dtype=float)
            lat_centers = np.asarray(ds["latitude"].values, dtype=float)
            _ = _build_exit_histogram(exits, lon_centers, lat_centers)

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

            if not reused_exit_csv:
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
            release_meteo=release_meteo,
        )
        if release_meteo:
            out.attrs["release_meteo_source"] = os.path.abspath(grib_file)

        encoding = _set_netcdf_compression(out, compression_level=4)
        out.to_netcdf(out_file, encoding=encoding)
        print("Wrote footprint products (compressed): {}".format(out_file))
    finally:
        ds.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
