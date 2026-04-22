#!/usr/bin/env python3
"""
Generic backward FLEXPART batch runner.

Loads domain configuration from domains_info/{DOMAIN}.txt
Loads receptor location from site_info.json
Generates FLEXPART input files and runs simulation.

Usage:
    ./run_backward_batch.py --domain EASTASIA --receptor GSN --num-particles 20000 --days 20 --end-time 2018022700 [options]
"""

import argparse
import json
import os
import sys
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re

# Configuration
REPO_ROOT = Path(__file__).parent.parent
SITE_INFO_FILE = REPO_ROOT / "site_domains" / "site_info.json"
DOMAINS_DIR = REPO_ROOT / "site_domains" / "domains_info"
OPTIONS_DIR = REPO_ROOT / "options"
FLEXPART_EXE_DEFAULT = REPO_ROOT / "src" / "FLEXPART"
FLEXPART_EXE_ETA = REPO_ROOT / "src" / "FLEXPART_ETA"
GFS_DATA_DIR = Path("/net/fs01/data/AGAGE/meteorology/gfs_grib")
GFS_AVAILABLE = GFS_DATA_DIR / "AVAILABLE"
POSTPROCESS_SCRIPT = REPO_ROOT / "run_scripts" / "postprocess_footprint.py"


def _gf_timestamp_from_name(filename):
    """Parse GFyymmddhh filename into a datetime, or return None if invalid."""
    match = re.match(r"^GF(\d{8})$", filename)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%y%m%d%H")
    except ValueError:
        return None


def generate_available_from_gf_files(gfs_data_dir, available_path, start_time=None, end_time=None):
    """
    Build an AVAILABLE file from GFyymmddhh files found in gfs_data_dir.

    If start/end are provided, prefer entries in that interval. If none are found,
    fall back to all discovered GF files.
    """
    entries = []
    for p in gfs_data_dir.iterdir():
        if not p.is_file():
            continue
        ts = _gf_timestamp_from_name(p.name)
        if ts is None:
            continue
        entries.append((ts, p.name))

    if not entries:
        raise FileNotFoundError(f"No GFyymmddhh files found in: {gfs_data_dir}")

    entries.sort(key=lambda x: x[0])

    selected = entries
    if start_time is not None and end_time is not None:
        # Extend the window by 6 h on each side so FLEXPART always has at least
        # one wind field before the simulation start and one after the end time.
        # Without this buffer the field bracketing the release/termination time
        # may be excluded, causing "NO METEO FIELDS AVAILABLE".
        buffer = timedelta(hours=6)
        in_window = [x for x in entries
                     if (start_time - buffer) <= x[0] <= (end_time + buffer)]
        if in_window:
            selected = in_window

    lines = [
        "XXXXXX EMPTY LINES XXXXXXXXX\n",
        "XXXXXX EMPTY LINES XXXXXXXX\n",
        "YYYYMMDD HHMMSS   name of the file(up to 80 characters)\n",
    ]

    for ts, fname in selected:
        lines.append(f"{ts.strftime('%Y%m%d')} {ts.strftime('%H')}0000      {fname:<12s} ON DISK\n")

    with open(available_path, "w") as f:
        f.writelines(lines)

    return len(selected)


def load_site_info(location_code):
    """
    Load location data from site_info.json.
    
    JSON structure:
        {
          "LOCATION": {
            "NETWORK": {
              "longitude": float,
              "latitude": float,
              "height_station_masl": int,
              "long_name": str,
              ...
            }
          }
        }
    
    Returns dict with keys: longitude, latitude, release_height_agl, long_name
    Raises ValueError if location not found.
    """
    with open(SITE_INFO_FILE) as f:
        data = json.load(f)
    
    if location_code not in data:
        raise ValueError(f"Location '{location_code}' not found in {SITE_INFO_FILE}")
    
    location_data = data[location_code]
    
    # Get first network (each location has at least one network with station info)
    station_info = None
    for network_name, network_data in location_data.items():
        if isinstance(network_data, dict) and 'longitude' in network_data:
            station_info = network_data
            break
    
    if not station_info:
        raise ValueError(f"No valid station data found for location '{location_code}'")
    
    # Prefer receptor inlet height above ground level (magl) from site metadata.
    release_height_agl = None
    for val in station_info.get('height', []):
        match = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*m\s*$", str(val), flags=re.IGNORECASE)
        if match:
            release_height_agl = float(match.group(1))
            break

    if release_height_agl is None:
        release_height_agl = 10.0

    return {
        'longitude': station_info.get('longitude'),
        'latitude': station_info.get('latitude'),
        'release_height_agl': release_height_agl,
        'height_station_masl': station_info.get('height_station_masl'),
        'height_name': station_info.get('height_name', []),
        'long_name': station_info.get('long_name', location_code),
    }


def load_domain_config(domain_name):
    """
    Load domain configuration from domains_info/{DOMAIN}.txt.
    
    Returns dict with keys: nX, nY, dX, dY, Xmin, Ymin
    Raises ValueError if file not found.
    """
    domain_file = DOMAINS_DIR / f"{domain_name}.txt"
    if not domain_file.exists():
        raise ValueError(f"Domain file not found: {domain_file}")
    
    config = {}
    with open(domain_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip("'\"")
                config[key] = value
    
    # Validate required keys
    required = ['Grid1_nX', 'Grid1_nY', 'Grid1_dX', 'Grid1_dY', 'Grid1_Xmin', 'Grid1_Ymin']
    for key in required:
        if key not in config:
            raise ValueError(f"Missing required key '{key}' in {domain_file}")
    
    return {
        'nX': int(config['Grid1_nX']),
        'nY': int(config['Grid1_nY']),
        'dX': float(config['Grid1_dX']),
        'dY': float(config['Grid1_dY']),
        'Xmin': float(config['Grid1_Xmin']),
        'Ymin': float(config['Grid1_Ymin']),
    }


def select_species_number(options_dir):
    """Pick an available species file, preferring valid tracer definitions."""
    species_dir = Path(options_dir) / 'SPECIES'
    if not species_dir.exists():
        raise FileNotFoundError(f"Missing SPECIES directory: {species_dir}")

    preferred = [24, 50]

    # 1) Preferred IDs that are already valid.
    for sid in preferred:
        p = species_dir / f"SPECIES_{sid:03d}"
        if p.exists() and _is_valid_species_namelist(p):
            return sid

    # 2) Any valid numeric species file.
    valid_candidates = []
    for p in species_dir.iterdir():
        m = re.match(r"^SPECIES_(\d{3})$", p.name)
        if m and p.is_file() and _is_valid_species_namelist(p):
            valid_candidates.append(int(m.group(1)))
    if valid_candidates:
        return sorted(valid_candidates)[0]

    # 3) If none are valid, return an existing preferred ID so caller can repair it.
    for sid in preferred:
        if (species_dir / f"SPECIES_{sid:03d}").exists():
            return sid

    # 4) Last resort: any numeric species file (to be repaired by caller).
    raw_candidates = []
    for p in species_dir.iterdir():
        m = re.match(r"^SPECIES_(\d{3})$", p.name)
        if m and p.is_file():
            raw_candidates.append(int(m.group(1)))

    if not raw_candidates:
        raise FileNotFoundError(
            f"No numeric species files found in {species_dir} (expected SPECIES_###)"
        )

    return sorted(raw_candidates)[0]


def _is_valid_species_namelist(species_file):
    """Return True if species file looks like a valid namelist with closing '/'."""
    try:
        lines = Path(species_file).read_text().splitlines()
    except OSError:
        return False

    has_header = any(line.strip().upper().startswith('&SPECIES_PARAMS') for line in lines)
    has_close = any(line.strip() == '/' for line in lines)
    return has_header and has_close


def ensure_valid_species_file(options_dir, species_num):
    """Ensure selected numeric species file exists and is valid namelist format."""
    species_dir = Path(options_dir) / 'SPECIES'
    species_file = species_dir / f"SPECIES_{species_num:03d}"

    if species_file.exists() and _is_valid_species_namelist(species_file):
        return

    fallback_content = (
        "&SPECIES_PARAMS\n"
        " PSPECIES=\"AIRTRACER\",\n"
        " PWEIGHTMOLAR=29.0,\n"
        " /\n"
    )
    species_file.write_text(fallback_content)
    print(f"Wrote fallback species file: {species_file}")


def generate_releases_file(location_data, num_particles, start_time, end_time, outfile, species_num):
    """Generate RELEASES file for backward simulation with instantaneous release."""
    lon = location_data['longitude']
    lat = location_data['latitude']
    release_height_agl = location_data['release_height_agl']

    # Instantaneous release at end_time (release time).
    idate1 = end_time.strftime('%Y%m%d')
    itime1 = end_time.strftime('%H%M%S')
    idate2 = end_time.strftime('%Y%m%d')
    itime2 = end_time.strftime('%H%M%S')

    content = f"""&RELEASES_CTRL
 NSPEC=1,
 SPECNUM_REL={species_num},
 /
&RELEASE
 IDATE1={idate1},ITIME1={itime1},
 IDATE2={idate2},ITIME2={itime2},
 LON1={lon},LON2={lon},
 LAT1={lat},LAT2={lat},
 Z1={release_height_agl},Z2={release_height_agl},
 ZKIND=1,
 PARTS={num_particles},
 MASS=1.0,
 COMMENT='{location_data['long_name']}'
/
"""

    with open(outfile, 'w') as f:
        f.write(content)
    print(f"Generated RELEASES: {outfile} (SPECNUM_REL={species_num})")


def generate_outgrid_file(domain_config, outfile):
    """
    Generate OUTGRID file based on domain configuration.

    GFS data is normalised to -180/180 by FLEXPART internally (NCEP convention),
    so the output grid right edge must not exceed 180°.  If the domain definition
    crosses the dateline we clip the grid and warn.
    """
    config = domain_config
    nX = config['nX']
    xmin = config['Xmin']
    dx = config['dX']

    xmax = xmin + nX * dx
    if xmax > 180.0:
        nX_clipped = int((180.0 - xmin) / dx)
        xmax_clipped = xmin + nX_clipped * dx
        print(
            f"WARNING: domain right edge {xmax:.3f}° exceeds GFS model domain (180°). "
            f"Clipping NUMXGRID from {nX} to {nX_clipped} (xmax={xmax_clipped:.3f}°)."
        )
        nX = nX_clipped

    content = f"""&OUTGRID
 OUTLON0=    {xmin:.3f},
 OUTLAT0=    {config['Ymin']:.3f},
 NUMXGRID=   {nX},
 NUMYGRID=   {config['nY']},
 DXOUT=      {dx:.3f},
 DYOUT=      {config['dY']:.3f},
 OUTHEIGHTS= 100.0, 500.0, 1000.0, 50000.0,
 /
"""

    with open(outfile, 'w') as f:
        f.write(content)
    print(f"Generated OUTGRID: {outfile}")


def write_pathnames_file(options_dir, output_dir, meteo_dir, available_file, pathnames_file):
    """Write FLEXPART pathnames file expected in the current working directory."""
    lines = [
        f"{Path(options_dir).resolve()}/\n",
        f"{Path(output_dir).resolve()}/\n",
        f"{Path(meteo_dir).resolve()}/\n",
        f"{Path(available_file).resolve()}\n",
    ]
    with open(pathnames_file, 'w') as f:
        f.writelines(lines)
    print(f"Generated pathnames: {pathnames_file}")


def _set_command_value(lines, key, value):
    """Set a single namelist key value in COMMAND lines, preserving comments."""
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^,]+)(,.*)$", re.IGNORECASE)
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if m:
            lines[i] = f"{m.group(1)}{value}{m.group(3)}\n"
            return True
    return False


def _get_command_int(lines, key):
    """Get integer namelist value from COMMAND lines, or None if missing/unparseable."""
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^,]+),", re.IGNORECASE)
    for line in lines:
        m = pattern.match(line)
        if m:
            try:
                return int(float(m.group(1).strip()))
            except ValueError:
                return None
    return None


def update_command_file(
    command_file,
    start_time,
    end_time,
    domain_config,
    nxshift_override=None,
    ipout_override=None,
    lsubgrid_value=1,
):
    """Update COMMAND for backward sensitivity runs and consistent sampling settings."""
    with open(command_file) as f:
        lines = f.readlines()

    domain_xmax = domain_config['Xmin'] + domain_config['nX'] * domain_config['dX']
    # GFS GRIB files start at 0°E and cover 0→360° natively, so NXSHIFT=0 is
    # always correct for GFS regardless of whether the domain crosses 180°E.
    # (ECMWF default would be 359; this script targets GFS only.)
    if nxshift_override is None:
        nxshift_value = "0"
    else:
        nxshift_value = str(nxshift_override)

    updates = [
        ("LDIRECT", "-1"),
        ("IOUT", "9"),
        ("IBDATE", start_time.strftime('%Y%m%d')),
        ("IBTIME", start_time.strftime('%H%M%S')),
        ("IEDATE", end_time.strftime('%Y%m%d')),
        ("IETIME", end_time.strftime('%H%M%S')),
        ("NXSHIFT", nxshift_value),
        # FLEXINVERT-style default: enable subgrid terrain parameterization.
        ("LSUBGRID", str(int(lsubgrid_value))),
        # Use receptor mixing-ratio units but keep initial-condition mode off.
        ("IND_RECEPTOR", "2"),
        ("LINIT_COND", "0"),
        # Surface-layer only output (0–100 m agl) to match flexinvertplus footprints.
        ("SFC_ONLY", "1"),
        # FLEXINVERT-style output layout: one inversion-formatted output per release.
        ("LINVERSIONOUT", "1"),
    ]

    # Also try legacy key; no-op if absent.
    _set_command_value(lines, "SURF_ONLY", "1")

    # Default to IPOUT=2 so partoutput_*.nc is available for domain-exit diagnostics.
    if ipout_override is None:
        ipout_override = 2
    updates.append(("IPOUT", str(ipout_override)))

    for key, value in updates:
        _set_command_value(lines, key, value)

    lsynctime = _get_command_int(lines, "LSYNCTIME")
    lrecoutsample = _get_command_int(lines, "LRECOUTSAMPLE")
    if lsynctime and lrecoutsample and (lrecoutsample % lsynctime != 0):
        _set_command_value(lines, "LRECOUTSAMPLE", str(lsynctime))
        print(
            f"Adjusted LRECOUTSAMPLE to {lsynctime} to be a multiple of LSYNCTIME={lsynctime}"
        )

    with open(command_file, 'w') as f:
        f.writelines(lines)

    print(f"Set NXSHIFT={nxshift_value} (domain xmax={domain_xmax:.3f})")
    print("Set IOUT=9 for backward gridded sensitivity output (grid_time_*.nc)")
    print(f"Set IPOUT={ipout_override}")
    print("Set IND_RECEPTOR=2, LINIT_COND=0, SFC_ONLY/SURF_ONLY=1 (surface-layer only)")
    print("Set LINVERSIONOUT=1 (FLEXINVERT-style inversion output format)")
    print(f"Set LSUBGRID={int(lsubgrid_value)}")
    print(f"Updated COMMAND: {command_file}")


def run_flexpart(options_dir, gfs_data_dir, gfs_available_file, flexpart_exe):
    """
    Run FLEXPART executable from options_dir.
    """
    if not flexpart_exe.exists():
        raise FileNotFoundError(f"FLEXPART executable not found: {flexpart_exe}")
    
    cmd = [str(flexpart_exe)]
    
    # Verify required options exist
    required_files = ['COMMAND', 'RELEASES', 'OUTGRID', 'OUTGRID_NEST']
    for fname in required_files:
        fpath = options_dir / fname
        if not fpath.exists():
            raise FileNotFoundError(f"Missing required options file: {fpath}")
    
    # Verify GFS data and AVAILABLE file
    if not gfs_data_dir.exists():
        raise FileNotFoundError(f"GFS data directory not found: {gfs_data_dir}")
    if not gfs_available_file.exists():
        raise FileNotFoundError(f"AVAILABLE file not found: {gfs_available_file}")
    
    print(f"\nRunning FLEXPART from: {options_dir}")
    print(f"Executable: {flexpart_exe}")
    print(f"GFS data: {gfs_data_dir}")
    
    # Change to options directory and run
    original_cwd = os.getcwd()
    try:
        os.chdir(options_dir)
        result = subprocess.run(cmd)
        return result.returncode
    finally:
        os.chdir(original_cwd)


def run_postprocess(
    output_dir,
    receptor,
    domain,
    release_height_agl,
    end_time,
    lowest_magl=100.0,
    source_layer_thickness_m=100.0,
    postprocess_python=None,
):
    """Run automatic postprocessing to create an AGAGE-like footprint NetCDF."""
    if not POSTPROCESS_SCRIPT.exists():
        print(f"WARNING: postprocess script not found: {POSTPROCESS_SCRIPT}")
        return None

    grid_files = sorted(output_dir.glob("grid_time_*.nc"))
    if not grid_files:
        print(f"WARNING: no grid_time_*.nc files found in {output_dir}; skipping postprocess")
        return None

    exact_name = f"grid_time_{end_time.strftime('%Y%m%d%H%M%S')}.nc"
    exact_file = output_dir / exact_name
    if exact_file.exists():
        grid_file = exact_file
    else:
        grid_file = None
        best_delta = None
        for gf in grid_files:
            m = re.match(r"^grid_time_(\d{14})\.nc$", gf.name)
            if not m:
                continue
            try:
                ts = datetime.strptime(m.group(1), "%Y%m%d%H%M%S")
            except ValueError:
                continue
            delta = abs((ts - end_time).total_seconds())
            if best_delta is None or delta < best_delta:
                best_delta = delta
                grid_file = gf

        if grid_file is None:
            grid_file = grid_files[0]

        print(
            f"WARNING: exact release-time grid file '{exact_name}' not found; using nearest: {grid_file.name}"
        )

    # Keep receptor release height readable in filenames (e.g. 10 -> "10", 10.5 -> "10.5").
    try:
        h = float(release_height_agl)
        if h.is_integer():
            magl_label = str(int(h))
        else:
            magl_label = f"{h:.3f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        magl_label = str(release_height_agl)
    out_name = f"{receptor}-{magl_label}magl_FLEXPART_GFS_{domain}_inert_{end_time.strftime('%Y%m%d%H')}.nc"
    out_file = output_dir / out_name
    exit_csv = output_dir / out_name.replace(".nc", "_domain_exit_points.csv")

    postprocess_python = str(postprocess_python or sys.executable)

    cmd = [
        postprocess_python,
        "-X",
        "faulthandler",
        str(POSTPROCESS_SCRIPT),
        "--grid-file",
        str(grid_file),
        "--out-file",
        str(out_file),
        "--source-layer-thickness-m",
        str(source_layer_thickness_m),
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
        "GFS",
    ]

    print("\nRunning automatic postprocessing:")
    print(f"  Grid file: {grid_file}")
    print(f"  Output file: {out_file}")
    print(f"  Python: {postprocess_python}")

    # Some shared filesystems/HDF5 combinations are unstable with file locking.
    # Disable locking for this short-lived read/write postprocess subprocess.
    post_env = os.environ.copy()
    post_env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    max_attempts = 2
    for attempt in range(1, max_attempts + 1):
        result = subprocess.run(cmd, env=post_env)
        if result.returncode == 0:
            print(f"Postprocessed footprint file: {out_file}")
            return out_file

        print(f"WARNING: postprocess attempt {attempt}/{max_attempts} failed with exit code {result.returncode}")
        if attempt < max_attempts:
            print("Retrying postprocess once...")

    print("WARNING: automatic postprocessing failed after retry; consider rerunning postprocess as a separate job.")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Run FLEXPART backward simulation with domain and location parameters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic run
  ./run_backward_batch.py --domain EASTASIA --receptor ADR

  # Custom particle count and duration
  ./run_backward_batch.py --domain EASTASIA --receptor ADR \\
      --num-particles 50000 --days 30

  # Specifiy start time
  ./run_backward_batch.py --domain EASTASIA --receptor ADR \\
      --end-time 2026040100

  # Dry-run (show what would be done)
  ./run_backward_batch.py --domain EASTASIA --receptor ADR --dry-run
""",
    )
    
    parser.add_argument(
        '--domain',
        required=True,
        metavar='NAME',
        help='Domain name (e.g., EASTASIA). Config loaded from domains_info/{NAME}.txt',
    )
    parser.add_argument(
        '--receptor',
        required=True,
        metavar='CODE',
        help='Receptor code (e.g., ADR). Loaded from site_info.json',
    )
    parser.add_argument(
        '--num-particles',
        type=int,
        default=20000,
        metavar='N',
        help='Number of particles (default: %(default)s)',
    )
    parser.add_argument(
        '--days',
        type=int,
        default=20,
        metavar='N',
        help='Backward simulation duration in days (default: %(default)s)',
    )
    parser.add_argument(
        '--end-time',
        metavar='YYYYMMDDHH',
        help='End time (release time) in UTC (default: now)',
    )
    parser.add_argument(
        '--outdir',
        metavar='DIR',
        help='Output directory for FLEXPART run (default: ./backward_{DOMAIN}_{LOCATION}_{TIME})',
    )
    parser.add_argument(
        '--gfs-data',
        type=Path,
        default=GFS_DATA_DIR,
        metavar='DIR',
        help=f'GFS meteorological data directory (default: {GFS_DATA_DIR})',
    )
    parser.add_argument(
        '--gfs-available',
        type=Path,
        default=GFS_AVAILABLE,
        metavar='FILE',
        help=f'AVAILABLE file for GFS data (default: {GFS_AVAILABLE})',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing',
    )
    parser.add_argument(
        '--nxshift',
        type=int,
        default=None,
        metavar='N',
        help='Override NXSHIFT in COMMAND. If omitted, auto-select based on domain extent.',
    )
    parser.add_argument(
        '--ipout',
        type=int,
        choices=[0, 1, 2],
        default=None,
        metavar='N',
        help='Optional IPOUT override in COMMAND (default when omitted: 2=emit partoutput for exit diagnostics).',
    )
    parser.add_argument(
        '--lsubgrid',
        type=int,
        choices=[0, 1],
        default=1,
        metavar='N',
        help='Set LSUBGRID in COMMAND (default: %(default)s).',
    )
    parser.add_argument(
        '--no-postprocess',
        action='store_true',
        help='Disable automatic postprocessing of grid_time output into AGAGE-style footprint NetCDF.',
    )
    parser.add_argument(
        '--postprocess-lowest-magl',
        type=float,
        default=100.0,
        metavar='M',
        help='Low-level integration top (m agl) used in automatic postprocessing (default: %(default)s).',
    )
    parser.add_argument(
        '--postprocess-source-layer-thickness-m',
        type=float,
        default=100.0,
        metavar='M',
        help='Source-layer thickness (m) used to convert SRR to m2 s mol-1 in postprocessing (default: %(default)s).',
    )
    parser.add_argument(
        '--postprocess-python',
        type=Path,
        default=None,
        metavar='FILE',
        help='Python executable used for automatic postprocessing. Default: reuse the current interpreter.',
    )
    parser.add_argument(
        '--executable',
        type=Path,
        default=None,
        metavar='FILE',
        help='FLEXPART executable path. Default prefers src/FLEXPART (non-ETA), then falls back to src/FLEXPART_ETA.',
    )
    
    args = parser.parse_args()

    if args.executable is not None:
        flexpart_exe = args.executable
    elif FLEXPART_EXE_DEFAULT.exists():
        flexpart_exe = FLEXPART_EXE_DEFAULT
    else:
        flexpart_exe = FLEXPART_EXE_ETA
    
    # Load configurations
    print(f"Loading receptor: {args.receptor}")
    location = load_site_info(args.receptor)
    print(
        f"  {location['long_name']}: {location['latitude']:.2f}°N, {location['longitude']:.2f}°E, "
        f"release={location['release_height_agl']} m agl"
    )
    
    print(f"\nLoading domain: {args.domain}")
    domain = load_domain_config(args.domain)
    print(f"  Grid: {domain['nX']}x{domain['nY']}, cell: {domain['dX']:.3f}°x{domain['dY']:.3f}°")
    
    # Calculate times (backward)
    end_time = datetime.strptime(args.end_time, '%Y%m%d%H') if args.end_time else datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0, tzinfo=None)
    start_time = end_time - timedelta(days=args.days)
    
    print(f"\nSimulation window:")
    print(f"  Start (begin backward): {start_time.strftime('%Y-%m-%d %H:%M UTC')}")
    print(f"  End (release time):     {end_time.strftime('%Y-%m-%d %H:%M UTC')}")
    print(f"  Duration: {args.days} days backward")
    print("  Release mode: instantaneous at end time")
    print(f"  Particles: {args.num_particles}")
    
    # Create output directory
    if args.outdir:
        out_dir = Path(args.outdir)
    else:
        timestamp = end_time.strftime('%Y%m%d_%H%M')
        out_dir = Path(f"./backward_{args.domain}_{args.receptor}_{timestamp}")
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy full options template (files and directories), then override key files.
    template_dir = OPTIONS_DIR
    print(f"\nPreparing options in: {out_dir}")
    import shutil
    for src in template_dir.iterdir():
        dst = out_dir / src.name
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
        print(f"  Copied {src.name}")

    update_command_file(
        out_dir / 'COMMAND',
        start_time,
        end_time,
        domain,
        nxshift_override=args.nxshift,
        ipout_override=args.ipout,
        lsubgrid_value=args.lsubgrid,
    )
    
    # Generate RELEASES and OUTGRID
    species_num = select_species_number(out_dir)
    ensure_valid_species_file(out_dir, species_num)
    generate_releases_file(
        location,
        args.num_particles,
        start_time,
        end_time,
        out_dir / 'RELEASES',
        species_num,
    )
    generate_outgrid_file(domain, out_dir / 'OUTGRID')
    
    # Show what would happen
    print(f"\nReady to run FLEXPART:")
    print(f"  Options: {out_dir}")
    print(f"  GFS data: {args.gfs_data}")

    effective_available = args.gfs_available
    if not effective_available.exists():
        auto_available = out_dir / 'AVAILABLE'
        print(f"  Requested AVAILABLE not found: {effective_available}")
        print(f"  Auto-generating AVAILABLE from GF files in {args.gfs_data}")
        n_entries = generate_available_from_gf_files(
            args.gfs_data,
            auto_available,
            start_time=start_time,
            end_time=end_time,
        )
        effective_available = auto_available
        print(f"  Generated AVAILABLE: {effective_available} ({n_entries} entries)")
    else:
        print(f"  AVAILABLE: {effective_available}")

    # FLEXPART requires file "pathnames" in the current working directory.
    output_dir = out_dir / 'output'
    if output_dir.exists():
        # Clean stale output to prevent conflicts with previous runs.
        for f in output_dir.iterdir():
            if f.is_file():
                f.unlink()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_pathnames_file(
        options_dir=out_dir,
        output_dir=output_dir,
        meteo_dir=args.gfs_data,
        available_file=effective_available,
        pathnames_file=out_dir / 'pathnames',
    )
    
    if args.dry_run:
        print("\n(Dry-run mode: not executing)")
        return 0
    
    # Run FLEXPART
    print("\n" + "="*60)
    returncode = run_flexpart(out_dir, args.gfs_data, effective_available, flexpart_exe)
    print("="*60)
    
    if returncode == 0:
        print(f"\n✓ FLEXPART completed successfully")
        print(f"Output directory: {out_dir}")
        # Look for output files
        output_dir = out_dir / 'output'
        if output_dir.exists():
            output_files = list(output_dir.glob('*'))
            print(f"Output files: {len(output_files)} files generated")

            if not args.no_postprocess:
                run_postprocess(
                    output_dir=output_dir,
                    receptor=args.receptor,
                    domain=args.domain,
                    release_height_agl=location['release_height_agl'],
                    end_time=end_time,
                    lowest_magl=args.postprocess_lowest_magl,
                    source_layer_thickness_m=args.postprocess_source_layer_thickness_m,
                    postprocess_python=args.postprocess_python,
                )
            else:
                print("Automatic postprocessing disabled (--no-postprocess).")
    else:
        print(f"\n✗ FLEXPART failed with exit code {returncode}")
        if returncode == -4:
            print(
                "  Hint: exit -4 = SIGILL (illegal instruction). "
                "This usually means the binary was compiled with a CPU-specific "
                "target (for example -march=native) that is unsupported on this node."
            )
            print(
                "  Rebuild with a portable target, e.g. make ... arch=x86-64, "
                "then retry."
            )
        elif returncode == -11:
            print(
                "  Hint: exit -11 = SIGSEGV. Check the corresponding Slurm .err file "
                "for the Fortran backtrace to identify the crashing routine."
            )
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
