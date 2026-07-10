#!/usr/bin/env python3
"""
Daily batch backward FLEXPART runner.

Runs all 24 hourly releases for a single calendar day in one FLEXPART simulation.
FLEXPART writes one grid_time_*.nc per release group, so 24 files are produced.

Usage:
    ./run_backward_batch_daily.py --domain WESTUSA --receptor THD --date 20190101 [options]

Met-field I/O advantage:
    The 20-day backward met window is read only once and shared across all 24 releases,
    giving roughly a 20x reduction in I/O compared to 24 separate hourly runs.

Exit-point note:
    IPOUT defaults to 0 for daily runs.  If IPOUT > 0, all 24 releases' particles
    share the same partoutput_*.nc files; exit-point attribution in postprocessing
    is then combined across all releases for each hourly footprint.  For per-release
    exit diagnostics use hourly mode (run_backward_batch.py).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

# Locate repo root and import shared utilities from run_scripts/
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT  = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_ROOT / "run_scripts"))

from run_backward_batch import (
    GFS_DATA_DIR,
    GFS_AVAILABLE,
    OPTIONS_DIR,
    FLEXPART_EXE_DEFAULT,
    FLEXPART_EXE_ETA,
    generate_available_from_gf_files,
    load_site_info,
    load_domain_config,
    select_species_number,
    ensure_valid_species_file,
    generate_outgrid_file,
    write_pathnames_file,
    update_command_file,
    run_flexpart,
)


POSTPROCESS_DAILY_SCRIPT = SCRIPT_DIR / "postprocess_daily_footprints.py"


def generate_releases_file_daily(location_data, num_particles, date, outfile, species_num):
    """Generate RELEASES with 24 instantaneous hourly release groups for one day.

    FLEXPART (backward mode) produces one grid_time_*.nc per release group, so
    this generates 24 output files from a single simulation.
    """
    lon = location_data['longitude']
    lat = location_data['latitude']
    z   = location_data['release_height_agl']

    header = (
        "&RELEASES_CTRL\n"
        " NSPEC=1,\n"
        f" SPECNUM_REL={species_num},\n"
        " /\n"
    )

    blocks = [header]
    for hour in range(24):
        rel_dt = date + timedelta(hours=hour)
        idate  = rel_dt.strftime('%Y%m%d')
        itime  = rel_dt.strftime('%H%M%S')
        blocks.append(
            f"&RELEASE\n"
            f" IDATE1={idate},ITIME1={itime},\n"
            f" IDATE2={idate},ITIME2={itime},\n"
            f" LON1={lon},LON2={lon},\n"
            f" LAT1={lat},LAT2={lat},\n"
            f" Z1={z},Z2={z},\n"
            f" ZKIND=1,\n"
            f" PARTS={num_particles},\n"
            f" MASS=1.0,\n"
            f" COMMENT='{location_data['long_name']} {itime[:2]}:00'\n"
            f"/\n"
        )

    with open(outfile, 'w') as f:
        f.writelines(blocks)
    print(f"Generated RELEASES: {outfile} (24 hourly groups, SPECNUM_REL={species_num})")


def run_postprocess_all(
    output_dir,
    receptor,
    domain,
    release_height_agl,
    num_particles,
    gfs_data_dir,
    footprint_outheight_m,
    source_layer_thickness_m,
    postprocess_python=None,
):
    """Run daily postprocessing once, deriving particle exits once per day."""
    grid_files = sorted(output_dir.glob("grid_time_*.nc"))
    if not grid_files:
        print("WARNING: no grid_time_*.nc found in output; skipping postprocessing.")
        return []

    print(f"\nPostprocessing {len(grid_files)} grid_time file(s) with daily fast path...")

    try:
        h = float(release_height_agl)
        magl_label = str(int(h)) if h.is_integer() else f"{h:.3f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        magl_label = str(release_height_agl)

    python_exe = str(postprocess_python or sys.executable)
    post_env = os.environ.copy()
    post_env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    cmd = [
        python_exe,
        "-X",
        "faulthandler",
        str(POSTPROCESS_DAILY_SCRIPT),
        "--run-dir",
        str(output_dir.parent),
        "--output-dir",
        str(output_dir),
        "--source-layer-thickness-m",
        str(source_layer_thickness_m),
        "--footprint-outheight-m",
        str(footprint_outheight_m),
        "--partoutput",
        str(output_dir),
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
        "--meteo-grib-dir",
        str(gfs_data_dir),
    ]

    print("  RUN " + " ".join(str(part) for part in cmd))
    result = subprocess.run(cmd, env=post_env)
    if result.returncode != 0:
        print(f"  WARNING: daily postprocessing failed (exit {result.returncode})")
        return []

    out_files = sorted(output_dir.glob(f"*-{magl_label}magl_FLEXPART_CFSv2_{domain}_inert_*.nc"))
    print(f"Postprocessed {len(out_files)} daily footprint file(s) successfully.")
    return out_files


def main():
    parser = argparse.ArgumentParser(
        description="Run 24 hourly FLEXPART backward releases for one calendar day.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ./run_backward_batch_daily.py --domain WESTUSA --receptor THD --date 20190101
  ./run_backward_batch_daily.py --domain WESTUSA --receptor THD --date 20190101 --dry-run
  ./run_backward_batch_daily.py --domain WESTUSA --receptor THD --date 20190101 --no-postprocess
""",
    )
    parser.add_argument('--domain',        required=True, metavar='NAME',
                        help='Domain name (e.g., WESTUSA). Config from domains_info/{NAME}.txt')
    parser.add_argument('--receptor',      required=True, metavar='CODE',
                        help='Receptor code (e.g., THD). Loaded from site_info.json')
    parser.add_argument('--date',          required=True, metavar='YYYYMMDD',
                        help='Calendar day: all 24 hourly releases are generated for this day')
    parser.add_argument('--days',          type=int, default=20, metavar='N',
                        help='Backward simulation duration in days (default: %(default)s)')
    parser.add_argument('--num-particles', type=int, default=20000, metavar='N',
                        help='Particles per release group (default: %(default)s); '
                             '24 groups are run, so total particles = N * 24')
    parser.add_argument('--outdir',        metavar='DIR', default=None,
                        help='Output directory (default: ./backward_{DOMAIN}_{RECEPTOR}_{DATE})')
    parser.add_argument('--gfs-data',      type=Path, default=GFS_DATA_DIR, metavar='DIR')
    parser.add_argument('--gfs-available', type=Path, default=GFS_AVAILABLE, metavar='FILE')
    parser.add_argument('--dry-run',       action='store_true',
                        help='Print what would be done without executing')
    parser.add_argument('--nxshift',       type=int, default=None, metavar='N',
                        help='Override NXSHIFT; auto-selected if omitted')
    parser.add_argument('--ipout',         type=int, choices=[0, 1, 2], default=None, metavar='N',
                        help='IPOUT override (default: 1). All 24 releases share partoutput '
                             'files, but per-release exit attribution is achieved by slicing '
                             'particles by sequential index (h*N to (h+1)*N-1 for hour h).')
    parser.add_argument('--lsubgrid',      type=int, choices=[0, 1], default=1, metavar='N')
    parser.add_argument('--linit-cond',    type=int, choices=[0, 1, 2], default=1, metavar='N')
    parser.add_argument('--no-postprocess', action='store_true',
                        help='Skip automatic postprocessing after FLEXPART')
    parser.add_argument('--postprocess-footprint-outheight-m',
                        type=float, default=100.0, metavar='M')
    parser.add_argument('--postprocess-source-layer-thickness-m',
                        type=float, default=100.0, metavar='M')
    parser.add_argument('--postprocess-python', type=Path, default=None, metavar='FILE')
    parser.add_argument('--executable',    type=Path, default=None, metavar='FILE',
                        help='FLEXPART executable path')

    args = parser.parse_args()

    # Parse date
    try:
        date = datetime.strptime(args.date, '%Y%m%d')
    except ValueError:
        print(f"ERROR: --date must be YYYYMMDD, got: {args.date!r}")
        return 1

    # IPOUT: default 1 for daily — partoutput needed for exit-point analysis.
    # Per-release particle slicing is used in run_postprocess_all to isolate
    # each release group's particles (indices h*N to (h+1)*N-1).
    ipout = args.ipout if args.ipout is not None else 1

    # Simulation window:
    #   IBDATE/IBTIME = earliest backward extent (sim_start, 20 days before first release)
    #   IEDATE/IETIME = last release time (23:00 on the date)
    sim_start = date - timedelta(days=args.days)
    sim_end   = date + timedelta(hours=23)

    # Choose executable
    if args.executable is not None:
        flexpart_exe = args.executable
    elif FLEXPART_EXE_DEFAULT.exists():
        flexpart_exe = FLEXPART_EXE_DEFAULT
    else:
        flexpart_exe = FLEXPART_EXE_ETA

    # Load site/domain
    print(f"Loading receptor: {args.receptor}")
    location = load_site_info(args.receptor)
    print(
        f"  {location['long_name']}: {location['latitude']:.2f}°N, "
        f"{location['longitude']:.2f}°E, release={location['release_height_agl']} m agl"
    )

    print(f"\nLoading domain: {args.domain}")
    domain = load_domain_config(args.domain)
    print(f"  Grid: {domain['nX']}x{domain['nY']}, cell: {domain['dX']:.3f}°x{domain['dY']:.3f}°")

    print(f"\nDaily simulation: {args.date}")
    print(f"  24 releases: {date.strftime('%Y-%m-%d')} 00:00 – 23:00 UTC")
    print(f"  Backward window: {sim_start.strftime('%Y-%m-%d %H:%M')} -> {sim_end.strftime('%Y-%m-%d %H:%M')} UTC")
    print(f"  Particles per release: {args.num_particles}  (total in run: {args.num_particles * 24:,})")
    print(f"  IPOUT: {ipout}")

    # Output directory
    out_dir = Path(args.outdir) if args.outdir else Path(f"./backward_{args.domain}_{args.receptor}_{args.date}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Copy options template
    print(f"\nPreparing options in: {out_dir}")
    for src in OPTIONS_DIR.iterdir():
        dst = out_dir / src.name
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
        print(f"  Copied {src.name}")

    # COMMAND: backward from sim_start to last release (sim_end = date 23:00)
    update_command_file(
        out_dir / 'COMMAND',
        sim_start,
        sim_end,
        domain,
        nxshift_override=args.nxshift,
        ipout_override=ipout,
        lsubgrid_value=args.lsubgrid,
        linit_cond_value=args.linit_cond,
    )

    # RELEASES: 24 hourly groups
    species_num = select_species_number(out_dir)
    ensure_valid_species_file(out_dir, species_num)
    generate_releases_file_daily(
        location,
        args.num_particles,
        date,
        out_dir / 'RELEASES',
        species_num,
    )
    generate_outgrid_file(domain, out_dir / 'OUTGRID')

    # AVAILABLE file
    effective_available = args.gfs_available
    if not effective_available.exists():
        auto_available = out_dir / 'AVAILABLE'
        print(f"  Requested AVAILABLE not found: {effective_available}")
        print(f"  Auto-generating from GF files in {args.gfs_data}")
        n_entries = generate_available_from_gf_files(
            args.gfs_data,
            auto_available,
            start_time=sim_start,
            end_time=sim_end,
        )
        effective_available = auto_available
        print(f"  Generated AVAILABLE: {effective_available} ({n_entries} entries)")
    else:
        print(f"  AVAILABLE: {effective_available}")

    # Output subdirectory
    output_dir = out_dir / 'output'
    if output_dir.exists():
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
    print("\n" + "=" * 60)
    returncode = run_flexpart(out_dir, args.gfs_data, effective_available, flexpart_exe)
    print("=" * 60)

    if returncode != 0:
        print(f"\n✗ FLEXPART failed with exit code {returncode}")
        return 1

    print(f"\n✓ FLEXPART completed successfully")
    n_out = len(list(output_dir.glob('*')))
    print(f"Output: {n_out} file(s) in {output_dir}")

    if not args.no_postprocess:
        run_postprocess_all(
            output_dir=output_dir,
            receptor=args.receptor,
            domain=args.domain,
            release_height_agl=location['release_height_agl'],
            num_particles=args.num_particles,
            gfs_data_dir=args.gfs_data,
            footprint_outheight_m=args.postprocess_footprint_outheight_m,
            source_layer_thickness_m=args.postprocess_source_layer_thickness_m,
            postprocess_python=args.postprocess_python,
        )
    else:
        print("Automatic postprocessing disabled (--no-postprocess).")

    return 0


if __name__ == '__main__':
    sys.exit(main())
