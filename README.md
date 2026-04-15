![](https://www.flexpart.eu/chrome/site/flexpart_banner.png)
# Welcome to Flexpart - The Lagrangian particle dispersion model

> **Note:** This is a **forked and adapted version** of FLEXPART optimized for backward simulations of inert particles using GFS meteorological data. Currently configured primarily for use on the **MIT Svante HPC system**. The codebase includes tools for seamless GFS data acquisition and batch processing of backward trajectories.

The the main development site is @ University of Vienna.

Other references:

- [FLEXPART.eu](https://flexpart.eu)
- [FLEXPART@NILU](https://git.nilu.no/flexpart/flexpart)

### MIT/Svante install instructions

Build and install on an EDR compute node on Svante (not on login or fsXX nodes).

Start an interactive EDR shell:
`srun -n 1 -p edr --pty /bin/bash`

Create and activate a conda environment:
`conda create --name flexpart`
`conda activate flexpart`

Install required packages:
`conda install -c conda-forge eccodes gfortran_linux-64 netcdf-fortran`

Clone this repository (URL omitted here), then build:
`cd flexpart_gfs/src`
`make -f makefile_svante clean`
`make -f makefile_svante eta=no -j4`


### Running Slurm Array Jobs (Svante/HPC)

For multi-date production runs, use the Slurm array launcher script:

```bash
./run_scripts/run_slurm_array_backward.sh
```

How to use it:
- Edit configuration values in `run_scripts/slurm_array_config.sh`.
- Set at least: `START_DATE`, `END_DATE`, `DOMAIN`, `RECEPTOR`.
- Run the script from the repository root.

Default output root is:
- `/net/fs06/d2/$USER/flexpart_outs`

After each task completes, the workflow keeps only the final postprocessed footprint file for that release hour:
- `*_FLEXPART_GFS_<DOMAIN>_inert_<YYYYMMDDHH>.nc`

All other files in that task output folder are deleted automatically.

Advanced usage:
- Direct submit helper (date range to array spec):

```bash
./run_scripts/submit_slurm_array_backward.sh 2018020100 2018022800 20
```

- Full override example:

```bash
DOMAIN=EASTASIA \
RECEPTOR=GSN \
STEP_HOURS=1 \
BACKWARD_DAYS=20 \
NUM_PARTICLES=20000 \
POSTPROCESS_LOWEST_MAGL=100 \
POSTPROCESS_SOURCE_LAYER_THICKNESS_M=100 \
OUTROOT=/net/fs06/d2/$USER/flexpart_outs \
./run_scripts/submit_slurm_array_backward.sh 2018020100 2018022800 20
```


### What is this repository for?

* This repository contains versions of the Lagrangian model FLEXPART
* Development versions
* Issues on the FLEXPART model, [tickets](https://gitlab.phaidra.org/flexpart/flexpart/-/issues)/[mail](mailto:gitlab.phaidra+flexpart-flexpart-456-issue-@univie.ac.at)
* Feature requests for future versions

## Getting started with Flexpart

The model is written in Fortran. It needs to be compiled for the architecture that runs it. Please have a look at the instructions on building FLEXPART available [here](./documentation/docs/building.md) or [online](https://flexpart.img.univie.ac.at/docs). There is also a containerized version of FLEXPART available.

## Downloading Meteorological Data

To run FLEXPART backward simulations, you need meteorological data files. The repository includes a script to download GFS (Global Forecast System) data from NOAA archives:

### Using the GFS Download Tool

The `tools/download_gfs_archive_python.py` script fetches historical GFS data from NOAA and generates an AVAILABLE file for FLEXPART.

**Installation requirements:**
- Python 3.6+
- For AWS S3 downloads (2021-01-01 onward): `pip install boto3`

**Basic usage:**

```bash
# Download data for a date range and create AVAILABLE file
./tools/download_gfs_archive_python.py \
  --start 2018021700 \
  --end 2018022800 \
  --outdir /path/to/gfs_data \
  --available /path/to/gfs_data/AVAILABLE
```

**Options:**
- `--start YYYYMMDDHH`: Start time (required)
- `--end YYYYMMDDHH`: End time, inclusive (required)
- `--outdir DIR`: Output directory for meteorological files (default: `./inputs`)
- `--available FILE`: Path to AVAILABLE file (default: `./AVAILABLE`)
- `--step-hours N`: Time step in hours (default: 3)
- `--source {auto,aws,nomads,ncei}`: Data source (default: `auto`)
  - `auto`: AWS S3 for dates >= 2021-01-01, NCEI for earlier dates
  - `aws`: AWS OpenData (0.25°, 2021-01-01+)
  - `ncei`: NCEI archives (0.5°, historical)
  - `nomads`: NOAA NOMADS (recent data only)
- `--force`: Re-download files even if they exist
- `--retries N`: Number of retry attempts (default: 3)
- `--timeout SEC`: Download timeout in seconds (default: 60)
- `--dry-run`: Preview without downloading

**Example: Download recent data with auto-selected source**

```bash
./tools/download_gfs_archive_python.py \
  --start 2026030100 \
  --end 2026040200 \
  --outdir ./gfs_data \
  --available ./gfs_data/AVAILABLE
```

### Running FLEXPART Backward Simulations

Once meteorological data is available, use the batch runner script:

```bash
./run_scripts/run_backward_batch.py \
  --domain EASTASIA \
  --receptor GSN \
  --end-time 2018022800 \
  --num-particles 10000 \
  --days 20 \
  --gfs-data /path/to/gfs_data
```

For more options, run:
```bash
./run_scripts/run_backward_batch.py --help
```

### Running Slurm Array Jobs (Current Workflow)

Current default workflow on Svante is:
- Use a precompiled FLEXPART executable (typically `src/FLEXPART`).
- Submit date-range array jobs via the Slurm scripts in `run_scripts/`.
- Manage run settings in `run_scripts/slurm_array_config.sh`.

Recommended submission flow:

```bash
# 1) Edit run parameters
vi run_scripts/slurm_array_config.sh

# 2) Submit
./run_scripts/run_slurm_array_backward.sh
```

Optional direct submit (without config file):

```bash
DOMAIN=EASTASIA \
RECEPTOR=GSN \
STEP_HOURS=1 \
BACKWARD_DAYS=20 \
NUM_PARTICLES=20000 \
POSTPROCESS_LOWEST_MAGL=100 \
POSTPROCESS_SOURCE_LAYER_THICKNESS_M=100 \
OUTROOT=/net/fs06/d2/$USER/flexpart_outs \
FLEXPART_EXE=/home/lwestern/work/flexpart_gfs/src/FLEXPART \
./run_scripts/submit_slurm_array_backward.sh 2018020100 2018022823 30
```

Notes:
- The date range is inclusive.
- One array task is created per `STEP_HOURS` timestamp.
- `%30` in the submit command above is the maximum concurrent array tasks.


### Contribution guidelines

* The version contributed should compile on a reference version of the system and compiler. 
   - `FLEXPART 10.4` used as a reference gfortran 5.4 on Ubuntu 16.04
   - `FLEXPART 11` uses as a reference gfortran 8.5.0 on AlmaLinux 8/RockyLinux 8 or gfortran 11.4.1 on RockyLinux 9

* Code contribution including new features and bug fixes should be complemented with appropriate tests
   An essential test consists of a set of input files and directories that allow FLEXPART to run.
   A test can be accompanied by output files for verification
* Code review
* report issues via mail to [support](mailto:gitlab.phaidra+flexpart-flexpart-456-issue-@univie.ac.at)
* become an active developer and request a user account on gitlab.
