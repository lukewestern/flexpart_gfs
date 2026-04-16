![](https://www.flexpart.eu/chrome/site/flexpart_banner.png)
# Welcome to Flexpart - The Lagrangian particle dispersion model

> **Note:** This is a **forked and adapted version** of FLEXPART optimized for backward simulations of inert particles using GFS meteorological data. Currently configured primarily for use on the **MIT Svante HPC system**. The codebase includes tools for seamless GFS data acquisition and batch processing of backward trajectories.

The the main development site is @ University of Vienna.

Other references:

- [FLEXPART.eu](https://flexpart.eu)
- [FLEXPART@NILU](https://git.nilu.no/flexpart/flexpart)

### MIT/Svante install instructions

The FLEXPART binary must be built **on an EDR compute node** so it is linked against the correct GLIBC and runtime libraries for that architecture.  Building on a login node or fsXX node will produce a binary that fails at runtime with `GLIBC_2.2x not found`.

#### 1. Start an interactive EDR session

```bash
srun -n 1 -p edr --pty /bin/bash
```

#### 2. Build ecCodes locally with the system Fortran compiler

The EDR nodes ship GCC 6 as the default compiler, which is too old to build ecCodes 2.34+ (requires C++17).  Use GCC 11 from the module system, but it requires some missing runtime libraries that can be borrowed from the Julia module.

```bash
# Load required modules
source ~/.bashrc
module purge
module load cmake/3.26.4 gcc/11.3.0

# GCC 11's Fortran frontend needs libmpfr.so.6; find and expose it
EXTRA_LIBS="$(
  for lib in libmpfr.so.6 libmpc.so.3; do
    find /home/software /usr -name "$lib" 2>/dev/null | head -n 1 | xargs -r dirname
  done | awk 'NF && !seen[$0]++' | paste -sd:
)"
export LD_LIBRARY_PATH="${EXTRA_LIBS}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Verify GCC 11 is active
hash -r
gfortran --version   # should say GNU Fortran (GCC) 11.3.0
```

Clone/install ecbuild (needed by ecCodes 2.23 CMake system):

```bash
mkdir -p "$HOME/local/src" "$HOME/local/build"
cd "$HOME/local/src"
git clone --depth 1 --branch 3.8.2 https://github.com/ecmwf/ecbuild.git

mkdir -p "$HOME/local/build/ecbuild-build"
cd "$HOME/local/build/ecbuild-build"
cmake "$HOME/local/src/ecbuild" -DCMAKE_INSTALL_PREFIX="$HOME/local/ecbuild"
cmake --build . -j 8 && cmake --install .
```

Download and build ecCodes 2.23.0 with GCC 11:

```bash
GCC11=/home/software/rhel/8/gcc/11.3.0

cd "$HOME/local/src"
# Download tarball — check https://confluence.ecmwf.int/display/ECC/Releases for latest 2.23.x
curl -LO https://confluencehpc.ecmwf.int/eccodes-2.23.0-Source.tar.gz
# or use: wget https://github.com/ecmwf/eccodes/releases/...
tar xzf eccodes-2.23.0-Source.tar.gz
mv eccodes-2.23.0-Source eccodes-2.23.0

mkdir -p "$HOME/local/build/eccodes-gcc11-build"
cd "$HOME/local/build/eccodes-gcc11-build"

cmake "$HOME/local/src/eccodes-2.23.0" \
  -DCMAKE_INSTALL_PREFIX="$HOME/local/eccodes-gcc11" \
  -DCMAKE_C_COMPILER="$GCC11/bin/gcc" \
  -DCMAKE_CXX_COMPILER="$GCC11/bin/g++" \
  -DCMAKE_Fortran_COMPILER="$GCC11/bin/gfortran" \
  -DCMAKE_PREFIX_PATH="$HOME/local/ecbuild" \
  -DENABLE_FORTRAN=ON \
  -DBUILD_SHARED_LIBS=ON \
  -DENABLE_AEC=OFF \
  -DENABLE_NETCDF=OFF \
  -DENABLE_JPG=OFF \
  -DENABLE_PNG=OFF

cmake --build . -j 8
cmake --install .
```

> **Note:** GCC 11 on EDR nodes links against a newer libgfortran than the system GLIBC provides.  
> If the ecCodes build fails with `undefined reference to getentropy@GLIBC_2.25` during the CMake compiler test, GCC 11 is not runtime-compatible here.  
> In that case fall back to building ecCodes 2.23.0 with GCC 6 (`FC=/usr/bin/gfortran`) — the resulting `grib_api.mod` will be usable for the FLEXPART Fortran build.  
> The FLEXPART source already contains the necessary `#ifdef _OPENMP` guards so the rank-8 array limitation of GCC 6 is not triggered in serial builds.

#### 3. Set up build environment and compile FLEXPART

```bash
# Unload GCC 11 if loaded; use system GCC 6 for the FLEXPART link step
module purge
hash -r

export ECCODES_PREFIX="$HOME/local/eccodes-gcc11"   # or eccodes-gcc6 if fallback
export GF=/home/software/rhel/8/gcc/11.3.0/bin/gfortran  # or /usr/bin/gfortran for fallback

export CPATH="/usr/include:${ECCODES_PREFIX}/include:/usr/lib64/gfortran/modules"
export LIBRARY_PATH="/usr/lib64:${ECCODES_PREFIX}/lib64:${ECCODES_PREFIX}/lib"
unset CONDA_PREFIX
unset LD_LIBRARY_PATH   # avoid conda RPATH contamination

cd flexpart_gfs/src
make -f makefile_svante cleanall
make -f makefile_svante eta=no ncf=yes SERIAL=yes GFLAG= -j4 FC="$GF" F90="$GF"
```

`GFLAG=` disables the `-g` debug-info flag which triggers assembler errors with the older `as` on EDR nodes.  
`SERIAL=yes` disables OpenMP (correct for EDR runs).

#### 4. Validate the binary is runtime-compatible

```bash
# GLIBC requirement — should show 2.14 or lower, not 2.25/2.27
readelf -V src/FLEXPART | grep 'Name: GLIBC_' | sort -u

# Should NOT contain conda or home directory lib paths
readelf -d src/FLEXPART | grep -E 'RPATH|RUNPATH' || true

# All libraries should resolve from /usr/lib64 and $HOME/local/eccodes-gcc*/lib*
ldd src/FLEXPART | grep -E 'eccodes|netcdf|gfortran|libc\.so'
```

A passing result looks like:

```
Name: GLIBC_2.14  ...
Name: GLIBC_2.2.5 ...
RPATH: [/usr/lib64:/home/lwestern/local/eccodes-gcc.../lib64:...]
libnetcdff.so.6 => /usr/lib64/libnetcdff.so.6
libeccodes.so   => /home/.../eccodes-gcc.../lib64/libeccodes.so
libgfortran.so.3 => /usr/lib64/libgfortran.so.3
libc.so.6       => /usr/lib64/libc.so.6
```


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
