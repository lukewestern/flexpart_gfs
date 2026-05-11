![](https://www.flexpart.eu/chrome/site/flexpart_banner.png)
# Welcome to Flexpart - The Lagrangian particle dispersion model

> **Note:** This is a **forked and adapted version** of FLEXPART optimized for backward simulations of inert particles using CFSv2 meteorological data. Currently configured primarily for use on the **MIT Svante HPC system**.

The the main development site is @ University of Vienna.

Other references:

- [FLEXPART.eu](https://flexpart.eu)
- [FLEXPART@NILU](https://git.nilu.no/flexpart/flexpart)

### MIT/Svante install instructions

#### Quick setup script (recommended)

Most of the Svante setup/build flow can be automated with:

```bash
./run_scripts/setup_svante_environment.sh
```

What it does:
- builds ecbuild + ecCodes (compiler/toolchain configurable)
- builds FLEXPART with portable flags (`arch=x86-64`, `SERIAL=yes`)
- creates/updates a stable postprocess conda env (`flexpart-post`, Python 3.12)

Useful options:

```bash
# Preview commands only
./run_scripts/setup_svante_environment.sh --dry-run

# Skip ecCodes build and only rebuild FLEXPART + postprocess env
./run_scripts/setup_svante_environment.sh --skip-eccodes

# Force a specific toolchain if your build node differs
./run_scripts/setup_svante_environment.sh \
  --gcc-prefix /home/software/rhel/8/gcc/11.3.0 \
  --fortran-compiler /usr/bin/gfortran

# Skip postprocess env creation
./run_scripts/setup_svante_environment.sh --no-postprocess-env
```

Then run the Slurm workflow:

Edit
```bash
./run_scripts/slurm_array_config.sh
```
Then run 
```bash
./run_scripts/run_slurm_array_backward.sh
./run_scripts/run_slurm_postprocess_all.sh
```

> The detailed manual steps for troubleshooting and custom builds are at the end of this README.


### Running Slurm Jobs (Svante/HPC)

Current recommended production workflow is a 2-stage Slurm process:

1. Run FLEXPART array jobs (write grid outputs only)
2. Run one postprocess Slurm job over all completed outputs

This separation is intentional and avoids postprocess crashes from unstable Python runtime stacks on some nodes.

#### Stage 1: FLEXPART array jobs

Edit `run_scripts/slurm_array_config.sh` and set at least:
- `START_DATE`, `END_DATE`, `DOMAIN`, `RECEPTOR`
- `DISABLE_AUTO_POSTPROCESS="1"`
- `PRUNE_TO_GRID_FILES="1"`

Then submit:

```bash
./run_scripts/run_slurm_array_backward.sh
```

This creates one run directory per timestamp under:
- `/net/fs06/d2/$USER/flexpart_outs`

With `PRUNE_TO_GRID_FILES="1"`, each run keeps only `output/grid_time_*.nc` after FLEXPART completion.

#### Stage 2: postprocess all runs in one Slurm job

Submit:

```bash
./run_scripts/run_slurm_postprocess_all.sh
```

Recommended one-time setup for a stable postprocess environment:

```bash
conda create -n flexpart-post python=3.12 -y
conda activate flexpart-post
python -m pip install --upgrade pip
python -m pip install numpy pandas xarray netCDF4
```

Then set in `run_scripts/slurm_array_config.sh`:
- `POSTPROCESS_PYTHON_CMD="/home/$USER/.conda/envs/flexpart-post/bin/python"`
- `POSTPROCESS_DRIVER_PYTHON="python3"` (or explicit path if needed)

This runs `run_scripts/postprocess_all_outputs.py` over discovered NetCDF run outputs in `OUTROOT`:
- `output/grid_time_*.nc`

Default postprocess-all behavior:
- writes final footprint files to `FINAL_DIR` (default: `OUTROOT`)
- moves final files to top-level output directory
- deletes per-run FLEXPART folders after successful move

The final files are named:
- `*_FLEXPART_CFSv2_<DOMAIN>_inert_<YYYYMMDDHH>.nc`

Key postprocess-all config variables in `run_scripts/slurm_array_config.sh`:
- `POSTPROCESS_DRIVER_PYTHON` (driver interpreter, default `python3`)
- `POSTPROCESS_PYTHON_CMD` (python used by `postprocess_footprint.py` subprocess)
- `FINAL_DIR` (destination for final footprint files)
- `POSTPROCESS_OVERWRITE` (`1` to replace existing final files)
- `KEEP_RUN_DIRS` (`1` to keep run directories)
- `LIMIT` (`>0` to process first N discovered files)
- `DRY_RUN` (`1` for preview-only)

#### Direct submit helpers

```bash
# FLEXPART array
./run_scripts/submit_slurm_array_backward.sh 2018020100 2018022800 20

# Postprocess-all single job
./run_scripts/run_slurm_postprocess_all.sh
```


### What is this repository for?

* This repository contains versions of the Lagrangian model FLEXPART
* Development versions
* Issues on the FLEXPART model, [tickets](https://gitlab.phaidra.org/flexpart/flexpart/-/issues)/[mail](mailto:gitlab.phaidra+flexpart-flexpart-456-issue-@univie.ac.at)
* Feature requests for future versions

## Getting started with Flexpart

The model is written in Fortran. It needs to be compiled for the architecture that runs it. Please have a look at the instructions on building FLEXPART available [here](./documentation/docs/building.md) or [online](https://flexpart.img.univie.ac.at/docs). There is also a containerized version of FLEXPART available.

## Downloading Meteorological Data (Manual GDEX Workflow)

Meteorological data acquisition is now manual. A user must download CFSv2 files from GDEX.

### 1) Download raw CFSv2 files from GDEX

1. Log into https://gdex.ucar.edu/datasets/d094000/ to access the NCEP Climate Forecast System Version 2 (CFSv2) 6-hourly products.
2. Click **Data Access** then **Get a Subset**.
3. From the Parameter Selection options, select **FLEXPART Model Input: 6-hour Forecasts**.
4. Select the time range you want.
5. Click **Continue**.
6. Under **Gridded Products**, click **All available** and then **Update/Click Here** at the bottom.
7. Click **1-hour Forecast** through **6-hour Forecast**.
8. Click **Click Here** again.
9. Choose either file type, then click **Submit**.
10. Wait for the email with links to download the meteorological data.

Once the download is complete, process these data with the tools in the `/tools/` folder.

### 2) Prepare FLEXPART inputs and AVAILABLE file

The repository now contains preparation scripts in `tools/`:
- `tools/prepare_flexpart_gfs.sh` (serial)
- `tools/prepare_flexpart_gfs_parallel.sh` (parallel, recommended)

Example (parallel):

```bash
./tools/prepare_flexpart_gfs_parallel.sh \
  /path/to/raw_cfsv2 \
  /path/to/flexpart_inputs \
  /path/to/flexpart_inputs/AVAILABLE
```

Example (serial):

```bash
./tools/prepare_flexpart_gfs.sh \
  /path/to/raw_cfsv2 \
  /path/to/flexpart_inputs \
  /path/to/flexpart_inputs/AVAILABLE
```

Useful environment variables:
- `NPROC` (parallel script worker count)
- `REPACK_TO_SIMPLE` (default `1`)
- `OVERWRITE_EXISTING` (default `0`)
- `STEP_FILTER_MODE` (default `nonzero_max`)

### Running FLEXPART Backward Simulations

Once meteorological data is available, use the batch runner script:

```bash
./run_scripts/run_backward_batch.py \
  --domain EASTASIA \
  --receptor GSN \
  --end-time 2018022800 \
  --num-particles 10000 \
  --days 20 \
  --gfs-data /path/to/flexpart_inputs
```

For more options, run:
```bash
./run_scripts/run_backward_batch.py --help
```

### Slurm Workflow Summary

Minimal run sequence:

```bash
# 1) Configure
vi run_scripts/slurm_array_config.sh

# 2) Run FLEXPART array jobs (stage 1)
./run_scripts/run_slurm_array_backward.sh

# 3) Run postprocess-all job (stage 2)
./run_scripts/run_slurm_postprocess_all.sh
```

Notes:
- Date range (`START_DATE`..`END_DATE`) is inclusive.
- FLEXPART stage creates one array task per `STEP_HOURS` timestamp.
- Postprocess stage is one Slurm job that scans `OUTROOT` for all matching run outputs.


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

### Detailed install instructions

The FLEXPART binary must be built **on an EDR compute node** so it is linked against the correct GLIBC and runtime libraries for that architecture.  Building on a login node or fsXX node will produce a binary that fails at runtime with `GLIBC_2.2x not found`.

#### 1. Start an interactive EDR session

```bash
srun -n 1 -p edr --pty /bin/bash
```

#### 2. Build ecCodes locally with a compiler compatible with your build node

The exact compiler/runtime stack can vary by node/partition. Use the same node class for build and run. If your selected compiler and ecCodes Fortran library do not agree on `libgfortran` soname, FLEXPART may compile but fail or misbehave at runtime.

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

# Verify selected Fortran compiler
hash -r
gfortran --version
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

Download and build ecCodes 2.23.0 with your selected compiler:

```bash
GCC_PREFIX=/home/software/rhel/8/gcc/11.3.0

cd "$HOME/local/src"
# Download tarball — check https://confluence.ecmwf.int/display/ECC/Releases for latest 2.23.x
curl -LO https://confluencehpc.ecmwf.int/eccodes-2.23.0-Source.tar.gz
# or use: wget https://github.com/ecmwf/eccodes/releases/...
tar xzf eccodes-2.23.0-Source.tar.gz
mv eccodes-2.23.0-Source eccodes-2.23.0

mkdir -p "$HOME/local/build/eccodes-build"
cd "$HOME/local/build/eccodes-build"

cmake "$HOME/local/src/eccodes-2.23.0" \
  -DCMAKE_INSTALL_PREFIX="$HOME/local/eccodes-gcc11" \
  -DCMAKE_C_COMPILER="$GCC_PREFIX/bin/gcc" \
  -DCMAKE_CXX_COMPILER="$GCC_PREFIX/bin/g++" \
  -DCMAKE_Fortran_COMPILER="$GCC_PREFIX/bin/gfortran" \
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

> If GCC 11 is not runtime-compatible on your node, build ecCodes and FLEXPART with `/usr/bin/gfortran` (or another matching compiler) instead.

#### 3. Set up build environment and compile FLEXPART

```bash
# Ensure a clean shell and select the compiler used for FLEXPART
module purge
hash -r

export ECCODES_PREFIX="$HOME/local/eccodes-gcc11"   # or eccodes-gcc6 if fallback
export GF="$(command -v gfortran)"  # or explicit path

export CPATH="/usr/include:${ECCODES_PREFIX}/include:/usr/lib64/gfortran/modules"
export LIBRARY_PATH="/usr/lib64:${ECCODES_PREFIX}/lib64:${ECCODES_PREFIX}/lib"
unset CONDA_PREFIX
unset LD_LIBRARY_PATH   # avoid conda RPATH contamination

# ABI preflight: ecCodes Fortran and compiler must use same libgfortran soname
ldd "${ECCODES_PREFIX}/lib64/libeccodes_f90.so" | grep libgfortran
printf 'program x\nprint *,1\nend program x\n' > /tmp/fp_abi.f90
"${GF}" /tmp/fp_abi.f90 -o /tmp/fp_abi
ldd /tmp/fp_abi | grep libgfortran
rm -f /tmp/fp_abi /tmp/fp_abi.f90

cd flexpart_gfs/src
make -f makefile_svante cleanall
make -f makefile_svante eta=no ncf=yes SERIAL=yes arch=x86-64 GFLAG= -j4 FC="$GF" F90="$GF"
```

`GFLAG=` disables the `-g` debug-info flag which triggers assembler errors with the older `as` on EDR nodes.  
`SERIAL=yes` disables OpenMP (correct for EDR runs).
`arch=x86-64` avoids node-specific instructions so the binary is portable across EDR/FDR.

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
