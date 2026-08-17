#!/bin/bash
# Config overrides for run_scripts_daily/run_slurm_array_backward_daily.sh
# Use this file to avoid long inline env var command lines.

# Inclusive date range, format YYYYMMDD (no hour — entire day is always run)
START_DATE="20210101"
END_DATE="20230101"
MAX_CONCURRENT="40"

DOMAIN="WESTUSA"
RECEPTOR="THD"
BACKWARD_DAYS="20"
NUM_PARTICLES="20000"

# IPOUT=1: write partoutput at every output interval so per-release exit-point
# analysis is possible. The daily runner slices particles by sequential index
# (release h uses particles [h*N, (h+1)*N-1]) to isolate each hour's group.
IPOUT="1"

# Optional: force NXSHIFT. Leave empty to auto-select.
# NXSHIFT="25"
LSUBGRID="1"
LINIT_COND="1"

POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"

# Enable postprocessing to derive exit-point CSVs for each hourly release.
# Preserve partoutput_*.nc and postprocessed footprints.
DISABLE_AUTO_POSTPROCESS="1"
# Skip rerunning a day if its run directory already exists under OUTROOT.
SKIP_IF_RUN_DIR_EXISTS="1"
# If postprocessed daily NetCDF(s) for a date already exist in output/, skip rerunning that day.
SKIP_IF_POSTPROCESSED_EXISTS="1"
# PRUNE_OUTPUTS: after postprocessing, delete raw FLEXPART outputs (PARTFXTR, DRYDEP, etc.)
# and partoutput_*.nc, keeping only postprocessed footprints, exit CSVs, and grid_time.
PRUNE_OUTPUTS="1"

OUTROOT="/net/fs06/d2/${USER}/flexpart_outs_daily"

# Absolute path to compiled FLEXPART executable.
FLEXPART_EXE="/home/lwestern/work/flexpart_gfs/src/FLEXPART"

PYTHON_CMD=""
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
USE_PROJECT_VENV="0"
DEBUG_ENV="0"
