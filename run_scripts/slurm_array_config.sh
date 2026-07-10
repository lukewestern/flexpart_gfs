#!/bin/bash
# Optional overrides for run_scripts/run_slurm_array_backward.sh
# Use this file to avoid long inline env var command lines.

# Inclusive date range, format YYYYMMDDHH
START_DATE="2019010100"
END_DATE="2020010100"
MAX_CONCURRENT="20"

DOMAIN="WESTUSA"
RECEPTOR="THD"
STEP_HOURS="1"
BACKWARD_DAYS="20"
NUM_PARTICLES="20000"
IPOUT="1"
# Optional: force NXSHIFT. Leave empty to auto-select in run_backward_batch.py.
# NXSHIFT="25"
LSUBGRID="1"
LINIT_COND="1"
# Preferred name (legacy alias POSTPROCESS_LOWEST_MAGL is still accepted by scripts)
POSTPROCESS_FOOTPRINT_OUTHEIGHT_M="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
DISABLE_AUTO_POSTPROCESS="1"
PRUNE_TO_GRID_FILES="0"
OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"

# Absolute path to compiled FLEXPART executable.
FLEXPART_EXE="/home/lwestern/work/flexpart_gfs/src/FLEXPART"

# Optional runtime controls
PYTHON_CMD=""
POSTPROCESS_PYTHON_CMD="/home/${USER}/.conda/envs/flexpart-post/bin/python"
USE_PROJECT_VENV="0"
DEBUG_ENV="0"

# Postprocess-all (single Slurm job) options
POSTPROCESS_DRIVER_PYTHON="python3"
POSTPROCESS_WORKERS="8"
# Max concurrent month tasks for run_slurm_postprocess_monthly_array.sh
POSTPROCESS_MONTHLY_MAX_CONCURRENT="6"
FINAL_DIR="${OUTROOT}/${RECEPTOR,,}_hourly"
POSTPROCESS_OVERWRITE="0"
KEEP_RUN_DIRS="1"
LIMIT="0"
WRITE_MONTHLY="1"
MONTHLY_DIR="${OUTROOT}/${RECEPTOR,,}"
DRY_RUN="0"
