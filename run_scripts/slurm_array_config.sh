#!/bin/bash
# Optional overrides for run_scripts/run_slurm_array_backward.sh
# Use this file to avoid long inline env var command lines.

# Inclusive date range, format YYYYMMDDHH
START_DATE="2018020100"
END_DATE="2018020200"
MAX_CONCURRENT="30"

DOMAIN="EASTASIA"
RECEPTOR="GSN"
STEP_HOURS="1"
BACKWARD_DAYS="20"
NUM_PARTICLES="20000"
POSTPROCESS_LOWEST_MAGL="100"
POSTPROCESS_SOURCE_LAYER_THICKNESS_M="100"
OUTROOT="/net/fs06/d2/${USER}/flexpart_outs"

# Absolute path to compiled FLEXPART executable.
FLEXPART_EXE="/home/lwestern/work/flexpart_gfs/src/FLEXPART"

# Optional runtime controls
PYTHON_CMD=""
USE_PROJECT_VENV="0"
DEBUG_ENV="0"
