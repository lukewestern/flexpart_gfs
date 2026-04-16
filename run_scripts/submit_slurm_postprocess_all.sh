#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${SCRIPT_DIR}/slurm_postprocess_all.sh"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ ! -f "${JOB_SCRIPT}" ]]; then
  echo "ERROR: job script not found: ${JOB_SCRIPT}"
  exit 2
fi

echo "Submitting postprocess-all job"
echo "  FLEXPART_REPO_ROOT=${REPO_ROOT}"
echo "  OUTROOT=${OUTROOT:-/net/fs06/d2/${USER}/flexpart_outs}"
echo "  FINAL_DIR=${FINAL_DIR:-<same as OUTROOT>}"

sbatch \
  --export=ALL,FLEXPART_REPO_ROOT="${REPO_ROOT}" \
  "${JOB_SCRIPT}"
