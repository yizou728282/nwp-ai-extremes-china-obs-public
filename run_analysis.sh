#!/bin/bash
# Usage: run_analysis.sh SET [vars]   builds cubes (once per call), scores, stats, figures. Logs in logs/analysis_SET.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")/scripts"; S=$1; shift; V=${@:-t2 tmax tmin ws pr24}; PY=../.venv/bin/python
[ "${SKIP_BUILD:-0}" = 1 ] || $PY -u 05_build_dataset.py
for v in $V; do $PY -u 07_scores.py --var $v --set $S --workers ${WORKERS:-4} && $PY -u 08_stats.py --var $v --set $S; done
$PY -u 09_figures.py --set $S
echo analysis $S done
