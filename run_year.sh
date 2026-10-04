#!/bin/bash
# Usage: run_year.sh YEARS [models...]   e.g. run_year.sh 2023  |  run_year.sh 2021-2025
# Streams GFS-initialised AI models (sequential, 16 threads) and GFS (8 threads, in parallel). Resumable.
cd "$(dirname "$0")/scripts"
Y=$1; shift
MODELS=${@:-graphcast pangu fourcastnetv2}
PY=../.venv/bin/python
$PY -u 04_gfs_idx.py --years $Y --threads 8 >> ../logs/04_gfs_$Y.log 2>&1 &
G=$!
$PY -u 03_stream_mlwp.py $MODELS --years $Y --threads 16 >> ../logs/03_ai_$Y.log 2>&1
wait $G
echo done
