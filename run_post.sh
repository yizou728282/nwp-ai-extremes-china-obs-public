#!/bin/bash
# Follow-on chain after run_all.sh (pid $1) finishes: fill gaps (resumable), IFS-init streaming (resume),
# rebuild cubes incl. IFS-init models, initial-condition decomposition (ic set), final figures and tables.
# core/core21 are produced earlier by run_core_early.sh (waits for its marker).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -n "$1" ] && while kill -0 "$1" 2>/dev/null; do sleep 60; done
echo "$(date) post: start"
cd scripts; PY=../.venv/bin/python
$PY -u 04_gfs_idx.py --years 2021-2025 --threads 8 >> ../logs/04_gfs_fill.log 2>&1
$PY -u 03_stream_mlwp.py graphcast pangu fourcastnetv2 --years 2021-2025 --threads 16 >> ../logs/03_ai_fill.log 2>&1
$PY -u 03_stream_mlwp.py graphcast_ifs pangu_ifs fourcastnetv2_ifs --years 2022-2025 --threads 16 >> ../logs/03_ai_ifs.log 2>&1
echo "$(date) post: streaming complete"; du -sh ../data
cd ..
while [ ! -f logs/core_done.marker ]; do sleep 120; done
WORKERS=4 ./run_analysis.sh ic >> logs/analysis_ic.log 2>&1; echo "$(date) post: ic done"
.venv/bin/python scripts/09_figures.py --set core >> logs/analysis_core.log 2>&1
for S in core core21 ic; do .venv/bin/python scripts/10_report_tables.py $S > /dev/null 2>> logs/run_post.log; done
echo "$(date) post: ALL DONE"
