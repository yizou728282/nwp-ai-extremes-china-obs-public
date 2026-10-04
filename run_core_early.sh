#!/bin/bash
# Starts the core / core21 analysis as soon as all GFS-initialised streaming (2021-2025/08) is on disk,
# in parallel with the IFS-init streaming of run_all.sh. Writes logs/core_done.marker when finished.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
need="data/fc/fourcastnetv2/2021-12.npz data/fc/gfs/2021-12.npz data/fc/graphcast/2025-08.npz data/fc/fourcastnetv2/2025-08.npz data/fc/gfs/2025-08.npz data/fc/pangu/2024-12.npz"
ready() { for f in $need; do [ -f "$f" ] || return 1; done; return 0; }
until ready; do sleep 120; done
sleep 300   # let the GFS job of the last year finish writing
echo "$(date) core-early: start"
cd scripts && ../.venv/bin/python -u 04_gfs_idx.py --years 2021-2025 --threads 8 >> ../logs/04_gfs_fill.log 2>&1; cd ..
WORKERS=4 ./run_analysis.sh core >> logs/analysis_core.log 2>&1 && echo "$(date) core done"
SKIP_BUILD=1 WORKERS=4 ./run_analysis.sh core21 t2 tmax tmin ws >> logs/analysis_core21.log 2>&1 && echo "$(date) core21 done"
for S in core core21; do .venv/bin/python scripts/10_report_tables.py $S > /dev/null; done
touch logs/core_done.marker; echo "$(date) core-early: done"
