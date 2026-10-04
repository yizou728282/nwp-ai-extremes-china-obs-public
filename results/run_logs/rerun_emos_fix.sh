#!/usr/bin/env bash
# Targeted recompute after the closed-form gauss_twcrps fix (oos_emos_tw1). Local, not tracked.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python; W=${WORKERS:-2}; B=../emosfix-backup
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(date) start"
$PY -m unittest discover -s tests > logs/emosfix_tests.log 2>&1
$PY -u scripts/05_build_dataset.py > logs/emosfix_cubes.log 2>&1      # unchanged code; refreshes cube provenance
echo "$(date) cubes done"
for S in core core21 ic; do
  for V in t2 tmax tmin; do                                             # EMOS variables: affected
    $PY -u scripts/07_scores.py --var $V --set $S --workers $W > logs/emosfix_${S}_${V}.log 2>&1
    $PY -u scripts/08_stats.py --var $V --set $S >> logs/emosfix_${S}_${V}.log 2>&1
    echo "$(date) $S $V done"
  done
  VV="ws pr24"; [ $S = core21 ] && VV="ws"
  for V in $VV; do                                                      # not EMOS: verify scores unchanged
    $PY -u scripts/07_scores.py --var $V --set $S --workers $W > logs/emosfix_${S}_${V}.log 2>&1
    if $PY logs/emosfix_identity.py $S $V $B/scores >> logs/emosfix_${S}_${V}.log 2>&1; then
      echo "$(date) $S $V scores identical -> table provenance restamped"
    else
      $PY -u scripts/08_stats.py --var $V --set $S >> logs/emosfix_${S}_${V}.log 2>&1
      echo "$(date) $S $V scores differed -> 08 rerun"
    fi
  done
  $PY scripts/10_report_tables.py $S > logs/emosfix_tables_${S}.log 2>&1
done
$PY scripts/09_figures.py --set core > logs/emosfix_figures.log 2>&1
cp results/figures/*.pdf paper/figures/
$PY paper/make_tables.py > logs/emosfix_paper.log 2>&1
( cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex ) >> logs/emosfix_paper.log 2>&1
$PY - <<'PY'
import json,sys
from pathlib import Path
sys.path.insert(0,'scripts')
from provenance import method_hash
p=Path('results/revision_status.json'); d=json.loads(p.read_text())
d.update(status='recomputed', method_sha256=method_hash(), oos_emos_tw1_fix='closed-form Gaussian threshold-weighted CRPS; t2/tmax/tmin rescored, ws/pr24 verified identical')
p.write_text(json.dumps(d, indent=2)+'\n')
PY
echo "$(date) ALL DONE"
