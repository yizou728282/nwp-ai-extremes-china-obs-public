#!/usr/bin/env bash
# Recompute every result affected by the scientific corrections using existing forecast-point caches.
# Build missing forecast caches with README's download commands first (several hours).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-.venv/bin/python}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"$PY" - <<'PY'
from pathlib import Path
import pandas as pd
models = ('gfs', 'graphcast', 'pangu', 'fourcastnetv2', 'graphcast_ifs', 'pangu_ifs', 'fourcastnetv2_ifs')
missing = [m for m in models if not any((Path('data/fc') / m).glob('*.npz'))]
if missing:
    raise SystemExit('Missing forecast caches: ' + ', '.join(missing) + '\nFollow README forecast-download commands, then rerun bash run_corrected.sh.')
months = pd.period_range('2021-01', '2025-08', freq='M').astype(str)
for model in models:
    # Known archive gaps are within months, except entire IFS months absent in 2024.
    required = [m for m in months if (not model.endswith('_ifs') and (model != 'graphcast' or m >= '2022-01'))]
    absent = [m for m in required if not (Path('data/fc') / model / f'{m}.npz').exists()]
    if absent:
        raise SystemExit(f'Incomplete {model} monthly caches: ' + ', '.join(absent))
Path('logs').mkdir(exist_ok=True)
PY
"$PY" -m unittest discover -s tests -v
"$PY" - <<'PY'
import json
from pathlib import Path
p = Path('results/revision_status.json')
d = json.loads(p.read_text()); d['status'] = 'recomputation_in_progress'
p.write_text(json.dumps(d, indent=2) + '\n')
PY
"$PY" -u scripts/02_obs.py --years 2011-2025 --threads "${OBS_THREADS:-6}" > logs/corrected_observations.log 2>&1
"$PY" -u scripts/02_obs.py --summary >> logs/corrected_observations.log 2>&1
"$PY" -u scripts/02b_ghcnd.py > logs/corrected_ghcnd.log 2>&1
"$PY" -u scripts/05_build_dataset.py > logs/corrected_cubes.log 2>&1
for comparison in core core21 ic; do
    variables="t2 tmax tmin ws pr24"
    if [ "$comparison" = core21 ]; then variables="t2 tmax tmin ws"; fi
    for variable in $variables; do
        "$PY" -u scripts/07_scores.py --var "$variable" --set "$comparison" --workers "${WORKERS:-4}" > "logs/corrected_${comparison}_${variable}.log" 2>&1
        "$PY" -u scripts/08_stats.py --var "$variable" --set "$comparison" >> "logs/corrected_${comparison}_${variable}.log" 2>&1
    done
    "$PY" scripts/10_report_tables.py "$comparison" > "logs/corrected_tables_${comparison}.log" 2>&1
done
"$PY" scripts/09_figures.py --set core > logs/corrected_figures.log 2>&1
mkdir -p paper/figures
cp results/figures/*.pdf paper/figures/
"$PY" paper/make_tables.py
if [ -f paper/main.tex ] && command -v latexmk >/dev/null 2>&1; then
    ( cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex )
else
    echo "Scientific results and LaTeX tables are ready; manuscript compilation is optional and its source is not deposited in this public release."
fi
"$PY" - <<'PY'
import json,sys
from pathlib import Path
sys.path.insert(0, 'scripts')
from provenance import method_hash
p = Path('results/revision_status.json'); d = json.loads(p.read_text())
d.update(status='recomputed', method_sha256=method_hash(), manuscript_status='Recomputed tables and figures; inspect updated results and interpretation before submission.')
p.write_text(json.dumps(d, indent=2) + '\n')
PY
echo "Corrected numerical results and LaTeX tables built. Review updated results and interpretation."
