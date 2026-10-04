# nwp-ai-extremes-china-obs-public

This public submission snapshot is identified by release `submission-2026-10-04`. The scientific results were generated previously; publication of this snapshot does not constitute a new full numerical rerun. Large inputs and score arrays must be reconstructed. See README.md for the exact included scope.

Code, derived results and manuscript sources for:

> Zou, Y. *Do AI weather models forecast extremes better than GFS? A station-based verification over China with a common initialization source, 2022–2025.* (submission prepared for *Asia-Pacific Journal of Atmospheric Sciences*)

GFS-initialised GraphCast Operational, Pangu-Weather and FourCastNet v2-small from the NOAA/CIRA AI weather prediction (AIWP) archive are compared with operational GFS forecasts at 334 ISD stations (249 for 24-h precipitation). The shared GFS F000 input source controls initialization source; it does not isolate a pure architecture effect from training and preprocessing. Paired IFS-initialised runs measure initialization sensitivity conditional on each trained system. Deterministic forecasts are scored with potential CRPS (PCRPS, via isotonic distributional regression/EasyUQ) and threshold-weighted PCRPS. Independent-period IDR and EMOS, calendar-block bootstrap intervals and station-wise DM/BH-FDR tests complement these in-sample measures.

**Revision status.** The code was corrected on 2026-10-03 and all verification results were recomputed from the forecast-point caches with the corrected pipeline (`bash run_corrected.sh`). Tables, figures and the quantitative statements in the manuscript correspond to the current method hash recorded in `results/revision_status.json`. See [REVISION_NOTES.md](REVISION_NOTES.md).

With complete forecast-point caches present, run:

```bash
bash run_corrected.sh
```

This refreshes old observation caches, rebuilds thresholds/cubes, recomputes all comparison sets, regenerates figures/tables, and optionally compiles a privately supplied manuscript. It stops on missing caches or failed steps. Inspect the updated numerical interpretation before submission.

Variables: 2-m temperature, daily Tmax/Tmin (12–12 UTC), 10-m wind speed, 24-h precipitation (GraphCast only among the AI models). Period: 00 UTC initialisations 2021-01-01 to 2025-08-12, leads 6–240 h every 6 h.

## Repository layout

```
scripts/            processing and verification pipeline (01–10, common.py)
run_*.sh            drivers for streaming and analysis
data/interim/       station lists (superset, analysis set with climate zone and variable flags, coverage)
data/aux/           rebuilt by the scripts; auxiliary maps are not included
data/refs/          rebuilt by 03_refs.py; archive chunk indices are not included
results/tables/     summary scores, bootstrap CIs, DM/FDR results, key-result tables (csv/md)
paper/figures/     seven verified submission figures (PDF)
paper/              result table generator, generated LaTeX tables and verified figure PDFs; manuscript files excluded
```

## Data sources

All inputs are public and are downloaded by the scripts; no account is needed.

| Data | Source |
|---|---|
| AI model forecasts (GraphCast Operational: GRAP_v100; Pangu: PANG_v100; FourCastNet v2-small: FOUR_v200) | NOAA/CIRA AIWP archive, `s3://noaa-oar-mlwp-data` (https://registry.opendata.aws/aiwp/), Radford et al. (2025) |
| GFS 0.25° deterministic forecasts | `s3://noaa-gfs-bdp-pds` (https://registry.opendata.aws/noaa-gfs-bdp-pds/), read with byte ranges from the `.idx` files |
| Station observations | NOAA NCEI Integrated Surface Database (global-hourly) and station history/inventory |
| Daily Tmax/Tmin/precipitation (gap filling, cross-checks) | GHCN-Daily, `s3://noaa-ghcn-pds` |
| Climate zones | Köppen–Geiger 0.1° map, Beck et al. (2023), Sci. Data 10, 724 (CC BY 4.0) |
| Coastlines | Natural Earth 1:50m (public domain) |

## Requirements

Python 3.12 and the packages in `requirements.txt`; everything runs on CPU.

```bash
uv venv .venv -p 3.12 && uv pip install -p .venv/bin/python -r requirements.txt
# or: python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

This public data/code release excludes manuscript PDFs, manuscript TeX bodies, bibliography/template files and the manuscript source ZIP. `paper/make_tables.py` generates the numerical LaTeX result tables from the deposited CSVs; it does not require the private manuscript. The reconstruction driver skips manuscript compilation when the private source or a TeX compiler is absent. Scientific scoring, bootstrap and FDR calculations are unchanged.

## Reproducing the study

All commands are run from the repository root unless stated otherwise. Every streaming step is resumable.

```bash
cd scripts
../.venv/bin/python 01_select_stations.py                    # station superset + Köppen–Geiger zones
../.venv/bin/python 02_obs.py --years 2011-2025               # ISD observations (streamed, parsed in memory)
../.venv/bin/python 02_obs.py --summary                       # final analysis set (334 stations)
../.venv/bin/python 02b_ghcnd.py                              # GHCN-Daily gap filling and checks
../.venv/bin/python 03_refs.py GRAP_v100_GFS PANG_v100_GFS FOUR_v200_GFS --years 2021-2025   # required: chunk tables are not bundled
../.venv/bin/python 03_refs.py GRAP_v100_IFS PANG_v100_IFS FOUR_v200_IFS --years 2022-2025
cd ..

./run_year.sh 2021-2025                # GFS-initialised AI models + GFS at station points
cd scripts && ../.venv/bin/python -u 03_stream_mlwp.py graphcast_ifs pangu_ifs fourcastnetv2_ifs --years 2022-2025 && cd ..

./run_analysis.sh core                                    # build cubes, scores, statistics, figures
SKIP_BUILD=1 ./run_analysis.sh core21 t2 tmax tmin ws     # sensitivity set incl. 2021 (no GraphCast)
SKIP_BUILD=1 ./run_analysis.sh ic                         # adds IFS-initialised models (initial-condition decomposition)
for S in core core21 ic; do .venv/bin/python scripts/10_report_tables.py $S; done
```

`run_all.sh`, `run_core_early.sh` and `run_post.sh` chain the same steps for unattended runs (streaming year by year, gap filling, then analysis). Logs are written to `logs/`.

Comparison sets (`07_scores.py --set`):
* `core`: GFS + the three GFS-initialised AI models, initialisations 2022-01 to 2025-08 (the GraphCast archive starts on 2021-12-31). Training initialisations are in 2022–2023 and their target observation windows must end before 2024-01-01; test initialisations are in 2024-01 to 2025-08.
* `core21`: GFS, Pangu-Weather and FourCastNet v2 from 2021-01; training 2021–2023.
* `ic`: adds the IFS-initialised models, restricted to common initialisations (in-sample only).
* `val23`: historical 2023-only validation; not included in this submission snapshot.

### Expected runtime and disk

Historical measurements on an 8-core machine with 15 GB RAM and a fast connection to AWS us-east-1 (not new measurements for this public snapshot):

| Step | Wall-clock | Disk |
|---|---|---|
| ISD observations 2011–2025 (≈6.6 GB streamed, not stored) | ≈10 min | ≈150 MB (`data/interim`) |
| AI forecasts, GFS-initialised 2021–2025 | ≈3 h | |
| AI forecasts, IFS-initialised 2022–2025 | ≈2.5 h | ≈1.9 GB total (`data/fc`) |
| GFS 2021–2025 (runs in parallel with the AI streaming) | ≈2 h | |
| Analysis cubes (`05_build_dataset.py`) | minutes | ≈1.8 GB (`data/analysis`) |
| Scores and statistics per comparison set | ≈15 min | ≈70 MB per set (`results/scores`) |

The full pipeline takes roughly 7–9 h and about 4 GB of disk. Peak memory stays below a few GB with the default worker counts.

### Not included (rebuilt by the pipeline)

| Path | Size | Regenerate with |
|---|---|---|
| `data/fc/` station-point forecasts | ≈1.9 GB | `run_year.sh`, `03_stream_mlwp.py`, `04_gfs_idx.py` |
| `data/analysis/` analysis cubes and thresholds | ≈1.8 GB | `05_build_dataset.py` (first step of `run_analysis.sh`) |
| `data/interim/isd/`, `data/interim/ghcnd_daily.parquet` | ≈150 MB | `02_obs.py`, `02b_ghcnd.py` |
| `data/raw/` ISD history/inventory files | ≈20 MB | `01_select_stations.py` |
| `results/scores/` per-station and daily scores | ≈240 MB | `07_scores.py` (via `run_analysis.sh`) |

## Building the numerical result tables

```bash
.venv/bin/python paper/make_tables.py
```

The builder refuses results whose provenance does not match the current method hash. It writes `paper/tables/*.tex` from the deposited CSVs. With the manuscript excluded, this is a numerical table-generation step; no manuscript is published or needed. If the author separately supplies private `paper/main_src.tex`, the builder can also inline those tables locally. Regression checks and the observation sample audit are available:

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/02c_precip_audit.py
```

## Known limitations

* ISD global-hourly ends on 2025-08-24, so initialisations end on 2025-08-12.
* The IFS-initialised archive has substantial gaps in 2024.
* No freely available IFS HRES forecasts at station points; GFS was upgraded to v16 on 2021-03-22.
* Tmax/Tmin forecasts use four 6-hourly samples for every model, including GFS. PCRPS removes constant/monotone calibration errors, but not sampling effects that change discrimination.
* ISD missing precipitation amounts require explicit source evidence for zero recovery. Accumulation-boundary amounts require agreement with the original SYNOP 24-hour group. Parser-version sidecars prevent reuse of old caches.
* Precipitation climatology can mix ISD 00–00 UTC and GHCN-Daily 12–12 UTC totals where ISD is sparse.
* Equal-station summary means and pooled station-case bootstrap skill have different weighting; mean station RMSE is not pooled RMSE.
* BH correction is applied across stations separately per model/reference/variable/lead/score, not across every comparison in the study. In-sample uncertainty conditions on fitted IDR maps.
* Among the AI models only GraphCast forecasts precipitation.
* The station set contains a few stations on Taiwan Island in addition to mainland China.

## License

Code: MIT (see `LICENSE`). No new license is asserted for third-party observations or forecasts. Auxiliary map datasets are not bundled: the Köppen–Geiger map is provided by its authors under CC BY 4.0 (Beck et al., 2023); Natural Earth data are in the public domain.
