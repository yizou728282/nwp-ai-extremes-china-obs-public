# Code and derived data for the China station verification study

**Article:** Do AI weather models forecast extremes better than GFS? A station-based verification over China with a common initialization source, 2022-2025  
**Author:** Yi Zou  
**Prepared for:** Asia-Pacific Journal of Atmospheric Sciences (initial submission)  
**Scientific pipeline:** 2026-10-03-verification-v2  
**Method SHA-256:** `6ad4b64b41a0eeb72307249c48995698c9c51b26fa89998a8f61cf271a847d9a`

## Obtain the submission version

- [Dated release: submission-2026-10-04](https://github.com/yizou728282/nwp-ai-extremes-china-obs-public/releases/tag/submission-2026-10-04)
- [Download all code and derived data](https://github.com/yizou728282/nwp-ai-extremes-china-obs-public/releases/download/submission-2026-10-04/Code_and_derived_data.zip)
- [Readable supplementary tables and station inventory](ESM_1.pdf)
- [Detailed data sources and reconstruction instructions](REPOSITORY_README.md)
- [Per-file SHA-256 inventory](PUBLICATION_MANIFEST.json)

This public snapshot contains processing and verification code, regression tests and fixtures, station lists, full derived core/core21/ic result tables, method metadata, available run/provenance records, and figure sources. `ESM_1.pdf` is the readable supplementary review file; the CSV/JSON results in this repository retain the detailed machine-readable values and uncertainty estimates. Manuscript PDFs, manuscript TeX bodies and the manuscript source ZIP are retained locally and excluded from this public repository and its release assets. The dated release identifies the version corresponding to this submission. This is a submission snapshot, not a claim of journal acceptance.

It excludes large forecast-point and observation caches, auxiliary map datasets, archive chunk indices, analysis cubes and case-level score arrays. A full numerical rerun requires reconstructing those inputs from the public sources using the supplied scripts. The historical `val23` analysis is not included. Station eligibility for the reported analysis is recorded in `data/interim/stations_analysis.csv`; some run logs contain preliminary counts from earlier processing stages. The Figure 4 snapshot is display-only and does not contain effect magnitudes, p values or q values.

The archived scientific results were not recomputed when preparing this public delivery version. The method-defining verification code, result CSVs, method metadata and seven figure PDFs are preserved byte for byte. The table builder and reconstruction driver were adjusted only to make private manuscript compilation optional; their numerical computations are unchanged. See `PUBLIC_RELEASE_NOTES.md` for the scope and validation.

## Environment and regression checks

Use Python 3.12 and requirements.txt. Create .venv and install those pinned dependencies. Run:

    .venv/bin/python -m unittest discover -s tests -v

The test suite checks precipitation parsing, calibration cutoffs, analytical Gaussian threshold scores, calendar-gap handling and provenance guards. It does not independently validate the complete forecast archive.

## Rebuild workflow

First rebuild the auxiliary data and public observation/forecast inputs following REPOSITORY_README.md. In particular, run 01_select_stations.py to fetch auxiliary datasets, 02_obs.py and 02b_ghcnd.py for observations, 03_refs.py for AIWP chunk indices, and the streaming drivers for forecast-point caches. Create logs/ before running the download drivers. Then run:

    bash run_corrected.sh

The corrected driver refuses missing forecast-point caches and stale method metadata. It rebuilds observations, thresholds, cubes, scores, bootstrap/FDR tables and paper tables. The public table builder writes the LaTeX result tables without requiring manuscript source. Manuscript compilation is optional and requires separately supplied private manuscript/template files and a TeX distribution. Some download URLs or dependencies may change in the future.

## Contents and interpretation

core: 2022-2025, GFS plus three GFS-initialized AI systems.
core21: adds 2021, omitting GraphCast.
ic: paired GFS/IFS initialization on common dates.
summary CSV: equal-station statistics; bootstrap skill: pooled station-case score sums.
Metadata JSON: method hash and input/provenance records.
results/revision_status.json: numerical prose checked against the current method hash.

The code is covered by the accompanying repository MIT LICENSE. No new license is asserted for third-party observations or forecasts.


## Publication figure layouts

Run `python scripts/09_publication_figures.py` after obtaining the auxiliary coastline data to redraw the seven submission figures from the unchanged summary tables. Model markers, line styles and signed labels provide monochrome identification.


The plotting script automatically creates `paper/submission/apjas/` before
writing its `Figure_source_checks.json` record.

The supplied Figure 4 snapshot can be used directly when redrawing the figures.
To regenerate that display-only snapshot from the archived original vector PDF,
install the optional extraction tools in addition to `requirements.txt`, then run
the extraction script from this archive's root:

    .venv/bin/python -m pip install -r requirements-figure-extraction.txt
    .venv/bin/python scripts/extract_figure4_snapshot.py

`requirements-figure-extraction.txt` lists the tested `pypdf` and `pdfplumber`
versions. These tools are needed only for extracting the snapshot; they are not
needed to use the snapshot already included here.

Figure 4 uses `data/interim/fig4_plot_snapshot.csv`, a display-only extraction from the archived original vector figure `paper/figure_sources/fig4_original.pdf`. `scripts/extract_figure4_snapshot.py` documents the extraction. All 2,289 marker positions match existing station identities; significance/sign counts agree with all 18 reported better/worse fractions at their four-decimal CSV precision. Original RGB fills are retained. The snapshot supplies no numerical effect magnitudes, p values or q values and is not a substitute for case-level verification data. Its provenance JSON records source and snapshot hashes.
