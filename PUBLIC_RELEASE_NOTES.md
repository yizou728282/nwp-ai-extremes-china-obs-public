# Public delivery version - 2026-10-04

Code and machine-readable derived data are provided through [the public repository](https://github.com/yizou728282/nwp-ai-extremes-china-obs-public) and release `submission-2026-10-04`. `ESM_1.pdf` contains the readable supplementary key-result tables and station inventory.

## Included and excluded material

Included: processing and verification scripts, regression tests and fixtures, station inventories, full derived core/core21/ic result tables, method metadata, available run/provenance records, figure PDFs, generated LaTeX result tables, the table generator, and the supplementary PDF.

The full manuscript, its TeX bodies, its bibliography/template files, and its source ZIP are excluded from both repository contents and release assets. The author retains those materials locally for journal submission. Large forecast-point/observation caches, auxiliary maps, archive chunk indices, analysis cubes and case-level score arrays are also excluded and must be rebuilt from the original public sources. Figure 4's display-only snapshot does not contain effect magnitudes, p values or q values.

## Scientific identity and delivery changes

Baseline review archive SHA-256: `61d6481db12bc9da3d63e40346b924172323dea56431989e97c740efee81aa2e`.
Scientific pipeline: `2026-10-03-verification-v2`.
Method SHA-256: `6ad4b64b41a0eeb72307249c48995698c9c51b26fa89998a8f61cf271a847d9a`.

All 46 scientific result CSVs, 14 method metadata JSON files, station inventories, method-defining scripts, regression tests, available run records and seven verified figure PDFs preserve their reviewed bytes. From the original 190-file archive, 180 files are preserved byte for byte, six manuscript/template files are excluded, and four files are adapted: README.md, REPOSITORY_README.md, paper/make_tables.py and run_corrected.sh.

The table builder now writes numerical LaTeX tables without requiring the private manuscript. The reconstruction driver compiles a manuscript only when its private source and a TeX compiler are present. These conditional delivery steps do not alter scoring, bootstrap, FDR or numerical table expressions. In a separate check, the public table builder completed without manuscript files and regenerated all nine existing table files with identical text; Windows line endings differed. The supplied table files themselves retain their reviewed bytes. The Bash driver's syntax check passed.

The supplementary PDF has 25 pages. All 1,170 result cells match the CSV strings in row order; all 334 station IDs are included. Every supplementary page was rendered and visually reviewed. No scientific scores were recomputed for this packaging change.

## Limits of the evidence

The earlier preflight reported 15 focused regression checks passing. The full synthetic pipeline check was not completed on the Windows review host; no new full source-data download or multi-year numerical rerun was performed for this release. Historical run records do not imply a new rerun. Bash reconstruction commands target a POSIX/Linux environment. Code retains the MIT license; no new license is asserted for third-party observations or forecasts.
