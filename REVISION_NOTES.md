# Scientific correction, 2026-10-03

Base repository commit: `701515adb272b1034ab19bae9c811212043d0350`.

The parser, statistics and manuscript have been corrected, and the full verification has been recomputed with the corrected pipeline from the existing forecast-point caches (`bash run_corrected.sh`). The repository still excludes `data/fc`, `data/analysis`, observation caches and `results/scores`; the tracked tables and figures were regenerated from them. All quantitative statements in the manuscript were reviewed against the regenerated tables, and `results/revision_status.json` records the reviewed method hash.

## Corrections

| Issue | Change | Required regeneration |
|---|---|---|
| Bootstrap discarded days with zero GFS loss | Retain valid zero losses; use availability counts; preserve missing dates on the calendar axis | All bootstrap intervals and their significance flags |
| ISD precipitation depth `9999` treated as zero | Missing amounts remain missing unless FM-12 section 333 explicitly contains `70000`; flag-3/4 numeric totals require agreement with the original 24-hour SYNOP group | Observation caches, precipitation eligibility, thresholds, scores and inference |
| Future training labels available only after the first test initialization | Purge training targets whose observation windows end on/after 2024-01-01; account for the extra 12 h in Tmax/Tmin targets | All out-of-sample scores |
| HAC removed missing calendar days | Preserve date spacing in loss differences | All station DM/FDR results |
| Instantaneous inference saved only daily leads | Save all six-hourly leads and add valid-hour panels; retain daily-lead EMOS | Six-hourly score/inference output and diagnostic figure |
| Heavy-rain conclusions inferred from overall OOS scores | Add the complete wet-day-q90 OOS comparison; distinguish overall and tail skill, and point rankings from paired significance | Recompute the new table and review conclusions |
| “Identical analyses” and “pure model effect” overclaimed | Specify GraphCast Operational, FourCastNet v2-small, GFS F000 products and conditional initialization sensitivity | Manuscript and documentation corrected |
| Figure descriptions and climate counts mismatched actual plots | Correct row/column captions, day-5 heatmap labels, variable-specific counts, significance markers and IC interval shading | Figures regenerated from recomputed results |
| Silent reuse of legacy outputs | Observation version sidecars, method hashes for cubes/scores/tables, strict paper-build gate and fail-fast rerun driver | Rebuild with `bash run_corrected.sh` |
| Gaussian threshold-weighted CRPS used numerical integration that reversed direction when the threshold lay beyond the effective forecast support, giving negative or undefined EMOS tail scores | Closed-form twCRPS for the Gaussian predictive distribution (upper and lower tails, finite intervals by difference), clipped at zero | All EMOS tail scores (`oos_emos_tw1`) |
| Rerun cache check required GraphCast monthly caches from 2021-01, although the archive has no 00 UTC GraphCast runs in 2021 | Require GraphCast caches from 2022-01, consistent with the *core* set | Rerun driver only |

## Executed validation

- Sixteen regression tests pass, covering valid zero losses, calendar gaps, calibration timing, real and synthetic precipitation reports, conflicting records, model availability and CRPS and Gaussian threshold-weighted CRPS (including thresholds far outside the forecast support) against independent numerical quadrature. A synthetic paired-cube integration test exercises scoring, OOS purging, both tail scores, aggregation and provenance end to end.
- A reproducible parser audit downloads two station IDs per climate group for 2011, 2023 and 2025: ten stations and thirty station-years. Source URLs and SHA-256 hashes are recorded in `results/validation/precip_parser_audit.csv`.
- In that sample, 316 legacy zero days have no valid zero evidence and are removed; 44 zero days are recovered from explicit `70000`. Another 575 accumulation-flagged numeric groups are retained because the source 24-hour group confirms their amount. These counts describe the sample only; they are not a national correction estimate.
- The audit also removes 416 legacy positive days lacking acceptable complete-period evidence. This can change station completeness and tail quantiles, so precipitation skill cannot be updated from the old summary CSVs alone.
- The Beijing fixture contains three authentic 2023 FM-12 reports with missing AA totals and explicit section-333 zero amounts.
- The default paper-build and rerun commands reject missing/stale input caches rather than silently mixing old and new results. A labelled historical draft can be built explicitly with `python paper/make_tables.py --allow-legacy`.

## Recompute

```bash
# Restore or download all forecast-point caches using README's source commands.
bash run_corrected.sh
```

The driver refreshes ISD caches whose parser version is old, runs GHCN-Daily checks, rebuilds analysis cubes and thresholds, recalculates core/core21/ic scores and inference, copies regenerated figures into the paper, and builds all tables and the PDF. The full archive download is a separate several-hour step; this driver uses existing forecast caches. All six-hourly leads increase the output volume compared with the original daily-lead inference files.

The paper shows a prose-review notice whenever `manuscript_numeric_review_method_sha256` in `results/revision_status.json` differs from the current method hash (`provenance.method_hash()`). For the current results the quantitative prose, climate-group counts and significance statements were updated from the regenerated tables (precipitation now uses 249 stations; the out-of-sample wet-day-q90 tail skill is significant only at day 1), and the field was set to the reviewed hash. Any later method change requires a new rerun and prose review.

## Source basis

- NOAA ISD format documentation: https://www.ncei.noaa.gov/pub/data/noaa/isd-format-document.pdf (AA amount/condition/quality fields).
- WMO Manual on Codes, No. 306, regulation 12.4.9: https://mgm.gov.tr/FILES/resmi-istatistikler/WMO_NO_306.pdf (SYNOP 24-hour totals in tenths of millimetres; distinguishes `70000`, trace and missing groups).
- NOAA/CIRA AIWP archive: https://registry.opendata.aws/aiwp/ and Radford et al. (2025), https://doi.org/10.1175/BAMS-D-24-0057.1.
- Related Chinese station precipitation study: Yan et al. (2025), https://doi.org/10.1038/s41598-025-98944-7.
