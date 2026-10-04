# Real ISD regression fixture

`isd_beijing_2023.csv` preserves the DATE, REPORT_TYPE, TMP, WND, AA1--AA3 and REM fields of the FM-12 reports at 00 UTC on 2023-04-06, 2023-04-21 and 2023-05-04 for ISD station 54511099999. Other fields and reports have been omitted.

Source: https://www.ncei.noaa.gov/data/global-hourly/access/2023/54511099999.csv

Each report contains `24,9999,9,9` and an explicit `70000` in SYNOP section 333. The fixture verifies source-backed zero recovery and the previous-day label; it does not imply that other missing amounts are zero.
