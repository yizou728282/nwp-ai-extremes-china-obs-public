"""Audit the legacy missing-as-zero precipitation rule on reproducibly selected ISD files.

This validates parsing on a sample; it does not recompute forecast verification scores.
Usage: python scripts/02c_precip_audit.py [--years 2011,2023,2025] [--per-zone 2]
"""
import argparse
import hashlib
import importlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from common import ROOT, RESULTS, session

obs = importlib.import_module("02_obs")

def audit_one(job):
    sid, year, zone = job
    url = f"https://www.ncei.noaa.gov/data/global-hourly/access/{year}/{sid}.csv"
    r = session().get(url, timeout=300)
    row = dict(station_id=sid, year=year, zone=zone, source_url=url,
               parser_version=obs.PARSER_VERSION, status=r.status_code)
    if r.status_code == 404:
        return row
    r.raise_for_status()
    row["source_sha256"] = hashlib.sha256(r.content).hexdigest()
    d = pd.read_csv(io.StringIO(r.text), dtype=str)
    d["time"] = pd.to_datetime(d.DATE)
    at00 = d[(d.time.dt.hour == 0) & (d.time.dt.minute == 0)]
    new, _, counts = obs.precipitation(at00)
    # The previous rule is reproduced only here to measure the change, not used for verification.
    old = {}
    for _, record in at00.iterrows():
        for key in ("AA1", "AA2", "AA3", "AA4"):
            field = record.get(key)
            if not isinstance(field, str):
                continue
            p = field.split(",")
            if len(p) != 4 or p[0] != "24" or p[3] not in obs.GOOD_QC:
                continue
            try:
                value = 0.0 if p[1] == "9999" else int(p[1]) / 10
            except ValueError:
                continue
            if 0 <= value <= 600:
                day = record.time.normalize() - pd.Timedelta(days=1)
                old[day] = max(old.get(day, 0), value)
    row.update(counts)
    row["legacy_valid_days"] = len(old)
    row["legacy_zero_days"] = sum(v == 0 for v in old.values())
    row["legacy_zero_days_retained"] = sum(v == 0 and new.get(day) == 0 for day, v in old.items())
    row["legacy_zero_days_removed"] = sum(v == 0 and day not in new for day, v in old.items())
    row["legacy_positive_days_removed"] = sum(v > 0 and day not in new for day, v in old.items())
    row["changed_numeric_days"] = sum(day in new and not np.isclose(v, new[day]) for day, v in old.items())
    return row

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", default="2011,2023,2025")
    ap.add_argument("--per-zone", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    st = pd.read_csv(ROOT / "data/interim/stations_analysis.csv", dtype={"station_id": str})
    sample = st.sort_values(["zone", "station_id"]).groupby("zone", sort=True).head(a.per_zone)
    jobs = [(sid, int(y), zone) for sid, zone in sample[["station_id", "zone"]].itertuples(index=False, name=None)
            for y in a.years.split(",")]
    with ThreadPoolExecutor(a.threads) as ex:
        rows = list(ex.map(audit_one, jobs))
    out = RESULTS / "validation"; out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "precip_parser_audit.csv", index=False)
    totals = pd.DataFrame(rows).select_dtypes(include="number").drop(columns=["year", "status"]).sum().astype(int).to_dict()
    summary = dict(parser_version=obs.PARSER_VERSION, station_years=len(jobs),
                   stations=sample.station_id.tolist(), selection="first station IDs in each climate group",
                   scope="parser validation sample, not a forecast-score rerun", totals=totals)
    (out / "precip_parser_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
