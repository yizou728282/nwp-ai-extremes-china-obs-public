"""ISD global-hourly station observations, 2011-2025, for the station superset.

Per station-year (resumable: skips only outputs with the current parser version), streamed in memory:
  data/interim/isd/syn/{year}/{sid}.parquet    valid_time, t2 (degC), ws (m/s)  at 00/06/12/18 UTC (+/-10 min, FM-12 preferred)
  data/interim/isd/daily/{year}/{sid}.parquet  date, tmax, tmin (KA 24 h at 12 UTC -> 12-12 UTC day ending on `date`),
                                               pr24 (AA 24 h total reported at 00 UTC of date+1 -> 00-00 UTC day `date`)
Precipitation: AA depth 9999 is missing, never a zero flag. A missing 24-hour AA amount can be
recovered as zero only if an FM-12 report explicitly contains 70000 in SYNOP section 333.
Numeric AA totals carrying accumulation-boundary flags are retained only when the same
unambiguous SYNOP 24-hour group independently confirms their amount.
Unresolved missing amounts, bad/incomplete conditions and conflicting duplicates remain missing.
Per-station/year JSON sidecars retain parser version, source hash and precipitation audit counts.

Usage:  python 02_obs.py [--years 2011-2025] [--threads 6]
        python 02_obs.py --summary        # coverage table + final analysis set (data/interim/stations_analysis.csv)
        python 02_obs.py --check-ghcnd    # precipitation cross-check against GHCN-D
"""
import argparse, hashlib, io, json, re, time, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from common import INTERIM, LOGS, session, log_json, peak_rss_mb

GOOD_QC = set("01459ACIMPRU")
OUT = INTERIM / "isd"
PARSER_VERSION = "isd-aa24-explicit-zero-v2"
# NOAA ISD format: 2 = trace, 9 = condition unavailable, E = estimated.
# Accumulated/deleted/missing-period boundaries and incomplete amounts are excluded.
PR_CONDITIONS = {"2", "9", "E"}

def synop_amount_24(remarks, report_type):
    """Decode an unambiguous FM-12 section-333 total in 0.1 mm (WMO 306, 12.4.9).

    Trace (79999), censored maximum (79998), ambiguous groups and other report
    types are not decoded. The study's 600 mm gross limit is applied here too.
    """
    if str(report_type).strip() != "FM-12" or not isinstance(remarks, str):
        return np.nan
    if not remarks.startswith("SYN"):
        return np.nan
    tokens = remarks.replace("=", " ").split()
    if tokens.count("333") != 1:
        return np.nan
    section = tokens[tokens.index("333") + 1:]
    for i, token in enumerate(section):
        if token in {"444", "555"}:
            section = section[:i]; break
    groups = [t for t in section if re.fullmatch(r"7[0-9/]{4}", t)]
    if len(groups) != 1 or not re.fullmatch(r"7[0-9]{4}", groups[0]):
        return np.nan
    code = int(groups[0][1:])
    return code / 10.0 if code <= 6000 else np.nan

def synop_zero_24(remarks, report_type):
    return synop_amount_24(remarks, report_type) == 0.0

def precipitation(at00):
    """Resolve daily totals from independent, non-conflicting evidence at 00 UTC."""
    candidates = {}
    audit = dict(n_00_reports=len(at00), n_aa24_groups=0, n_missing_aa24_groups=0,
                 n_rejected_condition_groups=0, n_rejected_quality_groups=0,
                 n_zero_recovered_reports=0, n_unresolved_missing_reports=0,
                 n_conflicting_days=0, n_estimated_groups=0, n_flagged_amounts_confirmed_groups=0)
    for _, row in at00.iterrows():
        values = []; recoverable = False; missing = False
        raw_amount = synop_amount_24(row.get("REM"), row.get("REPORT_TYPE"))
        explicit_zero = raw_amount == 0.0
        for key in ("AA1", "AA2", "AA3", "AA4"):
            field = row.get(key)
            if not isinstance(field, str):
                continue
            p = field.split(",")
            if len(p) != 4 or p[0] != "24":
                continue
            audit["n_aa24_groups"] += 1
            if p[1] == "9999":
                missing = True; audit["n_missing_aa24_groups"] += 1
            if p[2] in {"3", "4"} and p[1].isdigit() and p[1] != "9999" and p[3] in GOOD_QC:
                # Some Chinese records attach a boundary flag to a valid SYNOP 24-hour amount.
                # Preserve only amounts confirmed by the independent source group, not the flag.
                value = int(p[1]) / 10.0
                if np.isfinite(raw_amount) and np.isclose(value, raw_amount, rtol=0, atol=1e-8):
                    values.append((value, "synop333-confirmed-aa24"))
                    audit["n_flagged_amounts_confirmed_groups"] += 1
                    continue
            if p[2] not in PR_CONDITIONS:
                audit["n_rejected_condition_groups"] += 1; continue
            if p[3] not in GOOD_QC:
                audit["n_rejected_quality_groups"] += 1; continue
            if p[1] == "9999":
                recoverable |= p[2] == "9"
                continue
            try:
                value = int(p[1]) / 10.0
            except ValueError:
                continue
            if 0 <= value <= 600:
                values.append((value, "aa24"))
                audit["n_estimated_groups"] += int(p[2] == "E")
        if recoverable and explicit_zero:
            values.append((0.0, "synop33370000"))
            audit["n_zero_recovered_reports"] += 1
        if missing and not values:
            audit["n_unresolved_missing_reports"] += 1
        day = row["time"].normalize() - pd.Timedelta(days=1)
        candidates.setdefault(day, []).extend(values)
    pr, sources = {}, {}
    for day, values in candidates.items():
        distinct = {v for v, _ in values}
        if len(distinct) == 1:
            pr[day] = distinct.pop()
            sources[day] = "+".join(sorted({s for _, s in values}))
        elif len(distinct) > 1:
            audit["n_conflicting_days"] += 1
    audit["n_valid_pr_days"] = len(pr)
    audit["n_zero_recovered_days"] = sum("synop33370000" in s for s in sources.values())
    return pr, sources, audit

def qc_num(val, q, scale, missing):
    x = pd.to_numeric(val, errors="coerce") / scale
    x[(val == missing) | ~q.isin(GOOD_QC)] = np.nan
    return x

def parse(text, sid, year, return_audit=False):
    cols = ("DATE", "REPORT_TYPE", "TMP", "WND", "KA1", "KA2", "KA3", "KA4", "AA1", "AA2", "AA3", "AA4", "REM")
    d = pd.read_csv(io.StringIO(text), dtype=str, usecols=lambda c: c in cols)
    d["time"] = pd.to_datetime(d["DATE"])
    tmp = d["TMP"].fillna("+9999,9")
    d["t2"] = qc_num(tmp.str.slice(0, 5), tmp.str.slice(6, 7), 10.0, "+9999")
    d.loc[(d.t2 < -60) | (d.t2 > 55), "t2"] = np.nan
    w = d["WND"].fillna("999,9,9,9999,9").str.split(",", expand=True)
    d["ws"] = qc_num(w[3], w[4], 10.0, "9999")
    d.loc[w[2] == "C", "ws"] = 0.0
    d.loc[(d.ws < 0) | (d.ws > 60), "ws"] = np.nan
    # synoptic-hour values
    d["vt"] = d["time"].dt.round("6h")
    d["dt"] = (d["time"] - d["vt"]).abs()
    s = d[(d["dt"] <= pd.Timedelta("10min")) & (d.t2.notna() | d.ws.notna())].copy()
    s["pref"] = (s["REPORT_TYPE"].str.strip() != "FM-12").astype(int)
    s = s.sort_values(["vt", "pref", "dt"])
    syn = pd.DataFrame({"t2": s.groupby("vt")["t2"].first(), "ws": s.groupby("vt")["ws"].first()}).reset_index().rename(columns={"vt": "valid_time"})
    syn = syn[syn.valid_time.dt.year == year].astype({"t2": "float32", "ws": "float32"})
    # daily extremes (KA, 12 UTC, 24 h)
    at12 = d[(d.time.dt.hour == 12) & (d.time.dt.minute == 0)]
    tx, tn = {}, {}
    for k in ("KA1", "KA2", "KA3", "KA4"):
        if k not in at12:
            continue
        g = at12[k].dropna()
        if g.empty:
            continue
        p = g.str.split(",", expand=True)
        if p.shape[1] < 4:
            continue
        t = qc_num(p[2], p[3], 10.0, "+9999")
        ok = (pd.to_numeric(p[0], errors="coerce") == 240) & t.between(-60, 55)
        dates = at12.loc[g.index, "time"].dt.normalize()
        for i in np.where(ok)[0]:
            dd = dates.iloc[i]
            if p[1].iloc[i] == "M":
                tx[dd] = max(tx.get(dd, -99), t.iloc[i])
            elif p[1].iloc[i] == "N":
                tn[dd] = min(tn.get(dd, 99), t.iloc[i])
    # daily precipitation (AA, 00 UTC, 24 h) -> previous day
    at00 = d[(d.time.dt.hour == 0) & (d.time.dt.minute == 0)]
    pr, sources, audit = precipitation(at00)
    idx = sorted(set(tx) | set(tn) | set(pr))
    daily = pd.DataFrame({"date": pd.to_datetime(idx)})
    daily["tmax"] = np.array([tx.get(i, np.nan) for i in idx], dtype="float32")
    daily["tmin"] = np.array([tn.get(i, np.nan) for i in idx], dtype="float32")
    daily["pr24"] = np.array([pr.get(i, np.nan) for i in idx], dtype="float32")
    daily["pr24_source"] = [sources.get(i, "missing") for i in idx]
    return (syn, daily, audit) if return_audit else (syn, daily)

def one(sid, year):
    fs, fd = OUT / "syn" / str(year) / f"{sid}.parquet", OUT / "daily" / str(year) / f"{sid}.parquet"
    fm = OUT / "audit" / str(year) / f"{sid}.json"
    if fs.exists() and fd.exists() and fm.exists():
        try:
            if json.loads(fm.read_text())["parser_version"] == PARSER_VERSION:
                return None
        except (ValueError, KeyError):
            pass
    t0 = time.time()
    url = f"https://www.ncei.noaa.gov/data/global-hourly/access/{year}/{sid}.csv"
    for k in range(4):
        try:
            r = session().get(url, timeout=300)
            break
        except Exception:
            time.sleep(5 * (k + 1))
    else:
        return {"sid": sid, "year": year, "status": "error"}
    if r.status_code == 404:
        syn = pd.DataFrame({"valid_time": pd.Series(dtype="datetime64[ns]"), "t2": pd.Series(dtype="float32"), "ws": pd.Series(dtype="float32")})
        daily = pd.DataFrame({"date": pd.Series(dtype="datetime64[ns]"), "tmax": pd.Series(dtype="float32"), "tmin": pd.Series(dtype="float32"), "pr24": pd.Series(dtype="float32")})
        nb = 0
        audit = {}
    else:
        r.raise_for_status()
        nb = len(r.content)
        syn, daily, audit = parse(r.text, sid, year, return_audit=True)
    fs.parent.mkdir(parents=True, exist_ok=True); fd.parent.mkdir(parents=True, exist_ok=True)
    syn.to_parquet(fs, index=False); daily.to_parquet(fd, index=False)
    meta = {"sid": sid, "year": year, "status": "ok" if nb else "404", "MB": round(nb / 1e6, 2), "sec": round(time.time() - t0, 1),
            "parser_version": PARSER_VERSION, "source_url": url, "source_sha256": hashlib.sha256(r.content).hexdigest(),
            "n_syn": len(syn), "n_tmax": int(daily.tmax.notna().sum()), "n_pr": int(daily.pr24.notna().sum()), **audit}
    fm.parent.mkdir(parents=True, exist_ok=True)
    fm.write_text(json.dumps(meta, indent=2) + "\n")
    return meta

def load(kind, years, sids=None):
    fr = []
    for y in years:
        for f in sorted((OUT / kind / str(y)).glob("*.parquet")):
            if sids is None or f.stem in sids:
                audit = OUT / "audit" / str(y) / f"{f.stem}.json"
                if not audit.exists() or json.loads(audit.read_text()).get("parser_version") != PARSER_VERSION:
                    raise RuntimeError(f"Legacy ISD cache: {f}; rerun 02_obs.py --years 2011-2025")
                x = pd.read_parquet(f); x.insert(0, "station_id", f.stem); fr.append(x)
    if not fr:
        raise FileNotFoundError(f"No {kind} observations for the requested stations/years")
    return pd.concat(fr, ignore_index=True)

def summary():
    sup = pd.read_csv(INTERIM / "stations_superset.csv", dtype={"station_id": str})
    cand = set(pd.read_csv(INTERIM / "stations.csv", dtype={"station_id": str}).station_id)
    syn = load("syn", range(2011, 2026)); daily = load("daily", range(2011, 2026))
    s00 = syn[(syn.valid_time.dt.hour == 0) & syn.t2.notna()]
    s00 = s00.assign(year=s00.valid_time.dt.year)
    cnt = s00.groupby(["station_id", "year"]).size().unstack(fill_value=0)
    ev_days = (pd.Timestamp("2025-08-22") - pd.Timestamp("2021-01-01")).days + 1
    ev = s00[s00.valid_time.between("2021-01-01", "2025-08-22")].groupby("station_id").size() / ev_days
    clim_ok = (cnt.reindex(columns=range(2011, 2021), fill_value=0) >= 0.8 * 365).sum(axis=1)
    dd = daily.assign(year=daily.date.dt.year)
    tx_ev = dd[dd.date.between("2021-01-01", "2025-08-22")].groupby("station_id").tmax.count() / ev_days
    pr_ev = dd[dd.date.between("2021-01-01", "2025-08-22")].groupby("station_id").pr24.count() / ev_days
    pr_cl = dd[dd.year.between(2011, 2020)].groupby("station_id").pr24.count() / 3653
    tx_cl = dd[dd.year.between(2011, 2020)].groupby("station_id").tmax.count() / 3653
    cov = pd.DataFrame({"frac00_eval": ev, "clim_years_ok": clim_ok, "tmax_eval": tx_ev, "tmax_clim": tx_cl, "pr_eval": pr_ev, "pr_clim": pr_cl}).reindex(sup.station_id).fillna(0)
    cov["candidate"] = cov.index.isin(cand)
    cov["analysis"] = cov.candidate & (cov.frac00_eval >= 0.8) & (cov.clim_years_ok >= 8)
    cov["use_tmax"] = cov.analysis & (cov.tmax_eval >= 0.7) & (cov.tmax_clim >= 0.7)
    cov["use_pr"] = cov.analysis & (cov.pr_eval >= 0.7) & (cov.pr_clim >= 0.7)
    cov.round(3).to_csv(INTERIM / "station_coverage.csv")
    st = sup.set_index("station_id").join(cov[["analysis", "use_tmax", "use_pr"]])
    st = st[st.analysis].reset_index()
    st.to_csv(INTERIM / "stations_analysis.csv", index=False)
    out = {"superset": len(sup), "candidates": len(cand), "analysis": int(cov.analysis.sum()), "use_tmax": int(cov.use_tmax.sum()),
           "use_pr": int(cov.use_pr.sum()), "zones": st.zone.value_counts().to_dict()}
    log_json("02_summary.json", out); print(out)

def check_ghcnd(n=40):
    st = pd.read_csv(INTERIM / "stations_analysis.csv", dtype={"station_id": str})
    rows = []
    for sid in st.station_id[:: max(1, len(st) // n)]:
        wmo = sid[:5]
        r = session().get(f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/CHM000{wmo}.csv", timeout=120)
        if r.status_code != 200:
            continue
        g = pd.read_csv(io.StringIO(r.text), dtype=str)
        g = g[g.ELEMENT == "PRCP"]
        g = g[g.Q_FLAG.fillna("").str.strip() == ""]
        g = pd.DataFrame({"date": pd.to_datetime(g.DATE), "ghcnd": pd.to_numeric(g.DATA_VALUE) / 10.0})
        g = g[g.ghcnd.between(0, 600)]
        g = g[g.date.between("2015-01-01", "2024-12-31")]
        isd = load("daily", range(2015, 2025), {sid})[["date", "pr24"]]
        m = g.merge(isd, on="date").dropna()
        if len(m) < 300:
            continue
        rows.append({"sid": sid, "n": len(m), "corr": m.ghcnd.corr(m.pr24), "bias_mm": (m.pr24 - m.ghcnd).mean(),
                     "wet_agree": ((m.ghcnd >= 0.1) == (m.pr24 >= 0.1)).mean(), "mean_ghcnd": m.ghcnd.mean(), "mean_isd": m.pr24.mean()})
    df = pd.DataFrame(rows); df.round(3).to_csv(LOGS / "02_ghcnd_precip_check.csv", index=False)
    print(df.describe().round(3))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", default="2011-2025"); ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--summary", action="store_true"); ap.add_argument("--check-ghcnd", action="store_true")
    ap.add_argument("--stations", default="superset")
    a = ap.parse_args()
    if a.summary:
        summary(); raise SystemExit
    if a.check_ghcnd:
        check_ghcnd(); raise SystemExit
    y0, y1 = map(int, a.years.split("-")) if "-" in a.years else (int(a.years), int(a.years))
    st = pd.read_csv(INTERIM / ("stations_superset.csv" if a.stations == "superset" else "stations.csv"), dtype={"station_id": str})
    jobs = [(s, y) for y in range(y1, y0 - 1, -1) for s in st.station_id]
    t0 = time.time(); meta = []; done = 0
    with ThreadPoolExecutor(a.threads) as ex:
        futs = [ex.submit(one, s, y) for s, y in jobs]
        for f in as_completed(futs):
            try:
                m = f.result()
            except Exception as e:
                m = {"status": "error", "err": repr(e)[:200]}
            done += 1
            if m:
                meta.append(m)
            if done % 200 == 0:
                print(f"{done}/{len(jobs)} {time.time() - t0:.0f}s rss={peak_rss_mb()}MB", flush=True)
    meta = pd.DataFrame(meta)
    if len(meta):
        meta.to_csv(LOGS / f"02_obs_{y0}_{y1}_{int(time.time())}.csv", index=False)
        print(meta.groupby("status").size().to_dict(), "MB", round(meta.get("MB", pd.Series([0])).sum(), 1), "wall", round(time.time() - t0), "s")
        if (meta.status == "error").any():
            raise SystemExit("Observation downloads failed; rerun to resume before analysis.")
