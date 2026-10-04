"""Assemble mmap-able analysis cubes per target variable:
   data/analysis/<var>/{inits,leads,stations}.npy, obs.npy and <model>.npy, each [n_init, n_lead, n_station] float32
Targets (all forecasts sampled identically for every model; GFS uses the same 6-hourly instantaneous fields):
  t2    instantaneous 2 m temperature, leads 6..240 h every 6 h          obs: ISD synoptic t2 at valid time
  ws    10 m wind speed sqrt(u10^2+v10^2), leads 6..240 h                 obs: ISD synoptic wind speed
  tmax  day k=1..9: max of t2 at leads 24k-6, 24k, 24k+6, 24k+12          obs: KA 24 h Tmax, 12-12 UTC day
  tmin  same, min                                                         obs: KA 24 h Tmin
  pr24  day k=1..10: sum of 6-h precipitation over leads 24k-18..24k      obs: AA 24 h total, 00-00 UTC day
        (only GraphCast and GFS forecast precipitation)
Also data/analysis/clim_<var>.parquet: per-station 2011-2020 climatological quantiles (thresholds).
Usage: python 05_build_dataset.py [--models gfs,graphcast,...]
"""
import argparse, numpy as np, pandas as pd
from common import FC, INTERIM, ROOT, STUDY_START, STUDY_END, NLEAD
from provenance import write_metadata
obsmod = __import__("importlib").import_module("02_obs")
AN = ROOT / "data" / "analysis"
MODELS_ALL = ["gfs", "graphcast", "pangu", "fourcastnetv2", "graphcast_ifs", "pangu_ifs", "fourcastnetv2_ifs"]

def load_model(model, inits, sids):
    """Return dict var -> [n_init, 41, n_station] on the common init axis (NaN where missing)."""
    files = sorted((FC / model).glob("*.npz"))
    if not files:
        return None
    out = {}
    ipos = pd.Series(np.arange(len(inits)), index=inits)
    for f in files:
        z = np.load(f, allow_pickle=True)   # own files (station_id saved as object array)
        col = pd.Series(np.arange(len(z["station_id"])), index=z["station_id"]).reindex(sids).values
        assert not np.isnan(col.astype(float)).any()
        rows = ipos.reindex(pd.to_datetime(z["inits"])).values
        ok = ~np.isnan(rows.astype(float))
        for v in z.files:
            if v in ("inits", "station_id"):
                continue
            if v not in out:
                out[v] = np.full((len(inits), NLEAD, len(sids)), np.nan, "f4")
            out[v][rows[ok].astype(int)] = z[v][ok][:, :, col.astype(int)]
    return out

def qc_field(a, lo, hi):
    a = a.copy(); a[(a < lo) | (a > hi)] = np.nan; return a

def main(models):
    missing = [m for m in models if not any((FC / m).glob("*.npz"))]
    if missing:
        raise FileNotFoundError("Missing forecast caches: " + ", ".join(missing))
    st = pd.read_csv(INTERIM / "stations_analysis.csv", dtype={"station_id": str})
    sids = st.station_id.values; ns = len(sids)
    inits = pd.date_range(STUDY_START, STUDY_END, freq="D")
    ni = len(inits)
    # ---------- observations ----------
    syn = obsmod.load("syn", range(2020, 2026), set(sids)).drop_duplicates(["station_id", "valid_time"])
    daily = obsmod.load("daily", range(2020, 2026), set(sids))
    daily = daily.groupby(["station_id", "date"], as_index=False).max()    # year-boundary duplicates
    # Tmax/Tmin: ISD KA groups (12-12 UTC day) filled with GHCN-D TMAX/TMIN (identical day convention and
    # values where both exist; ISD KA coverage drops in 2024-2025)
    gd = pd.read_parquet(INTERIM / "ghcnd_daily.parquet")
    daily = daily.merge(gd[["station_id", "date", "tmax", "tmin"]], on=["station_id", "date"], how="outer", suffixes=("", "_g"))
    daily["tmax"] = daily.tmax.fillna(daily.tmax_g); daily["tmin"] = daily.tmin.fillna(daily.tmin_g)
    scol = pd.Series(np.arange(ns), index=sids)
    def obs_cube(df, tcol, vcol, offsets_h):
        """obs[i, l, s] = value at inits[i] + offsets_h[l]."""
        cube = np.full((ni, len(offsets_h), ns), np.nan, "f4")
        x = df[[tcol, "station_id", vcol]].dropna()
        key = pd.Series(x[vcol].values.astype("f4"), index=pd.MultiIndex.from_arrays([x[tcol].values, scol[x.station_id].values]))
        for l, off in enumerate(offsets_h):
            vt = inits + pd.to_timedelta(off, unit="h")
            idx = pd.MultiIndex.from_product([vt, np.arange(ns)])
            cube[:, l, :] = key.reindex(idx).values.reshape(ni, ns)
        return cube
    leads6 = np.arange(6, 241, 6)
    days9 = np.arange(1, 10); days10 = np.arange(1, 11)
    cubes = {
        "t2": ("t2", obs_cube(syn, "valid_time", "t2", leads6), leads6),
        "ws": ("ws", obs_cube(syn, "valid_time", "ws", leads6), leads6),
        # tmax/tmin on date D = init + k days (window 12 UTC D-1 .. 12 UTC D)
        "tmax": ("tmax", obs_cube(daily, "date", "tmax", 24 * days9), days9),
        "tmin": ("tmin", obs_cube(daily, "date", "tmin", 24 * days9), days9),
        # pr24 on date D = init + (k-1) days (00-00 UTC)
        "pr24": ("pr24", obs_cube(daily, "date", "pr24", 24 * (days10 - 1)), days10),
    }
    for v, (_, o, leads) in cubes.items():
        d = AN / v; d.mkdir(parents=True, exist_ok=True)
        np.save(d / "obs.npy", o); np.save(d / "leads.npy", leads)
        np.save(d / "inits.npy", inits.values.astype("datetime64[s]")); np.save(d / "stations.npy", np.array(list(sids), dtype="U16"))
        print(v, "obs finite", int(np.isfinite(o).sum()))
    # ---------- forecasts ----------
    for m in models:
        f = load_model(m, inits, sids)
        if f is None:
            print(m, "no forecasts yet"); continue
        t2 = qc_field(f["t2"], -70, 60)
        ws = qc_field(np.hypot(f["u10"], f["v10"]), 0, 80)
        np.save(AN / "t2" / f"{m}.npy", t2[:, 1:, :])
        np.save(AN / "ws" / f"{m}.npy", ws[:, 1:, :])
        tx = np.stack([np.max(t2[:, [4*k-1, 4*k, 4*k+1, 4*k+2], :], axis=1) for k in days9], axis=1)  # NaN if any missing
        tn = np.stack([np.min(t2[:, [4*k-1, 4*k, 4*k+1, 4*k+2], :], axis=1) for k in days9], axis=1)
        np.save(AN / "tmax" / f"{m}.npy", tx.astype("f4")); np.save(AN / "tmin" / f"{m}.npy", tn.astype("f4"))
        if "apcp" in f:
            ap = np.clip(qc_field(f["apcp"], -1, 500), 0, None)
            pr = np.stack([np.where(np.isfinite(ap[:, 4*k-3:4*k+1, :]).all(axis=1), np.nansum(ap[:, 4*k-3:4*k+1, :], axis=1), np.nan)
                           for k in days10], axis=1)
            np.save(AN / "pr24" / f"{m}.npy", pr.astype("f4"))
        print(m, {v: int(np.isfinite(np.load(AN / v / f"{m}.npy", mmap_mode="r")).sum()) for v in ("t2", "ws", "tmax")}, flush=True)
    # ---------- climatology 2011-2020 (thresholds; independent of the 2021-2025 evaluation) ----------
    csyn = obsmod.load("syn", range(2011, 2021), set(sids)); cd = obsmod.load("daily", range(2010, 2021), set(sids))
    cd = cd.groupby(["station_id", "date"], as_index=False).max()
    cd = cd[cd.date.dt.year.between(2011, 2020)]
    # ISD KA groups only exist from late 2019 and 24-h AA totals are incomplete before 2019 -> combine with GHCN-D
    g11 = gd[gd.date.dt.year.between(2011, 2020)]
    cd = cd.merge(g11, on=["station_id", "date"], how="outer", suffixes=("", "_g"))
    cd["tmax"] = cd.tmax.fillna(cd.tmax_g); cd["tmin"] = cd.tmin.fillna(cd.tmin_g)
    csyn["hour"] = csyn.valid_time.dt.hour
    rows = []
    for sid in sids:
        a = csyn[csyn.station_id == sid]; b = cd[cd.station_id == sid]
        r = {"station_id": sid, "n_syn": len(a), "n_tmax": int(b.tmax.notna().sum()), "n_tmin": int(b.tmin.notna().sum()),
             "min_month_tmax": int(b[b.tmax.notna()].date.dt.month.value_counts().reindex(range(1, 13), fill_value=0).min()),
             "n_pr_isd": int(b.pr24.notna().sum())}
        pr = b.pr24 if r["n_pr_isd"] >= 1500 else b.pr24.fillna(b.prcp)   # supplement with GHCN-D (12-12 UTC window)
        r["pr_source"] = "isd" if r["n_pr_isd"] >= 1500 else "isd+ghcnd"; r["n_pr"] = int(pr.notna().sum())
        for q in (0.05, 0.10, 0.90, 0.95):
            qq = f"q{int(q*100):02d}"
            for h in (0, 6, 12, 18):          # t2 / ws thresholds are valid-hour specific
                ah = a[a.hour == h]
                r[f"t2_h{h:02d}_{qq}"] = ah.t2.quantile(q); r[f"ws_h{h:02d}_{qq}"] = ah.ws.quantile(q)
            ok = r["n_tmax"] >= 1000 and r["min_month_tmax"] >= 50
            r[f"tmax_{qq}"] = b.tmax.quantile(q) if ok else np.nan; r[f"tmin_{qq}"] = b.tmin.quantile(q) if ok else np.nan
        okp = r["n_pr"] >= 1500
        wet = pr[pr >= 0.1]
        r["pr24_wetfrac"] = (pr[pr.notna()] >= 0.1).mean(); r["pr24_wet_q90"] = wet.quantile(0.9) if okp else np.nan
        r["pr24_wet_q75"] = wet.quantile(0.75) if okp else np.nan
        r["pr24_q95"] = pr.quantile(0.95) if okp else np.nan; r["pr24_q99"] = pr.quantile(0.99) if okp else np.nan
        rows.append(r)
    cl = pd.DataFrame(rows); cl.to_csv(AN / "climatology_2011_2020.csv", index=False)
    # station usage flags per variable: >= 50% of evaluation days observed and thresholds available
    for v, col in (("tmax", "tmax_q90"), ("tmin", "tmin_q10"), ("pr24", "pr24_wet_q90")):
        o = np.load(AN / v / "obs.npy", mmap_mode="r")
        cov = np.isfinite(o[:, 0, :]).mean(axis=0)
        st[f"use_{v}"] = (cov >= 0.5) & np.isfinite(cl[col].values)
    st.to_csv(INTERIM / "stations_analysis.csv", index=False)
    print("climatology written;", {c: int(st[c].sum()) for c in st if c.startswith("use_")})
    write_metadata(AN / "provenance.json", models=models, observation_parser=obsmod.PARSER_VERSION,
                   station_count=len(st), variable_counts={v: int(st[f"use_{v}"].sum()) for v in ("tmax", "tmin", "pr24")})

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--models", default=",".join(MODELS_ALL))
    main(ap.parse_args().models.split(","))
