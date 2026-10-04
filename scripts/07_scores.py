"""Probabilistic verification of deterministic forecasts at stations.

For each (station, lead, model) on a common sample (all models of the comparison set + obs available):
  * raw deterministic: bias, MAE, RMSE
  * in-sample potential CRPS (PCRPS, Gneiting et al. 2026): EasyUQ/IDR of obs on the forecast fitted on the
    whole evaluation sample; CRPS of the fitted step CDFs.  Threshold-weighted PCRPS (twPCRPS) with
    indicator weights at the station's 2011-2020 climatological quantiles (independent of 2021-2025).
  * out-of-sample: IDR (CDFs at new x linearly interpolated between neighbouring fitted CDFs) and, for
    temperature variables, Gaussian EMOS with seasonal terms (mu = a + b x + c sin + d cos; log sigma =
    e + f sin + g cos; minimum-CRPS fit), trained on TRAIN period, scored on TEST period.
  * reference: empirical climatology of the same (in-sample) or training (out-of-sample) obs.
  * Diebold-Mariano test per station (AI vs reference model; Newey-West HAC variance) on per-day losses.
  * per-day station sums of each score (by climate zone) for the date-block bootstrap in 08_stats.py.
Comparison sets:
  core    gfs, graphcast, pangu, fourcastnetv2 ; inits 2022-01-01..2025-08-12 ; train 2022-2023, test 2024-01..2025-08
  core21  gfs, pangu, fourcastnetv2 ; inits 2021-01-01..2025-08-12 ; train 2021-2023 (sensitivity)
  ic      gfs + 3 AI GFS-init + 3 AI IFS-init ; common inits (in-sample only; IFS-init 2024 has gaps)
Usage: python 07_scores.py --var t2 --set core [--workers 4]
"""
import os
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")   # one BLAS thread per worker (shared box)
import time, argparse, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import isotonic_regression, minimize
try:   # direct call of scipy's PAVA kernel (same algorithm as isotonic_regression, without per-call overhead)
    from scipy.optimize._pava_pybind import pava as _pava
except ImportError:
    _pava = None
from scipy.stats import norm
from common import ROOT, RESULTS, log_json, peak_rss_mb
from provenance import require_current, write_metadata

AN = ROOT / "data" / "analysis"
OUT = RESULTS / "scores"; OUT.mkdir(parents=True, exist_ok=True)
SETS = {
    "core": dict(models=["gfs", "graphcast", "pangu", "fourcastnetv2"], start="2022-01-01", train_end="2023-12-31", ref="gfs"),
    "core21": dict(models=["gfs", "pangu", "fourcastnetv2"], start="2021-01-01", train_end="2023-12-31", ref="gfs"),
    # pipeline validation on 2023 only (train Jan-Jun, test Jul-Dec; OOS numbers are NOT meaningful)
    "val23": dict(models=["gfs", "graphcast", "pangu", "fourcastnetv2"], start="2023-01-01", train_end="2023-06-30", ref="gfs", end="2023-12-31"),
    "ic": dict(models=["gfs", "graphcast", "pangu", "fourcastnetv2", "graphcast_ifs", "pangu_ifs", "fourcastnetv2_ifs"],
               start="2022-01-01", train_end=None, ref="gfs"),
}
# thresholds: (name, climatology column template, tail)
THRESH = {
    "t2": [("hi95", "t2_h{h:02d}_q95", "upper"), ("lo05", "t2_h{h:02d}_q05", "lower")],
    "ws": [("hi90", "ws_h{h:02d}_q90", "upper"), ("hi95", "ws_h{h:02d}_q95", "upper")],
    "tmax": [("hi90", "tmax_q90", "upper"), ("hi95", "tmax_q95", "upper")],
    "tmin": [("lo10", "tmin_q10", "lower"), ("lo05", "tmin_q05", "lower")],
    "pr24": [("wet90", "pr24_wet_q90", "upper"), ("ge25", 25.0, "upper")],
}
EMOS_VARS = {"t2", "tmax", "tmin"}
SAVE_SCORES = ["pcrps", "tw1", "tw2", "oos_idr", "oos_tw1", "oos_emos", "oos_emos_tw1", "abs_err"]

# ---------------------------------------------------------------- EasyUQ / IDR
def idr_fit(x, y):
    """IDR with one real covariate. Returns ux (sorted unique x), z (sorted unique y), F[len(ux), len(z)], inv (x->ux index)."""
    z = np.unique(y)
    ux, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    # cumulative indicator sums per unique x for all z at once: sort y, count per ux
    yr = np.searchsorted(z, y)                                  # rank of each y in z
    H = np.zeros((len(ux), len(z)))
    np.add.at(H, (inv, yr), 1.0)
    H = np.cumsum(H, axis=1) / cnt[:, None]                      # empirical P(Y<=z_k | x = ux_i)
    F = np.empty_like(H)
    n = len(ux); wrev = np.ascontiguousarray(cnt[::-1], dtype=np.float64)
    Hr = np.ascontiguousarray(H[::-1].T)                         # decreasing in x == increasing in reversed order
    for k in range(len(z) - 1):
        if _pava is not None:
            xk = Hr[k].copy(); _pava(xk, wrev.copy(), np.full(n + 1, -1, dtype=np.intp)); F[:, k] = xk[::-1]
        else:
            F[:, k] = isotonic_regression(H[:, k], weights=cnt, increasing=False).x
    F[:, -1] = 1.0
    return ux, z, np.clip(F, 0, 1), inv

def idr_predict(ux, F, xnew):
    """Linear interpolation between fitted CDFs of neighbouring training x; constant beyond the range."""
    j = np.searchsorted(ux, xnew)
    j0 = np.clip(j - 1, 0, len(ux) - 1); j1 = np.clip(j, 0, len(ux) - 1)
    d = ux[j1] - ux[j0]
    w = np.where(d > 0, (xnew - ux[j0]) / np.where(d > 0, d, 1), 0.0)
    return F[j0] * (1 - w[:, None]) + F[j1] * w[:, None]

def crps_step(z, F, y, a=-np.inf, b=np.inf):
    """CRPS of step CDFs (jumps at z, F[:,k] on [z_k, z_k+1)) restricted to the weight region [a, b]
    (indicator-weighted twCRPS). Handles y outside the support (out-of-sample)."""
    if len(z) == 1:
        return np.abs(np.clip(y, a, b) - np.clip(z[0], a, b))
    lo = np.clip(z[:-1], a, b); hi = np.clip(z[1:], a, b)
    w = hi - lo
    I = (y[:, None] <= z[None, :-1])
    s = ((F[:, :-1] - I) ** 2 * w[None, :]).sum(axis=1)
    # y strictly inside an interval (z_k, z_k+1): indicator switches to 1 at y within that interval
    k = np.searchsorted(z, y, side="right") - 1
    inside = (k >= 0) & (k < len(z) - 1)
    kk = np.clip(k, 0, len(z) - 2)
    part = np.clip(np.minimum(z[kk + 1], b) - np.clip(y, a, b), 0, None)
    part = np.where(inside & (y > z[kk]), part, 0.0)
    Fk = F[np.arange(len(y)), kk]
    s += (1 - 2 * Fk) * part
    # below support: F=0 on (-inf, z0): loss 1 for z >= y
    s += np.where(y < z[0], np.clip(np.minimum(z[0], b) - np.maximum(y, a), 0, None), 0.0)
    # above support: F=1 on [zm, inf): loss 1 for z < y
    s += np.where(y > z[-1], np.clip(np.minimum(y, b) - np.maximum(z[-1], a), 0, None), 0.0)
    return s

# ---------------------------------------------------------------- Gaussian EMOS
def gauss_crps(mu, sig, y):
    zz = (y - mu) / sig
    return sig * (zz * (2 * norm.cdf(zz) - 1) + 2 * norm.pdf(zz) - 1 / np.sqrt(np.pi))

def _g_std(u):
    """Integral of Phi(s)^2 over (-inf, u] for the standard normal CDF Phi (closed form)."""
    return u * norm.cdf(u) ** 2 + 2 * norm.pdf(u) * norm.cdf(u) - norm.cdf(np.sqrt(2) * u) / np.sqrt(np.pi)

def _gauss_tw_upper(mu, sig, y, t):
    """Closed-form CRPS of N(mu, sig^2) with weight 1{z > t}: integral over (t, inf) of (F(z) - 1{y <= z})^2."""
    u, uy = (t - mu) / sig, (y - mu) / sig
    # y <= t: int_t^inf (1-F)^2 = sig * G(-u);  y > t: int_t^y F^2 + int_y^inf (1-F)^2
    above = sig * (_g_std(uy) - _g_std(u) + _g_std(-uy))
    below = sig * _g_std(-u)
    s = np.where(y > t, above, below)
    s = np.where(np.isneginf(t), gauss_crps(mu, sig, y), s)
    return np.clip(s, 0.0, None)

def gauss_twcrps(mu, sig, y, a=-np.inf, b=np.inf):
    """Indicator-weighted CRPS of N(mu, sig^2) over the weight region (a, b), in closed form.

    Upper tail (a, inf) follows from splitting the threshold integral at max(y, a); the lower tail
    (-inf, b) is the upper-tail score of the reflected forecast and observation. A finite interval is
    the difference of two upper-tail scores. Scores are non-negative for any threshold position."""
    mu, sig, y = np.broadcast_arrays(np.asarray(mu, float), np.asarray(sig, float), np.asarray(y, float))
    a = np.broadcast_to(np.asarray(a, float), mu.shape); b = np.broadcast_to(np.asarray(b, float), mu.shape)
    if np.all(np.isinf(b)) and np.all(np.isinf(a)):
        return gauss_crps(mu, sig, y)
    if np.all(np.isinf(b)):
        return _gauss_tw_upper(mu, sig, y, a)
    if np.all(np.isinf(a)):
        return _gauss_tw_upper(-mu, sig, -y, -b)
    return np.clip(_gauss_tw_upper(mu, sig, y, a) - _gauss_tw_upper(mu, sig, y, b), 0.0, None)

def emos_fit(x, y, doy):
    S, C = np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)
    A = np.column_stack([np.ones_like(x), x, S, C])
    beta = np.linalg.lstsq(A, y, rcond=None)[0]
    r = y - A @ beta
    p0 = np.r_[beta, np.log(r.std() + 1e-3), 0, 0]
    B = np.column_stack([np.ones_like(x), S, C])
    def f(p):                                   # mean CRPS and analytic gradient
        mu = A @ p[:4]; sig = np.exp(B @ p[4:])
        zz = (y - mu) / sig; Phi = norm.cdf(zz); phi = norm.pdf(zz)
        crps = sig * (zz * (2 * Phi - 1) + 2 * phi - 1 / np.sqrt(np.pi))
        dmu = -(2 * Phi - 1); dsig = 2 * phi - 1 / np.sqrt(np.pi)
        g = np.r_[A.T @ dmu, B.T @ (dsig * sig)] / len(y)
        return crps.mean(), g
    res = minimize(f, p0, jac=True, method="L-BFGS-B")
    return res.x

def emos_predict(p, x, doy):
    S, C = np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)
    return p[0] + p[1] * x + p[2] * S + p[3] * C, np.exp(p[4] + p[5] * S + p[6] * C)

# ---------------------------------------------------------------- DM test
def dm_test(d, lag):
    d = np.asarray(d, float)
    ok = np.isfinite(d); n = int(ok.sum())
    if n < 50:
        return np.nan, np.nan
    m = d[ok].mean(); e = np.where(ok, d - m, 0.0)
    v = e @ e / n
    for k in range(1, min(lag, len(d) - 1) + 1):
        v += 2 * (1 - k / (lag + 1)) * (e[k:] @ e[:-k]) / n
    if v <= 0:
        return np.nan, np.nan
    stat = m / np.sqrt(v / n)
    return stat, 2 * norm.sf(abs(stat))

# ---------------------------------------------------------------- worker
G = {}
def expected_models(var, setname):
    models = SETS[setname]["models"]
    return [m for m in models if m == "gfs" or m.startswith("graphcast")] if var == "pr24" else models

def check_inputs(var, setname):
    d = AN / var
    required = [d / f"{m}.npy" for m in expected_models(var, setname)]
    required += [d / f"{k}.npy" for k in ("obs", "leads", "inits", "stations")]
    missing = [str(f) for f in required if not f.exists()]
    if missing:
        raise FileNotFoundError("Missing comparison inputs; rebuild analysis cubes:\n" + "\n".join(missing))

def init_worker(var, setname):
    cfg = SETS[setname]
    d = AN / var
    G["var"] = var; G["cfg"] = cfg
    G["obs"] = np.load(d / "obs.npy", mmap_mode="r")
    G["fc"] = {m: np.load(d / f"{m}.npy", mmap_mode="r") for m in expected_models(var, setname)}
    G["leads"] = np.load(d / "leads.npy"); G["inits"] = pd.to_datetime(np.load(d / "inits.npy"))
    G["clim"] = pd.read_csv(AN / "climatology_2011_2020.csv", dtype={"station_id": str}).set_index("station_id")
    G["sids"] = np.load(d / "stations.npy")
    sa = pd.read_csv(ROOT / "data" / "interim" / "stations_analysis.csv", dtype={"station_id": str}).set_index("station_id")
    G["use"] = sa[f"use_{var}"].reindex(G["sids"]).values if f"use_{var}" in sa else np.ones(len(G["sids"]), bool)

def lead_hours(var, lead):
    return int(lead) if var in ("t2", "ws") else 24 * int(lead)

def save_lead(var, lead):
    return True

def observation_end(var, lead, inits):
    """Time when the target observation is available (no reporting delay assumed)."""
    hours = int(lead) if var in ("t2", "ws") else 24 * int(lead)
    if var in ("tmax", "tmin"):
        hours += 12       # 12 UTC on the labelled date, not its midnight
    return inits + pd.to_timedelta(hours, unit="h")

def oos_masks(inits, period, train_end, var, lead):
    """Freeze calibration before the first test initialization; purge future labels."""
    cutoff = pd.Timestamp(train_end) + pd.Timedelta(days=1)
    train = period & (inits < cutoff) & (observation_end(var, lead, inits) < cutoff)
    test = period & (inits >= cutoff)
    return train, test

def thresholds(var, sid, lead):
    out = []
    c = G["clim"].loc[sid]
    for name, col, tail in THRESH[var]:
        if isinstance(col, str):
            h = lead_hours(var, lead) % 24
            t = c[col.format(h=h)]
        else:
            t = col
        out.append((name, t, tail))
    return out

def work(s):
    var, cfg = G["var"], G["cfg"]
    models = list(G["fc"]); ref = cfg["ref"]
    sid = G["sids"][s]; inits = G["inits"]
    period = (inits >= pd.Timestamp(cfg["start"])) & (inits <= pd.Timestamp(cfg.get("end", "2100-01-01")))
    doy = inits.dayofyear.values.astype(float)
    rows, dms, perday = [], [], {}
    if not G["use"][s]:
        return s, rows, dms, perday
    for li, lead in enumerate(G["leads"]):
        y_all = np.asarray(G["obs"][:, li, s], dtype=float)
        X = np.stack([np.asarray(G["fc"][m][:, li, s], dtype=float) for m in models])
        mask = period & np.isfinite(y_all) & np.isfinite(X).all(axis=0)
        n = int(mask.sum())
        if n < 150:
            continue
        thr = thresholds(var, sid, lead)
        if any(not np.isfinite(t) for _, t, _ in thr):
            continue
        bounds = [((t, np.inf) if tail == "upper" else (-np.inf, t)) for _, t, tail in thr]
        y = y_all[mask]
        keep = save_lead(var, lead)
        losses = {}
        # climatological reference (in-sample)
        zc = np.unique(y); Fc = np.searchsorted(np.sort(y), zc, side="right") / len(y)
        Fc = np.broadcast_to(Fc, (len(y), len(zc)))
        rec0 = {"station_id": sid, "lead": int(lead), "n": n}
        rows.append({**rec0, "model": "clim", "pcrps": crps_step(zc, Fc, y).mean(),
                     **{f"tw{k+1}": crps_step(zc, Fc, y, *bounds[k]).mean() for k in range(len(bounds))}})
        # out-of-sample masks
        if cfg["train_end"] is not None:
            train, test = oos_masks(inits, period, cfg["train_end"], var, lead)
            mtr, mte = mask & train, mask & test
            oos_ok = mtr.sum() >= 150 and mte.sum() >= 100
        else:
            oos_ok = False
        if oos_ok:
            ytr, yte = y_all[mtr], y_all[mte]
            zt = np.unique(ytr); Ft = np.searchsorted(np.sort(ytr), zt, side="right") / len(ytr)
            Ft = np.broadcast_to(Ft, (len(yte), len(zt)))
            rows[-1].update({"oos_idr": crps_step(zt, Ft, yte).mean(), "oos_tw1": crps_step(zt, Ft, yte, *bounds[0]).mean(),
                             "oos_tw2": crps_step(zt, Ft, yte, *bounds[1]).mean(), "n_test": int(mte.sum()), "n_train": int(mtr.sum())})
        for k, m in enumerate(models):
            x = X[k][mask]
            e = x - y
            ux, z, F, inv = idr_fit(x, y)
            Fo = F[inv]
            c = crps_step(z, Fo, y); tw = [crps_step(z, Fo, y, *b) for b in bounds]
            r = {**rec0, "model": m, "bias": e.mean(), "mae": np.abs(e).mean(), "rmse": np.sqrt((e ** 2).mean()),
                 "pcrps": c.mean(), "tw1": tw[0].mean(), "tw2": tw[1].mean()}
            L = {"pcrps": (mask, c), "tw1": (mask, tw[0]), "tw2": (mask, tw[1]), "abs_err": (mask, np.abs(e))}
            if oos_ok:
                xtr, xte = X[k][mtr], X[k][mte]
                ux2, z2, F2, _ = idr_fit(xtr, ytr)
                P = idr_predict(ux2, F2, xte)
                co = crps_step(z2, P, yte); cot = crps_step(z2, P, yte, *bounds[0]); cot2 = crps_step(z2, P, yte, *bounds[1])
                r.update({"oos_idr": co.mean(), "oos_tw1": cot.mean(), "oos_tw2": cot2.mean(),
                          "oos_mae": np.abs(xte - yte).mean(), "n_test": int(mte.sum()), "n_train": int(mtr.sum())})
                L["oos_idr"] = (mte, co); L["oos_tw1"] = (mte, cot); L["oos_tw2"] = (mte, cot2)
                if var in EMOS_VARS and (var != "t2" or lead % 24 == 0):
                    p = emos_fit(xtr, ytr, doy[mtr])
                    mu, sg = emos_predict(p, xte, doy[mte])
                    ce = gauss_crps(mu, sg, yte)
                    ct = gauss_twcrps(mu, sg, yte, *bounds[0])
                    r.update({"oos_emos": ce.mean(), "oos_emos_tw1": ct.mean()})
                    L["oos_emos"] = (mte, ce); L["oos_emos_tw1"] = (mte, ct)
            rows.append(r)
            losses[m] = L
        if keep:
            hlag = max(int(np.ceil(lead_hours(var, lead) / 24)), int(np.ceil(n ** (1 / 3))))
            for m in models:
                for sc, (mk, v) in losses[m].items():
                    full = np.full(len(inits), np.nan, "f4"); full[mk] = v
                    perday[(int(lead), m, sc)] = full
            refs = [ref] + ([m.replace("_ifs", "") for m in models if m.endswith("_ifs")] if any(m.endswith("_ifs") for m in models) else [])
            for m in models:
                for rm in set(refs):
                    if m == rm or rm not in losses or (rm != ref and m != rm + "_ifs"):
                        continue
                    for sc in losses[m]:
                        if sc not in losses[rm]:
                            continue
                        mk = losses[m][sc][0]
                        if not np.array_equal(mk, losses[rm][sc][0]):
                            raise ValueError("DM comparison requires paired cases")
                        d = losses[m][sc][1] - losses[rm][sc][1]
                        calendar_d = np.full(len(inits), np.nan); calendar_d[mk] = d
                        st_, p_ = dm_test(calendar_d, hlag)
                        dms.append({"station_id": sid, "lead": int(lead), "model": m, "ref": rm, "score": sc,
                                    "mean_diff": float(d.mean()), "ref_mean": float(losses[rm][sc][1].mean()), "dm": st_, "p": p_})
    return s, rows, dms, perday

def main(var, setname, workers):
    require_current(AN / "provenance.json")
    check_inputs(var, setname)   # fail in the parent, not repeatedly in Pool initializers
    t0 = time.time()
    sids = np.load(AN / var / "stations.npy")
    st = pd.read_csv(ROOT / "data" / "interim" / "stations_analysis.csv", dtype={"station_id": str}).set_index("station_id")
    zones = st.loc[sids, "zone"].values
    inits = pd.to_datetime(np.load(AN / var / "inits.npy"))
    rows, dms = [], []
    sums, cnts = {}, {}
    with Pool(workers, initializer=init_worker, initargs=(var, setname)) as pool:
        for k, (s, r, d, pdy) in enumerate(pool.imap_unordered(work, range(len(sids)))):
            rows += r; dms += d
            for key, v in pdy.items():
                for zz in ("All", zones[s]):
                    kk = (zz,) + key
                    if kk not in sums:
                        sums[kk] = np.zeros(len(inits)); cnts[kk] = np.zeros(len(inits))
                    f = np.isfinite(v); sums[kk][f] += v[f]; cnts[kk][f] += 1
            if k % 50 == 0:
                print(f"{var}/{setname}: {k}/{len(sids)} stations {time.time()-t0:.0f}s rss {peak_rss_mb()}MB", flush=True)
    df = pd.DataFrame(rows); df.insert(1, "zone", st.loc[df.station_id, "zone"].values)
    df.to_parquet(OUT / f"{setname}_{var}_station.parquet", index=False)
    pd.DataFrame(dms).to_parquet(OUT / f"{setname}_{var}_dm.parquet", index=False)
    keys = sorted(sums)
    np.savez_compressed(OUT / f"{setname}_{var}_daily.npz", keys=np.array(["|".join(map(str, k)) for k in keys]),
                        sums=np.array([sums[k] for k in keys], "f8"), cnts=np.array([cnts[k] for k in keys], "f4"),
                        inits=inits.values.astype("datetime64[s]"))
    log_json(f"07_{setname}_{var}.json", {"stations": len(sids), "rows": len(df), "wall_s": round(time.time() - t0), "peak_rss_MB": peak_rss_mb()})
    write_metadata(OUT / f"{setname}_{var}_meta.json", comparison_set=setname, variable=var,
                   leads=np.load(AN / var / "leads.npy").tolist(),
                   train_initialization_end=SETS[setname]["train_end"], future_labels_purged=True)
    print(f"done {var}/{setname} in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--var", required=True); ap.add_argument("--set", default="core")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); main(a.var, a.set, a.workers)
