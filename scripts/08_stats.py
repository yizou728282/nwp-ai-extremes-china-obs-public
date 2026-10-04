"""Aggregate scores, date-block bootstrap CIs and Diebold-Mariano tests with Benjamini-Hochberg FDR.
Inputs : results/scores/{set}_{var}_{station.parquet, dm.parquet, daily.npz}
Outputs: results/tables/{set}_{var}_summary.csv   mean scores per lead x model (+ skill vs GFS and vs climatology)
         results/tables/{set}_{var}_bootstrap.csv skill vs reference with 95% moving-block bootstrap CI, per zone
         results/tables/{set}_{var}_dm_fdr.csv    fraction of stations significantly better / worse (BH-FDR 5%)
         results/scores/{set}_{var}_dm_fdr.parquet station-level q-values (for maps)
Usage: python 08_stats.py --var t2 --set core
"""
import argparse, numpy as np, pandas as pd
from common import RESULTS
from provenance import require_current, write_metadata

SC = RESULTS / "scores"; TB = RESULTS / "tables"; TB.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261002)

def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = np.isfinite(p)
    pv = p[ok]; n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv); r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n); qq[o] = np.minimum(r, 1); q[ok] = qq
    return q

def block_boot(a, b, block=14, B=1000, available=None):
    """Paired, calendar-day block bootstrap of pooled skill.

    Zero losses are valid. Missing days remain on the calendar axis with zero
    contribution, rather than being joined to neighbouring observed days.
    ``available`` distinguishes an observed zero loss from an empty daily sum.
    Only the unused leading/trailing part of the calendar axis is trimmed.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.ndim != 1 or a.shape != b.shape or block < 1 or B < 1:
        raise ValueError("paired one-dimensional losses and positive block/B required")
    ok = np.isfinite(a) & np.isfinite(b)
    if available is not None:
        available = np.asarray(available, bool)
        if available.shape != a.shape:
            raise ValueError("availability must match the calendar axis")
        ok &= available
    if ok.sum() < 2 * block:
        return np.nan, np.nan, np.nan
    first, last = np.flatnonzero(ok)[[0, -1]]
    a, b = np.where(ok, a, 0.0)[first:last+1], np.where(ok, b, 0.0)[first:last+1]
    n = len(a)
    if b.sum() <= 0:
        return np.nan, np.nan, np.nan
    nb = int(np.ceil(n / block)); starts = rng.integers(0, n - block + 1, size=(B, nb))
    idx = (starts[:, :, None] + np.arange(block)[None, None, :]).reshape(B, -1)[:, :n]
    den = b[idx].sum(1)
    valid = den > 0
    if valid.sum() < 0.95 * B:
        return 1 - a.sum() / b.sum(), np.nan, np.nan
    sk = 1 - a[idx[valid]].sum(1) / den[valid]
    return 1 - a.sum() / b.sum(), np.quantile(sk, 0.025), np.quantile(sk, 0.975)

def main(var, setname, ref="gfs"):
    require_current(SC / f"{setname}_{var}_meta.json")
    st = pd.read_parquet(SC / f"{setname}_{var}_station.parquet")
    scores = [c for c in ["pcrps", "tw1", "tw2", "oos_idr", "oos_tw1", "oos_tw2", "oos_emos", "oos_emos_tw1", "mae", "rmse", "bias", "oos_mae"] if c in st]
    # use only (station, lead) cells where every model has a value -> same stations for all
    agg = st.groupby(["lead", "model"])[scores].mean().reset_index()
    nst = st.groupby(["lead", "model"]).station_id.nunique().rename("n_stations").reset_index()
    agg = agg.merge(nst, on=["lead", "model"])
    for s in scores:
        if s in ("bias",):
            continue
        r = agg[agg.model == ref].set_index("lead")[s]; c = agg[agg.model == "clim"].set_index("lead")[s] if s in agg[agg.model == "clim"].dropna(axis=1, how="all") else None
        agg[f"{s}_ss_{ref}"] = 1 - agg[s] / agg.lead.map(r)
        if c is not None:
            agg[f"{s}_ss_clim"] = 1 - agg[s] / agg.lead.map(c)
    agg.round(5).to_csv(TB / f"{setname}_{var}_summary.csv", index=False)
    # bootstrap from daily sums
    z = np.load(SC / f"{setname}_{var}_daily.npz")
    keys = [k.split("|") for k in z["keys"]]; sums = z["sums"]
    D = {tuple(k): i for i, k in enumerate(keys)}
    rows = []
    for (zone, lead, model, score), i in D.items():
        if model == ref:
            continue
        refs = [ref] + ([model.replace("_ifs", "")] if model.endswith("_ifs") else [])
        for rm in refs:
            j = D.get((zone, lead, rm, score))
            if j is None:
                continue
            if not np.array_equal(z["cnts"][i], z["cnts"][j]):
                raise ValueError(f"unpaired daily samples: {zone}/{lead}/{model}/{rm}/{score}")
            available = z["cnts"][i] > 0
            s, lo, hi = block_boot(sums[i], sums[j], available=available)
            rows.append({"zone": zone, "lead": int(lead), "model": model, "ref": rm, "score": score, "skill": s, "ci_lo": lo, "ci_hi": hi,
                         "n_days": int((z["cnts"][i] > 0).sum())})
    bt = pd.DataFrame(rows).sort_values(["score", "zone", "model", "ref", "lead"])
    bt.round(5).to_csv(TB / f"{setname}_{var}_bootstrap.csv", index=False)
    # DM + FDR
    dm = pd.read_parquet(SC / f"{setname}_{var}_dm.parquet")
    dm["q"] = dm.groupby(["lead", "model", "ref", "score"])["p"].transform(bh)
    dm.to_parquet(SC / f"{setname}_{var}_dm_fdr.parquet", index=False)
    g = dm.groupby(["lead", "model", "ref", "score"])
    fr = pd.DataFrame({"n": g.size(),
                       "frac_better_fdr05": g.apply(lambda x: ((x.q < 0.05) & (x.mean_diff < 0)).mean()),
                       "frac_worse_fdr05": g.apply(lambda x: ((x.q < 0.05) & (x.mean_diff > 0)).mean()),
                       "frac_better_p05_unadj": g.apply(lambda x: ((x.p < 0.05) & (x.mean_diff < 0)).mean()),
                       "frac_worse_p05_unadj": g.apply(lambda x: ((x.p < 0.05) & (x.mean_diff > 0)).mean()),
                       "frac_better_raw": g.apply(lambda x: (x.mean_diff < 0).mean())}).reset_index()
    fr.round(4).to_csv(TB / f"{setname}_{var}_dm_fdr.csv", index=False)
    write_metadata(TB / f"{setname}_{var}_meta.json", comparison_set=setname, variable=var,
                   bootstrap_block_days=14, bootstrap_replicates=1000, valid_zero_losses_retained=True,
                   missing_calendar_days_retained=True, aggregate_weighting="pooled station-case sums",
                   summary_weighting="equal station means")
    print(f"{setname}/{var}: summary {len(agg)} rows, bootstrap {len(bt)}, dm {len(dm)}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--var", required=True); ap.add_argument("--set", default="core")
    a = ap.parse_args(); main(a.var, a.set)
