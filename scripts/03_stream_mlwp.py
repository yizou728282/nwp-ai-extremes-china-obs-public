"""Stream station-point forecasts from the AIWP archive (byte-range GETs of single 2-D chunks; whole
NetCDF files are never downloaded) for the station superset, bilinear interpolation.
Output per (model, year, month): data/fc/<model>/<YYYY-MM>.npz with
   inits [n] (datetime64[s]),  <var> [n, 41, n_station] float32   (t2 degC, u10/v10 m/s, apcp mm/6h)
Resumable: months already written are skipped.  A month is written only if every chunk succeeded or
was confirmed missing (404/zero-size -> NaN).
Usage: python 03_stream_mlwp.py MODEL [MODEL...] --years 2023 [--threads 16]
"""
import os, time, zlib, argparse, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from numcodecs import Shuffle
from common import ALL_AI, AI_VARS, FC, INTERIM, NLEAD, get_range, log_json, peak_rss_mb
import importlib; refs_mod = importlib.import_module("03_refs")

unshuffle = Shuffle(elementsize=4)

def weights(st):
    lat = st.lat.values; lon = np.mod(st.lon.values, 360.0)
    fi = (90 - lat) / 0.25; fj = lon / 0.25
    i0 = np.floor(fi).astype(int); j0 = np.floor(fj).astype(int)
    return i0, np.minimum(i0 + 1, 720), j0, (j0 + 1) % 1440, (fi - i0).astype("f4"), (fj - j0).astype("f4")

def interp(a, w):
    i0, i1, j0, j1, wi, wj = w
    return a[i0, j0] * (1 - wi) * (1 - wj) + a[i1, j0] * wi * (1 - wj) + a[i0, j1] * (1 - wi) * wj + a[i1, j1] * wi * wj

def run(model, year, threads):
    arch = ALL_AI[model]
    st = pd.read_csv(INTERIM / "stations_superset.csv", dtype={"station_id": str})
    w = weights(st); ns = len(st)
    refs = refs_mod.build(arch, year)
    if refs.empty:
        print(model, year, "no data"); return
    out_dir = FC / model; out_dir.mkdir(parents=True, exist_ok=True)
    refs["month"] = refs.init_time.dt.strftime("%Y-%m")
    for mon, g in refs.groupby("month"):
        f = out_dir / f"{mon}.npz"
        if f.exists():
            continue
        t0 = time.time()
        inits = np.sort(g.init_time.unique()); ipos = {t: k for k, t in enumerate(inits)}
        vars_ = AI_VARS[arch]
        arr = {v: np.full((len(inits), NLEAD, ns), np.nan, "f4") for v in vars_}
        tasks = list(g[["init_time", "var", "lead_idx", "url", "offset", "size"]].itertuples(index=False))
        def fetch(t):
            b = get_range(t.url, t.offset, t.size)
            if b is None:
                return t, None
            a = np.frombuffer(unshuffle.decode(zlib.decompress(b)), "<f4").reshape(721, 1440)
            return t, interp(a, w)
        nbytes = 0; nmiss = 0
        with ThreadPoolExecutor(threads) as ex:
            for t, v in ex.map(fetch, tasks):
                if v is None:
                    nmiss += 1; continue
                nbytes += t.size
                arr[t.var][ipos[t.init_time], t.lead_idx] = v
        if "t2" in arr:
            arr["t2"] -= 273.15
        if "apcp" in arr:
            arr["apcp"] *= 1000.0            # m per 6 h -> mm per 6 h (verified: global mean ~0.6 mm/6h)
        np.savez(f, inits=inits.astype("datetime64[s]"), station_id=np.array(st.station_id.tolist(), dtype="U16"), **arr)
        dt = time.time() - t0
        print(f"{model} {mon}: {len(inits)} inits, {len(tasks)} chunks, miss {nmiss}, {nbytes/1e6:.0f} MB, {dt:.0f}s, "
              f"{nbytes/1e6/dt:.0f} MB/s, rss {peak_rss_mb()} MB", flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("models", nargs="+"); ap.add_argument("--years", default="2023")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("AIX_THREADS", "16")))
    a = ap.parse_args()
    y0, y1 = map(int, a.years.split("-")) if "-" in a.years else (int(a.years),) * 2
    for y in range(y0, y1 + 1):
        for m in a.models:
            run(m, y, a.threads)
