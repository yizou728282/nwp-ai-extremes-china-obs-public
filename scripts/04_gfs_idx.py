"""Station-point GFS (0.25 deg, s3://noaa-gfs-bdp-pds) via .idx byte ranges: per (init, lead) the idx is
fetched once, then the 4 GRIB2 messages TMP 2 m, UGRD/VGRD 10 m and the 6-h APCP bucket.
GFS APCP buckets alternate (0-6, 0-12 | 12-18, 12-24 ...): the 6-h amount is derived as
  apcp6(f) = acc(f) - acc(f-6) when the message at f accumulates over more than 6 h.
Output per month: data/fc/gfs/<YYYY-MM>.npz  inits, t2 (degC), u10, v10, apcp (mm/6h) [n, 41, n_station]
Usage: python 04_gfs_idx.py --years 2023 [--threads 16]
"""
import os, re, time, argparse, numpy as np, pandas as pd, eccodes
from concurrent.futures import ThreadPoolExecutor
from common import GFS, FC, INTERIM, NLEAD, STUDY_START, STUDY_END, session, get_range, peak_rss_mb
w_mod = __import__("importlib").import_module("03_stream_mlwp")

KEYS = {"t2": ":TMP:2 m above ground:", "u10": ":UGRD:10 m above ground:", "v10": ":VGRD:10 m above ground:"}

def url(t, h):
    sub = "atmos/" if t >= pd.Timestamp("2021-03-23") else ""
    return f"{GFS}/gfs.{t:%Y%m%d}/00/{sub}gfs.t00z.pgrb2.0p25.f{h:03d}"

def decode(b, w):
    if b is None or len(b) < 8 or b[-4:] != b"7777":
        raise ValueError("truncated GRIB message")
    gid = eccodes.codes_new_from_message(b)
    try:
        ni, nj = eccodes.codes_get(gid, "Ni"), eccodes.codes_get(gid, "Nj")
        la1 = eccodes.codes_get(gid, "latitudeOfFirstGridPointInDegrees")
        a = eccodes.codes_get_values(gid).reshape(nj, ni)
    finally:
        eccodes.codes_release(gid)
    assert (ni, nj) == (1440, 721) and abs(la1 - 90) < 1e-6
    return w_mod.interp(a, w).astype("f4")

def fetch(task, w):
    t, h = task
    u = url(t, h)
    for k in range(4):
        try:
            r = session().get(u + ".idx", timeout=60); break
        except Exception:
            time.sleep(2 ** k)
    else:
        return t, h, None, "idx error", 0
    if r.status_code != 200:
        return t, h, None, "idx missing", 0
    lines = r.text.strip().split("\n")
    def find(pred):
        k = next((n for n, l in enumerate(lines) if pred(l)), None)
        if k is None:
            return None
        s = int(lines[k].split(":")[1]); e = int(lines[k + 1].split(":")[1]) if k + 1 < len(lines) else s + 3_000_000
        return s, e - s, lines[k]
    out, nb = {}, 0
    tag = ":anl:" if h == 0 else f":{h} hour fcst:"
    for v, key in KEYS.items():
        m = find(lambda l: key in l and tag in l)
        if m:
            for attempt in range(3):
                try:
                    b = get_range(u, m[0], m[1]); nb += len(b); out[v] = decode(b, w); break
                except Exception:
                    time.sleep(2)
    if h > 0:
        m = find(lambda l: ":APCP:surface:" in l and re.search(r":\d+-%d hour acc fcst:" % h, l))
        if m:
            a0 = int(re.search(r":(\d+)-%d hour acc" % h, m[2]).group(1))
            for attempt in range(3):
                try:
                    b = get_range(u, m[0], m[1]); nb += len(b); out["apcp_acc"] = (a0, decode(b, w)); break
                except Exception:
                    time.sleep(2)
    return t, h, out, None, nb

def run(year, threads):
    st = pd.read_csv(INTERIM / "stations_superset.csv", dtype={"station_id": str})
    w = w_mod.weights(st); ns = len(st)
    d = FC / "gfs"; d.mkdir(parents=True, exist_ok=True)
    inits_all = pd.date_range(max(f"{year}-01-01", STUDY_START), min(f"{year}-12-31", STUDY_END), freq="D")
    for mon, inits in pd.Series(inits_all, index=inits_all).groupby(inits_all.strftime("%Y-%m")):
        f = d / f"{mon}.npz"
        if f.exists():
            continue
        t0 = time.time(); inits = list(inits); ipos = {t: k for k, t in enumerate(inits)}
        arr = {v: np.full((len(inits), NLEAD, ns), np.nan, "f4") for v in ("t2", "u10", "v10", "apcp")}
        acc = {}
        tasks = [(t, h) for t in inits for h in range(0, 241, 6)]
        nb = 0; miss = []
        with ThreadPoolExecutor(threads) as ex:
            for t, h, out, err, b in ex.map(lambda x: fetch(x, w), tasks):
                nb += b
                if out is None:
                    miss.append(f"{t.date()}+{h} {err}"); continue
                for v in ("t2", "u10", "v10"):
                    if v in out:
                        arr[v][ipos[t], h // 6] = out[v]
                if "apcp_acc" in out:
                    acc[(t, h)] = out["apcp_acc"]
        for (t, h), (a0, x) in acc.items():          # kg m-2 == mm
            if h - a0 == 6:
                arr["apcp"][ipos[t], h // 6] = x
            elif (t, h - 6) in acc and acc[(t, h - 6)][0] == a0:
                arr["apcp"][ipos[t], h // 6] = np.maximum(x - acc[(t, h - 6)][1], 0)
        arr["t2"] -= 273.15
        np.savez(f, inits=np.array(inits, dtype="datetime64[s]"), station_id=np.array(st.station_id.tolist(), dtype="U16"), **arr)
        dt = time.time() - t0
        print(f"gfs {mon}: {len(inits)} inits, {len(tasks)} fields, missing {len(miss)} {miss[:3]}, {nb/1e6:.0f} MB, {dt:.0f}s, rss {peak_rss_mb()} MB", flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--years", default="2023")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("AIX_THREADS", "16")))
    a = ap.parse_args()
    y0, y1 = map(int, a.years.split("-")) if "-" in a.years else (int(a.years),) * 2
    for y in range(y0, y1 + 1):
        run(y, a.threads)
