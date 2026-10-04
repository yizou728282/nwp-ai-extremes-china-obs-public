"""Chunk tables (url, offset, size per init x lead x variable) for the AIWP archive, one parquet per
(archive, year) in data/refs/. Two sources:
  (a) the kerchunk parquet index (parquet/<ARCH>_combined_all.parq), which ends 2025-04-01;
  (b) "03b": for 00 UTC files on S3 that are missing from the index (all files after 2025-04-01 and
      any gaps), the chunk byte offsets are read directly from the NetCDF4/HDF5 file over HTTP with
      h5py (Dataset.id.get_chunk_info_by_coord), fetching only metadata blocks (~0.5-1 s per file).
Both give the same encoding (zlib + 4-byte shuffle) for the 2-D (lat, lon) chunks of the
(41, 721, 1440) variables.  Also stores lat/lon grid checks.
Usage: python 03_refs.py ARCH [ARCH...] --years 2021-2025
"""
import io, re, json, sys, time, argparse, numpy as np, pandas as pd, fsspec, h5py
from concurrent.futures import ThreadPoolExecutor
from common import MLWP, REFS, AI_VARS, NLEAD, STUDY_START, STUDY_END, session, get_range, log_json

S3 = "s3://noaa-oar-mlwp-data"
_index_cache = {}

def s3_keys(prefix):
    import urllib.parse, xml.etree.ElementTree as ET
    keys, token = [], None
    while True:
        q = f"{MLWP}/?list-type=2&prefix={prefix}" + (f"&continuation-token={urllib.parse.quote(token, safe='')}" if token else "")
        x = ET.fromstring(session().get(q, timeout=60).content)
        ns = {"s": x.tag.split("}")[0].strip("{")}
        keys += [c.find("s:Key", ns).text for c in x.findall("s:Contents", ns)]
        t = x.find("s:NextContinuationToken", ns)
        if t is None:
            return keys
        token = t.text

def index(arch):
    """Kerchunk index: init times and lazily loaded ref record files per variable."""
    if arch in _index_cache:
        return _index_cache[arch]
    base = f"{MLWP}/parquet/{arch}_combined_all.parq"
    meta = json.loads(session().get(f"{base}/.zmetadata", timeout=60).text)
    rec = meta["record_size"]
    d = pd.read_parquet(io.BytesIO(session().get(f"{base}/init_time/refs.0.parq", timeout=60).content)).iloc[0]
    b = d["raw"] if d["raw"] is not None else get_range(d["path"].replace(S3, MLWP), int(d["offset"]), int(d["size"]))
    inits = pd.to_datetime(np.frombuffer(b, "<M8[ns]"))
    ix = {"base": base, "rec": rec, "inits": inits, "refs": {}}
    _index_cache[arch] = ix
    return ix

def index_refs(ix, var, r):
    k = (var, r)
    if k not in ix["refs"]:
        ix["refs"][k] = pd.read_parquet(io.BytesIO(session().get(f"{ix['base']}/{var}/refs.{r}.parq", timeout=120).content))
    return ix["refs"][k]

def h5_chunks(url, vars_):
    """03b: chunk offsets straight from the HDF5 file."""
    fs = fsspec.filesystem("http", block_size=2**20, cache_type="bytes")
    rows = []
    with fs.open(url) as fo, h5py.File(fo, "r") as f:
        lat = f["latitude"][:]; lon = f["longitude"][:]
        assert abs(lat[0] - 90) < 1e-4 and abs(lat[-1] + 90) < 1e-4 and abs(lon[0]) < 1e-4 and len(lon) == 1440, "grid changed"
        for v in vars_:
            if v not in f:
                continue
            ds = f[v]
            sh = ds.shape[-3:]
            assert sh == (NLEAD, 721, 1440) and ds.chunks[-2:] == (721, 1440), (v, ds.shape, ds.chunks)
            filt = ds.id.get_create_plist()
            nf = filt.get_nfilters(); codes = [filt.get_filter(i)[0] for i in range(nf)]
            assert set(codes) <= {1, 2}, codes          # 1=deflate, 2=shuffle
            for li in range(NLEAD):
                coord = (0,) * (ds.ndim - 3) + (li, 0, 0)
                ci = ds.id.get_chunk_info_by_coord(coord)
                rows.append((v, li, int(ci.byte_offset), int(ci.size)))
    return rows

def build(arch, year):
    out = REFS / f"{arch}_{year}.parquet"
    if out.exists():
        return pd.read_parquet(out)
    t0 = time.time()
    vars_ = AI_VARS[arch]
    want = pd.date_range(max(f"{year}-01-01", STUDY_START), min(f"{year}-12-31", STUDY_END), freq="D")
    ix = index(arch)
    pos = pd.Series(np.arange(len(ix["inits"])), index=ix["inits"]); pos = pos[~pos.index.duplicated()]
    rows, from_index = [], set()
    for t in want:
        if t not in pos.index:
            continue
        ok = True; tmp = []
        for v in vars_:
            for li in range(NLEAD):
                flat = int(pos[t]) * NLEAD + li
                r = index_refs(ix, v, flat // ix["rec"]).iloc[flat % ix["rec"]]
                if not isinstance(r["path"], str) or not r["size"]:
                    ok = False; break
                tmp.append((t, v, li, r["path"].replace(S3, MLWP), int(r["offset"]), int(r["size"]), "kerchunk"))
            if not ok:
                break
        if ok:
            rows += tmp; from_index.add(t)
    # files on S3 not covered by the index -> 03b
    keys = [k for k in s3_keys(f"{arch}/{year}/") if re.search(r"_\d{8}00_f000_f240_06\.nc$", k)]
    s3 = {pd.Timestamp(re.search(r"_(\d{8})00_f000", k).group(1)): k for k in keys}
    todo = [(t, s3[t]) for t in want if t not in from_index and t in s3]
    def one(tk):
        t, k = tk
        url = f"{MLWP}/{k}"
        for a in range(3):
            try:
                return [(t, v, li, url, off, size, "h5py") for v, li, off, size in h5_chunks(url, vars_)]
            except AssertionError:
                raise
            except Exception as e:
                err = e; time.sleep(3 * (a + 1))
        print("03b failed", url, err, flush=True); return []
    with ThreadPoolExecutor(8) as ex:
        for r in ex.map(one, todo):
            rows += r
    df = pd.DataFrame(rows, columns=["init_time", "var", "lead_idx", "url", "offset", "size", "source"])
    df.to_parquet(out, index=False)
    n_init = df.init_time.nunique() if len(df) else 0
    s = {"arch": arch, "year": year, "wanted": len(want), "inits": n_init, "from_index": len(from_index), "from_h5py": len(todo),
         "missing": [str(t.date()) for t in want if t not in set(df.init_time)], "sec": round(time.time() - t0, 1)}
    log_json(f"03_refs_{arch}_{year}.json", s)
    print({k: (v if k != "missing" else len(v)) for k, v in s.items()}, flush=True)
    return df

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("arch", nargs="+"); ap.add_argument("--years", default="2021-2025")
    a = ap.parse_args()
    y0, y1 = map(int, a.years.split("-")) if "-" in a.years else (int(a.years),) * 2
    for arch in a.arch:
        for y in range(y0, y1 + 1):
            build(arch, y)
