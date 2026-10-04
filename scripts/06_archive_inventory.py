"""Inventory of the AIWP archive for the full-study period: 00 UTC inits per year present in the
kerchunk index (with a non-null t2 reference at +24 h) vs. NetCDF files actually on S3.
Output: results/archive_inventory.csv
"""
import io, json, re, numpy as np, pandas as pd, xml.etree.ElementTree as ET
from common import MLWP, session, get_range, RESULTS

ARCH = ["GRAP_v100_GFS", "PANG_v100_GFS", "FOUR_v200_GFS", "GRAP_v100_IFS", "PANG_v100_IFS", "FOUR_v200_IFS"]
YEARS = range(2020, 2027)

def s3_keys(prefix):
    keys, token = [], None
    while True:
        q = f"{MLWP}/?list-type=2&prefix={prefix}" + (f"&continuation-token={__import__('urllib.parse').parse.quote(token, safe='')}" if token else "")
        x = ET.fromstring(session().get(q, timeout=60).content)
        ns = {"s": x.tag.split("}")[0].strip("{")}
        keys += [(c.find("s:Key", ns).text, int(c.find("s:Size", ns).text)) for c in x.findall("s:Contents", ns)]
        t = x.find("s:NextContinuationToken", ns)
        if t is None:
            return keys
        token = t.text

rows = []
for a in ARCH:
    base = f"{MLWP}/parquet/{a}_combined_all.parq"
    meta = json.loads(session().get(f"{base}/.zmetadata", timeout=60).text)
    rec = meta["record_size"]
    d = pd.read_parquet(io.BytesIO(session().get(f"{base}/init_time/refs.0.parq", timeout=60).content)).iloc[0]
    b = d["raw"] if d["raw"] is not None else get_range(d["path"].replace("s3://noaa-oar-mlwp-data", MLWP), int(d["offset"]), int(d["size"]))
    inits = pd.to_datetime(np.frombuffer(b, "<M8[ns]"))
    n = len(inits); nrec = -(-n * 41 // rec)
    refs = pd.concat([pd.read_parquet(io.BytesIO(session().get(f"{base}/t2/refs.{r}.parq", timeout=120).content)) for r in range(nrec)], ignore_index=True)
    ok24 = np.array([isinstance(refs["path"].iloc[i * 41 + 4], str) and refs["size"].iloc[i * 41 + 4] > 0 for i in range(n)])
    for y in YEARS:
        m = (inits.year == y) & (inits.hour == 0)
        keys = s3_keys(f"{a}/{y}/")
        k00 = [k for k, s in keys if re.search(r"_\d{8}00_f000", k)]
        sizes = [s for k, s in keys if re.search(r"_\d{8}00_f000", k)]
        ndays = 366 if y % 4 == 0 else 365
        rows.append({"archive": a, "year": y, "days": ndays, "index_00z_inits": int(m.sum()), "index_00z_with_t2_ref": int((m & ok24).sum()),
                     "s3_00z_files": len(k00), "s3_00z_GB_total": round(sum(sizes) / 1e9, 1)})
        print(rows[-1])
pd.DataFrame(rows).to_csv(RESULTS / "archive_inventory.csv", index=False)
