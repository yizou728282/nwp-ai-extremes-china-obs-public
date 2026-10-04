"""GHCN-Daily TMAX/TMIN/PRCP for analysis stations that have a GHCN-D record (CHM000<WMO>), 2011-2025.
Used (i) for 2011-2020 climatological thresholds of Tmax/Tmin/precipitation (ISD KA groups only exist
from late 2019 and 24-h AA totals are sparse before 2019), (ii) as a cross-check of the ISD daily values.
Output: data/interim/ghcnd_daily.parquet (station_id, date, tmax, tmin, prcp), logs/02b_ghcnd_check.csv
"""
import io, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from common import INTERIM, LOGS, session
obsmod = __import__("importlib").import_module("02_obs")

st = pd.read_csv(INTERIM / "stations_analysis.csv", dtype={"station_id": str})
def one(sid):
    r = session().get(f"https://noaa-ghcn-pds.s3.amazonaws.com/csv/by_station/CHM000{sid[:5]}.csv", timeout=120)
    if r.status_code != 200:
        return None
    g = pd.read_csv(io.StringIO(r.text), dtype=str, usecols=["DATE", "ELEMENT", "DATA_VALUE", "Q_FLAG"])
    g = g[g.ELEMENT.isin(["TMAX", "TMIN", "PRCP"]) & g.Q_FLAG.isna()]
    g["date"] = pd.to_datetime(g.DATE); g = g[g.date >= "2011-01-01"]
    g["v"] = pd.to_numeric(g.DATA_VALUE) / 10.0
    w = g.pivot_table(index="date", columns="ELEMENT", values="v").reset_index()
    w.columns = [c.lower() for c in w.columns]
    w.insert(0, "station_id", sid)
    return w.rename(columns={"prcp": "prcp"})
with ThreadPoolExecutor(6) as ex:
    out = [o for o in ex.map(one, st.station_id) if o is not None]
gd = pd.concat(out, ignore_index=True).astype({"tmax": "float32", "tmin": "float32", "prcp": "float32"}, errors="ignore")
gd.to_parquet(INTERIM / "ghcnd_daily.parquet", index=False)
gd["y"] = gd.date.dt.year
print("stations with GHCN-D:", gd.station_id.nunique())
print(gd.groupby("y")[["tmax", "tmin", "prcp"]].count())
# agreement with ISD (2020-2023) for lags -1/0/+1 days
isd = obsmod.load("daily", range(2020, 2024), set(gd.station_id)).groupby(["station_id", "date"], as_index=False).max()
rows = []
for lag in (-1, 0, 1):
    g2 = gd.assign(date=gd.date + pd.Timedelta(days=lag))
    m = g2.merge(isd, on=["station_id", "date"], suffixes=("_g", "_i"))
    rows.append({"lag_days": lag, "n": len(m), "tmax_mad": (m.tmax_g - m.tmax_i).abs().median(), "tmax_exact": ((m.tmax_g - m.tmax_i).abs() < 0.05).mean(),
                 "tmin_mad": (m.tmin_g - m.tmin_i).abs().median(), "tmin_exact": ((m.tmin_g - m.tmin_i).abs() < 0.05).mean(),
                 "pr_corr": m[["prcp", "pr24"]].dropna().corr().iloc[0, 1], "pr_exact": ((m.prcp - m.pr24).abs() < 0.05)[m.pr24 > 0].mean()})
print(pd.DataFrame(rows).round(3)); pd.DataFrame(rows).to_csv(LOGS / "02b_ghcnd_check.csv", index=False)
