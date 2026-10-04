"""Station selection for the full study.

superset (used for forecast extraction; cheap because every streamed field is global):
    mainland China ISD stations (WMO block 50-59, 73-136E, 17-54N), active 2021-2025.
analysis set (decided from the ISD inventory, refined later by 02 coverage checks):
    * every year 2021-2024 and Jan-Aug 2025: every month >= MIN_PER_DAY obs/day (00 UTC reliably present)
    * climatology 2011-2020: >= 8 years with an annual mean of >= MIN_PER_DAY obs/day
    * de-duplicated: of stations closer than 10 km, keep the one with more reports
Each station gets a Koppen-Geiger class (Beck et al. 2023, V2 1991-2020, 0.1 deg) and a broad zone.
Outputs: data/interim/stations_superset.csv, data/interim/stations.csv
"""
import os, numpy as np, pandas as pd, tifffile
from common import RAW, INTERIM, ROOT, session

MIN_PER_DAY = float(os.environ.get("AIX_MIN_PER_DAY", "4"))
MONTHS = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]

def fetch(name):
    p = RAW / name
    if not p.exists():
        r = session().get(f"https://www.ncei.noaa.gov/pub/data/noaa/{name}", timeout=120); r.raise_for_status(); p.write_bytes(r.content)
    return p

def koppen(lat, lon):
    aux = ROOT / "data" / "aux"
    img = tifffile.imread(aux / "koppen_geiger_0p1.tif")          # 1800 x 3600, 90N..90S, 180W..180E
    leg = {}
    for line in open(aux / "legend.txt", encoding="utf-8", errors="ignore"):
        parts = line.strip().split(":")
        if len(parts) >= 2 and parts[0].strip().isdigit():
            leg[int(parts[0])] = parts[1].strip().split()[0]
    i = np.clip(((90 - lat) / 0.1).astype(int), 0, img.shape[0] - 1)
    j = np.clip(((lon + 180) / 0.1).astype(int), 0, img.shape[1] - 1)
    codes = img[i, j]
    # stations on a coast pixel classified as water (0): search a 3x3 neighbourhood
    for k in np.where(codes == 0)[0]:
        win = img[max(i[k]-2, 0):i[k]+3, max(j[k]-2, 0):j[k]+3].ravel()
        win = win[win > 0]
        codes[k] = np.bincount(win).argmax() if len(win) else 0
    cls = np.array([leg.get(int(c), "NA") for c in codes])
    return cls

def zone(cls, elev):
    """Broad zones for stratification: Tropical (A), Arid (B), Temperate (C), Cold (D), Highland/Polar (E or elev>3000 m)."""
    z = np.array([{"A": "Tropical", "B": "Arid", "C": "Temperate", "D": "Cold", "E": "Highland"}.get(c[0], "NA") for c in cls], dtype=object)
    z[np.asarray(elev) > 3000] = "Highland"
    return z

hist = pd.read_csv(fetch("isd-history.csv"), dtype=str)
inv = pd.read_csv(fetch("isd-inventory.csv"), dtype={"USAF": str, "WBAN": str})
h = hist[(hist["CTRY"] == "CH") & hist["USAF"].str.match(r"^5[0-9]{4}0$")].copy()
for c in ("LAT", "LON", "ELEV(M)"):
    h[c] = pd.to_numeric(h[c], errors="coerce")
h = h.dropna(subset=["LAT", "LON"])
h = h[h["LON"].between(73, 136) & h["LAT"].between(17, 54) & (h["END"] >= "20250101") & (h["BEGIN"] <= "20201231")]
h["station_id"] = h["USAF"] + h["WBAN"]
h = h.drop_duplicates("station_id")
sup = h.rename(columns={"STATION NAME": "name", "LAT": "lat", "LON": "lon", "ELEV(M)": "elev"})[["station_id", "USAF", "WBAN", "name", "lat", "lon", "elev"]]
sup = sup.reset_index(drop=True)
sup["koppen"] = koppen(sup.lat.values, sup.lon.values)
sup["zone"] = zone(sup.koppen.values, sup.elev.fillna(0).values)
sup.to_csv(INTERIM / "stations_superset.csv", index=False)

inv["station_id"] = inv["USAF"] + inv["WBAN"]
inv = inv[inv.station_id.isin(sup.station_id)]
# NOTE: NCEI's isd-inventory.csv now only covers 2020 onward, so the 2011-2020 climatology criterion
# (>= 8 years with >= 80% of days having a report) is applied after parsing, in 02_obs.py --summary.
# Here: annual mean >= MIN_PER_DAY obs/day every year 2021-2024 and Jan-Jul 2025 (network-wide outages,
# e.g. Nov 2024, make a per-month criterion too strict), and station history beginning before 2011.
ok_eval = pd.Series(True, index=sup.station_id)
for y in range(2021, 2026):
    m = MONTHS if y < 2025 else MONTHS[:7]
    days = 365 if y < 2025 else 212
    iy = inv[inv.YEAR == y].set_index("station_id")[m].sum(axis=1)
    ok_eval &= (iy.reindex(sup.station_id).fillna(0) >= MIN_PER_DAY * days).values
begin = h.set_index("station_id")["BEGIN"]
ok_clim = (begin.reindex(sup.station_id) <= "20110101")
tot = inv[inv.YEAR.between(2021, 2024)].groupby("station_id")[MONTHS].sum().sum(axis=1).reindex(sup.station_id).fillna(0)
sel = sup[(ok_eval & ok_clim).values].copy()
sel["reports_2021_2024"] = tot.reindex(sel.station_id).values
# de-duplicate stations closer than 10 km
sel = sel.sort_values("reports_2021_2024", ascending=False).reset_index(drop=True)
keep = []
for k, r in sel.iterrows():
    d = 6371 * np.hypot(np.radians(sel.loc[keep, "lat"] - r.lat), np.radians(sel.loc[keep, "lon"] - r.lon) * np.cos(np.radians(r.lat))) if keep else np.array([])
    if not (np.asarray(d) < 10).any():
        keep.append(k)
sel = sel.loc[keep].sort_values("station_id").reset_index(drop=True)
sel.to_csv(INTERIM / "stations.csv", index=False)
print(f"superset {len(sup)}; eval-coverage {int(ok_eval.sum())}; clim-coverage {int(ok_clim.sum())}; analysis set {len(sel)}")
print(sel.zone.value_counts().to_dict()); print(sel.koppen.value_counts().head(10).to_dict())
