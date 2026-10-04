"""Shared paths and helpers for the ai-extremes-china pipeline."""
import os, time, json, pathlib, requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"
for d in (RAW, INTERIM, RESULTS, LOGS):
    d.mkdir(parents=True, exist_ok=True)

YEAR = int(os.environ.get("AIX_YEAR", "2023"))
LEADS_H = [int(x) for x in os.environ.get("AIX_LEADS", "24,72,120,168,240").split(",")]
N_STATIONS = int(os.environ.get("AIX_NSTATIONS", "50"))

MLWP = "https://noaa-oar-mlwp-data.s3.amazonaws.com"
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
AI_MODELS = {  # GFS-initialised AIWP archives (directory names)
    "graphcast": "GRAP_v100_GFS",
    "pangu": "PANG_v100_GFS",
    "fourcastnetv2": "FOUR_v200_GFS",
}
AI_MODELS_IFS = {  # IFS-initialised counterparts (initial-condition decomposition)
    "graphcast_ifs": "GRAP_v100_IFS",
    "pangu_ifs": "PANG_v100_IFS",
    "fourcastnetv2_ifs": "FOUR_v200_IFS",
}
ALL_AI = {**AI_MODELS, **AI_MODELS_IFS}
# variables streamed per AI archive (apcp only exists for GraphCast)
AI_VARS = {a: (["t2", "u10", "v10", "apcp"] if a.startswith("GRAP") else ["t2", "u10", "v10"]) for a in ALL_AI.values()}
NLEAD = 41                      # 0..240 h every 6 h
LEADS_ALL = list(range(0, 241, 6))
# Full-study period. ISD global-hourly stops on 2025-08-24 (ISD retired in favour of GHCNh), so the
# evaluation ends with inits on 2025-08-12 (10-day forecasts verify by 2025-08-22).
STUDY_START = "2021-01-01"
STUDY_END = "2025-08-12"
CLIM_YEARS = list(range(2011, 2021))   # independent climatology for thresholds
FC = ROOT / "data" / "fc"          # streamed station-point forecasts (npz per model-month)
REFS = ROOT / "data" / "refs"      # cached chunk tables
for d in (FC, REFS):
    d.mkdir(parents=True, exist_ok=True)

_session = None
def session():
    global _session
    if _session is None:
        s = requests.Session()
        a = requests.adapters.HTTPAdapter(pool_connections=32, pool_maxsize=32, max_retries=3)
        s.mount("https://", a)
        _session = s
    return _session

def get_range(url, start, length, timeout=60):
    """HTTP byte-range GET; returns bytes. Retries a few times."""
    hdr = {"Range": f"bytes={start}-{start + length - 1}"}
    for k in range(4):
        try:
            r = session().get(url, headers=hdr, timeout=timeout)
            if r.status_code in (200, 206):
                return r.content
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 ** k)
    raise RuntimeError(f"range GET failed: {url} {start}+{length}")

def log_json(name, obj):
    with open(LOGS / name, "w") as f:
        json.dump(obj, f, indent=2, default=str)

def peak_rss_mb():
    import resource
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
