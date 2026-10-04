"""Prevent legacy observation, score and table caches from masquerading as corrected results."""
import hashlib
import json
from datetime import datetime, timezone
from common import ROOT

PIPELINE_VERSION = "2026-10-03-verification-v2"
METHOD_FILES = ("02_obs.py", "02b_ghcnd.py", "05_build_dataset.py", "07_scores.py", "08_stats.py", "provenance.py")

def method_hash():
    digest = hashlib.sha256()
    for name in METHOD_FILES:
        digest.update(name.encode())
        digest.update((ROOT / "scripts" / name).read_bytes())
    return digest.hexdigest()

def write_metadata(path, **extra):
    path.parent.mkdir(parents=True, exist_ok=True)
    value = dict(pipeline_version=PIPELINE_VERSION, method_sha256=method_hash(),
                 generated_utc=datetime.now(timezone.utc).isoformat(), **extra)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return value

def require_current(path):
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise RuntimeError(f"Missing provenance for {path}; run ./run_corrected.sh") from error
    if value.get("pipeline_version") != PIPELINE_VERSION or value.get("method_sha256") != method_hash():
        raise RuntimeError(f"Stale verification output: {path}; run ./run_corrected.sh")
    return value
