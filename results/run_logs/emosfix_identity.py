"""Exit 0 iff the re-scored station/day/DM outputs equal the backed-up ones; then refresh the 08 table provenance."""
import json, sys, pathlib
import numpy as np, pandas as pd
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root / "scripts"))
from provenance import method_hash
S, V, B = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3])
SC = root / "results" / "scores"; TB = root / "results" / "tables"
ok = True
for f in sorted(B.glob(f"{S}_{V}_*")):
    new = SC / f.name
    if f.suffix == ".parquet":
        a, b = pd.read_parquet(f), pd.read_parquet(new)
        keys = [c for c in ("station_id", "lead", "model", "ref", "score") if c in a.columns]
        if keys:
            a = a.sort_values(keys).reset_index(drop=True); b = b.sort_values(keys).reset_index(drop=True)
        try: pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-12, atol=0)
        except AssertionError as e: print("DIFF", f.name, str(e)[:300]); ok = False
    elif f.suffix == ".npz":
        a, b = np.load(f, allow_pickle=True), np.load(new, allow_pickle=True)
        if sorted(a.files) != sorted(b.files) or not all((np.array_equal(a[k], b[k], equal_nan=True) if a[k].dtype.kind not in 'fc' else np.allclose(a[k], b[k], rtol=1e-12, atol=0, equal_nan=True)) for k in a.files):
            print("DIFF", f.name); ok = False
    print("checked", f.name)
if not ok: sys.exit(1)
m = TB / f"{S}_{V}_meta.json"; d = json.loads(m.read_text())
d.update(previous_method_sha256=d.get("method_sha256"), method_sha256=method_hash(),
         restamp_reason="scores re-run under new method hash were identical to the inputs of these tables")
m.write_text(json.dumps(d, indent=2) + "\n"); print("restamped", m.name)
