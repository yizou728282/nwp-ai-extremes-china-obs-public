"""Compact key-result tables (markdown + csv) for the report/manuscript: results/tables/key_results_<set>.md"""
import sys, numpy as np, pandas as pd
from common import RESULTS
TB = RESULTS / "tables"
S = sys.argv[1] if len(sys.argv) > 1 else "core"
LEADS = {"t2": [24, 72, 120, 168, 240], "ws": [24, 72, 120, 168, 240], "tmax": [1, 3, 5, 7, 9], "tmin": [1, 3, 5, 7, 9], "pr24": [1, 3, 5, 7, 10]}
out = []
for v, Ls in LEADS.items():
    f = TB / f"{S}_{v}_bootstrap.csv"
    if not f.exists():
        continue
    b = pd.read_csv(f); b = b[(b.zone == "All") & (b.ref == "gfs") & b.lead.isin(Ls)]
    dm = pd.read_csv(TB / f"{S}_{v}_dm_fdr.csv"); dm = dm[dm.ref == "gfs"]
    for sc in ["pcrps", "tw1", "tw2", "oos_idr", "oos_tw1", "oos_tw2", "oos_emos", "oos_emos_tw1"]:
        for m in sorted(b.model.unique()):
            d = b[(b.score == sc) & (b.model == m)].set_index("lead")
            if d.empty:
                continue
            row = {"var": v, "score": sc, "model": m}
            for L in Ls:
                if L in d.index:
                    r = d.loc[L]; star = "*" if (r.ci_lo > 0 or r.ci_hi < 0) else ""
                    q = dm[(dm.model == m) & (dm.score == sc) & (dm.lead == L)]
                    fb = f" [{100*q.frac_better_fdr05.iloc[0]:.0f}/{100*q.frac_worse_fdr05.iloc[0]:.0f}]" if len(q) else ""
                    row[f"L{L}"] = f"{100*r.skill:+.1f}{star}{fb}"
            out.append(row)
df = pd.DataFrame(out); df.to_csv(TB / f"key_results_{S}.csv", index=False)
with open(TB / f"key_results_{S}.md", "w") as fh:
    fh.write(f"# Key results ({S}): skill vs GFS in % (* = 95% block-bootstrap CI excludes 0) "
             "[% stations significantly better / worse, DM test, BH-FDR 5%]\n\n")
    for v in df["var"].unique():
        fh.write(f"## {v}\n\n" + df[df["var"] == v].drop(columns="var").dropna(axis=1, how="all").to_markdown(index=False) + "\n\n")
print(open(TB / f"key_results_{S}.md").read()[:3000])
