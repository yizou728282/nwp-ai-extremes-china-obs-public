"""Publication figures (PDF + PNG, 300 dpi) under results/figures/.
Usage: python 09_figures.py [--set core]
Note: only coastlines are drawn (Natural Earth 1:50m); no national boundaries are plotted.
"""
import json, argparse, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import ROOT, RESULTS
from provenance import require_current

TB = RESULTS / "tables"; SC = RESULTS / "scores"; FIG = RESULTS / "figures"; FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.bbox": "tight", "pdf.fonttype": 42})
COL = {"graphcast": "#1b9e77", "pangu": "#d95f02", "fourcastnetv2": "#7570b3", "gfs": "k",
       "graphcast_ifs": "#1b9e77", "pangu_ifs": "#d95f02", "fourcastnetv2_ifs": "#7570b3"}
LAB = {"graphcast": "GraphCast Operational", "pangu": "Pangu-Weather", "fourcastnetv2": "FourCastNet v2-small", "gfs": "GFS"}
VLAB = {"t2": "2 m temperature", "tmax": "Daily Tmax", "tmin": "Daily Tmin", "ws": "10 m wind speed", "pr24": "24 h precipitation"}
TWLAB = {"t2": ("> q95", "< q05"), "ws": ("> q90", "> q95"), "tmax": ("> q90", "> q95"), "tmin": ("< q10", "< q05"), "pr24": ("> wet-day q90", "≥ 25 mm")}
ZONES = ["Tropical", "Temperate", "Arid", "Cold", "Highland"]
ZCOL = {"Tropical": "#e41a1c", "Temperate": "#4daf4a", "Arid": "#ff7f00", "Cold": "#377eb8", "Highland": "#984ea3"}

def save(fig, name):
    if name != "fig1_stations":
        try:
            for comparison in ("core", "ic"):
                for variable in ("t2", "tmax", "tmin", "ws", "pr24"):
                    require_current(TB / f"{comparison}_{variable}_meta.json")
        except RuntimeError:
            y = -0.40 if name.startswith("fig_ic_decomposition") else -0.02
            fig.text(0.5, y, "Historical numerical results - recomputation pending", ha="center", fontsize=7, color="0.35")
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300)
    plt.close(fig); print("saved", FIG / name)

def xlead(var, leads):
    return (np.asarray(leads) / 24.0) if var in ("t2", "ws") else np.asarray(leads)

def coast(ax):
    gj = json.load(open(ROOT / "data" / "aux" / "ne_50m_coastline.geojson"))
    for f in gj["features"]:
        g = f["geometry"]; parts = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
        for p in parts:
            p = np.asarray(p)
            if ((p[:, 0] > 65) & (p[:, 0] < 140) & (p[:, 1] > 10) & (p[:, 1] < 58)).any():
                ax.plot(p[:, 0], p[:, 1], lw=0.4, color="0.4")
    ax.set_xlim(72, 136); ax.set_ylim(17, 54.5); ax.set_aspect(1.25)

def fig_stations():
    st = pd.read_csv(ROOT / "data" / "interim" / "stations_analysis.csv", dtype={"station_id": str})
    fig, ax = plt.subplots(figsize=(5.2, 3.8)); coast(ax)
    for z in ZONES:
        s = st[st.zone == z]
        ax.scatter(s.lon, s.lat, s=9, c=ZCOL[z], label=f"{z} (n={len(s)})", edgecolors="none")
    ax.legend(loc="lower left", frameon=False, title="Köppen–Geiger zone", title_fontsize=7)
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    save(fig, "fig1_stations")

def load_boot(setname, var):
    f = TB / f"{setname}_{var}_bootstrap.csv"
    return pd.read_csv(f) if f.exists() else None

def fig_skill(setname, vars_, scores, name, title):
    models = [m for m in ["graphcast", "pangu", "fourcastnetv2"]]
    fig, axes = plt.subplots(len(scores), len(vars_), figsize=(1.9 * len(vars_) + 0.4, 1.7 * len(scores) + 0.4), squeeze=False, sharex="col")
    for j, v in enumerate(vars_):
        b = load_boot(setname, v)
        for i, (sc, slab) in enumerate(scores):
            ax = axes[i, j]; ax.axhline(0, color="k", lw=0.6)
            if b is None:
                continue
            for m in models:
                d = b[(b.zone == "All") & (b.model == m) & (b.ref == "gfs") & (b.score == sc)].sort_values("lead")
                if d.empty:
                    continue
                x = xlead(v, d.lead)
                ax.plot(x, 100 * d.skill, color=COL[m], lw=1.2, label=LAB[m])
                ax.fill_between(x, 100 * d.ci_lo, 100 * d.ci_hi, color=COL[m], alpha=0.2, lw=0)
            if i == 0:
                ax.set_title(VLAB[v])
            if j == 0:
                ax.set_ylabel(f"{slab}\nskill vs GFS (%)")
            if i == len(scores) - 1:
                ax.set_xlabel("Lead time (days)")
            if sc.startswith("tw") or sc.endswith("tw1"):
                k = 0 if sc in ("tw1", "oos_tw1", "oos_emos_tw1") else 1
                ax.text(0.97, 0.95, TWLAB[v][k], transform=ax.transAxes, ha="right", va="top", fontsize=7)
    h = [Line2D([], [], color=COL[m], lw=1.5, label=LAB[m]) for m in models]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    if title:
        fig.suptitle(title, y=1.06, fontsize=9)
    fig.tight_layout(); save(fig, name)

def fig_dm_map(setname, panels, model="graphcast", name=None):
    """panels: list of rows, each row a list of (var, score, lead). Filled = BH-FDR significant (q<0.05)."""
    st = pd.read_csv(ROOT / "data" / "interim" / "stations_analysis.csv", dtype={"station_id": str}).set_index("station_id")
    nr, nc = len(panels), max(len(r) for r in panels)
    fig, axes = plt.subplots(nr, nc, figsize=(7.0, 1.75 * nr), squeeze=False, gridspec_kw={"wspace": 0.04, "hspace": 0.22})
    for i, row in enumerate(panels):
        for j, (var, score, L) in enumerate(row):
            ax = axes[i, j]
            dm = pd.read_parquet(SC / f"{setname}_{var}_dm_fdr.parquet")
            d = dm[(dm.model == model) & (dm.ref == "gfs") & (dm.score == score) & (dm.lead == L)].copy()
            d["imp"] = -100 * d.mean_diff / d.ref_mean
            lon, lat = st.loc[d.station_id, "lon"].values, st.loc[d.station_id, "lat"].values
            coast(ax); sig = (d.q < 0.05).values
            ax.scatter(lon[~sig], lat[~sig], s=5, c="0.75", marker="o", edgecolors="none", zorder=2)
            sc = ax.scatter(lon[sig], lat[sig], c=d.imp.values[sig], s=8, cmap="RdBu", vmin=-30, vmax=30, edgecolors="k", linewidths=0.15, zorder=3)
            dd = L / 24 if var in ("t2", "ws") else L
            ax.set_title(f"({chr(97 + i * nc + j)}) {VLAB[var]}, day {dd:g}", fontsize=7.5, loc="left")
            better = (sig & (d.imp.values > 0)).mean() * 100; worse = (sig & (d.imp.values < 0)).mean() * 100
            ax.text(0.98, 0.03, f"better {better:.0f}%\nworse {worse:.0f}%", transform=ax.transAxes, fontsize=6.5, va="bottom", ha="right")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True); sp.set_linewidth(0.5)
    cb = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.015); cb.set_label(f"{LAB[model]} PCRPS improvement over GFS (%)")
    save(fig, name or f"fig_dm_map_{setname}_{model}")

def fig_zone_heatmap(setname, vars_, score, leadsel, name=None):
    st = pd.read_csv(ROOT / "data" / "interim" / "stations_analysis.csv", dtype={"station_id": str})
    SHORT = {"t2": "T2m", "tmax": "Tmax", "tmin": "Tmin", "ws": "Wind", "pr24": "Precip."}
    widths = [3 if v != "pr24" else 1.3 for v in vars_]
    # Each variable has a different eligible station set; never share its count labels.
    fig, axes = plt.subplots(1, len(vars_), figsize=(9.0, 2.8), sharey=False, squeeze=False, gridspec_kw={"width_ratios": widths})
    axes = axes[0]
    for ax, v in zip(axes, vars_):
        b = load_boot(setname, v)
        models = ["graphcast", "pangu", "fourcastnetv2"] if v != "pr24" else ["graphcast"]
        M = np.full((len(ZONES) + 1, len(models)), np.nan); S = np.zeros_like(M, bool)
        for j, m in enumerate(models):
            for i, z in enumerate(ZONES + ["All"]):
                d = b[(b.zone == z) & (b.model == m) & (b.ref == "gfs") & (b.score == score) & (b.lead == leadsel[v])]
                if len(d):
                    M[i, j] = 100 * d.skill.iloc[0]; S[i, j] = (d.ci_lo.iloc[0] > 0) or (d.ci_hi.iloc[0] < 0)
        im = ax.imshow(M, cmap="RdBu", vmin=-20, vmax=20, aspect="auto")
        for (i, j), val in np.ndenumerate(M):
            if np.isfinite(val):
                ax.text(j, i, f"{val:.0f}{'*' if S[i, j] else ''}", ha="center", va="center", fontsize=6.5)
        ax.set_xticks(range(len(models))); ax.set_xticklabels([{"graphcast": "GC", "pangu": "PW", "fourcastnetv2": "FCN2"}[m] for m in models])
        use = st if v in ("t2", "ws") else st[st[f"use_{v}"]]
        nz = use.zone.value_counts()
        labels = ([f"{z} ({nz.get(z, 0)})" for z in ZONES] + [f"All ({len(use)})"] if ax is axes[0]
                  else [f"n={nz.get(z, 0)}" for z in ZONES] + [f"n={len(use)}"])
        ax.set_yticks(range(len(ZONES) + 1)); ax.set_yticklabels(labels, fontsize=6.5)
        dd = leadsel[v] / 24 if v in ("t2", "ws") else leadsel[v]
        ax.set_title(f"{SHORT[v]}, day {dd:g}", fontsize=8)
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.subplots_adjust(wspace=0.60)
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.01); cb.set_label(f"{'PCRPS' if score == 'pcrps' else 'twPCRPS'} skill vs GFS (%)")
    save(fig, name or f"fig_zone_{setname}_{score}")

def fig_ic(vars_, score="pcrps"):
    rows = []
    for v in vars_:
        b = load_boot("ic", v)
        if b is None:
            continue
        for arch in ["graphcast", "pangu", "fourcastnetv2"]:
            for L in sorted(b.lead.unique()):
                g = b[(b.zone == "All") & (b.score == score) & (b.lead == L)]
                a = g[(g.model == arch) & (g.ref == "gfs")]; i = g[(g.model == arch + "_ifs") & (g.ref == "gfs")]
                if len(a) and len(i):
                    rows.append({"var": v, "arch": arch, "lead": L, "gfs_init": a.skill.iloc[0], "ifs_init": i.skill.iloc[0],
                                 "gfs_lo": a.ci_lo.iloc[0], "gfs_hi": a.ci_hi.iloc[0], "ifs_lo": i.ci_lo.iloc[0], "ifs_hi": i.ci_hi.iloc[0]})
    if not rows:
        return
    df = pd.DataFrame(rows); df.to_csv(TB / f"ic_decomposition_{score}.csv", index=False)
    vs = [v for v in vars_ if v in set(df["var"])]
    fig, axes = plt.subplots(1, len(vs), figsize=(1.9 * len(vs) + 0.4, 2.2), squeeze=False)
    for ax, v in zip(axes[0], vs):
        ax.axhline(0, color="k", lw=0.6)
        for arch in ["graphcast", "pangu", "fourcastnetv2"]:
            d = df[(df["var"] == v) & (df.arch == arch)].sort_values("lead")
            x = xlead(v, d.lead)
            ax.plot(x, 100 * d.gfs_init, color=COL[arch], lw=1.2, label=f"{LAB[arch]} (GFS init)")
            ax.plot(x, 100 * d.ifs_init, color=COL[arch], lw=1.2, ls="--", label=f"{LAB[arch]} (IFS init)")
            ax.fill_between(x, 100 * d.gfs_lo, 100 * d.gfs_hi, color=COL[arch], alpha=0.15, lw=0)
            ax.fill_between(x, 100 * d.ifs_lo, 100 * d.ifs_hi, color=COL[arch], alpha=0.08, lw=0)
        ax.set_title(VLAB[v]); ax.set_xlabel("Lead time (days)")
    axes[0, 0].set_ylabel(f"{score.upper()} skill vs GFS (%)")
    h = [Line2D([], [], color=COL[a], lw=1.3, label=LAB[a]) for a in ["graphcast", "pangu", "fourcastnetv2"]]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.15))
    h = [Line2D([], [], color="k", lw=1.2, label="GFS initial conditions"), Line2D([], [], color="k", lw=1.2, ls="--", label="IFS initial conditions")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.27))
    fig.tight_layout(); save(fig, f"fig_ic_decomposition_{score}")

def fig_is_vs_oos(setname, vars_):
    fig, axes = plt.subplots(1, len(vars_), figsize=(1.9 * len(vars_) + 0.4, 2.2), squeeze=False)
    for ax, v in zip(axes[0], vars_):
        b = load_boot(setname, v)
        if b is None:
            continue
        ax.axhline(0, color="k", lw=0.6)
        for m in ["graphcast", "pangu", "fourcastnetv2"]:
            for sc, ls in (("pcrps", "-"), ("oos_idr", "--"), ("oos_emos", ":")):
                d = b[(b.zone == "All") & (b.model == m) & (b.ref == "gfs") & (b.score == sc)].sort_values("lead")
                if len(d):
                    ax.plot(xlead(v, d.lead), 100 * d.skill, color=COL[m], ls=ls, lw=1.1)
        ax.set_title(VLAB[v]); ax.set_xlabel("Lead time (days)")
    axes[0, 0].set_ylabel("CRPS skill vs GFS (%)")
    h = [Line2D([], [], color="k", ls=ls, label=l) for ls, l in (("-", "PCRPS (in-sample IDR)"), ("--", "IDR out-of-sample"), (":", "EMOS out-of-sample"))]
    h += [Line2D([], [], color=COL[m], label=LAB[m]) for m in ["graphcast", "pangu", "fourcastnetv2"]]
    fig.legend(handles=h, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.08), fontsize=6)
    fig.tight_layout(); save(fig, f"fig_insample_vs_oos_{setname}")

def fig_diurnal(setname):
    """Display all six-hourly leads, separated by valid UTC hour."""
    fig, axes = plt.subplots(2, 4, figsize=(9.5, 4.0), sharex=True, sharey="row")
    for i, var in enumerate(("t2", "ws")):
        b = load_boot(setname, var)
        if b is None:
            continue
        for j, hour in enumerate((0, 6, 12, 18)):
            ax = axes[i, j]; ax.axhline(0, color="k", lw=0.6)
            for model in ("graphcast", "pangu", "fourcastnetv2"):
                d = b[(b.zone == "All") & (b.ref == "gfs") & (b.model == model)
                      & (b.score == "pcrps") & (b.lead % 24 == hour)].sort_values("lead")
                if not d.empty:
                    x = xlead(var, d.lead)
                    ax.plot(x, 100 * d.skill, color=COL[model], lw=1.1)
                    ax.fill_between(x, 100 * d.ci_lo, 100 * d.ci_hi, color=COL[model], alpha=0.15, lw=0)
            if i == 0:
                ax.set_title(f"Valid hour {hour:02d} UTC")
            else:
                ax.set_xlabel("Lead time (days)")
            if j == 0:
                ax.set_ylabel(f"{VLAB[var]}\nPCRPS skill vs GFS (%)")
    h = [Line2D([], [], color=COL[m], label=LAB[m]) for m in ("graphcast", "pangu", "fourcastnetv2")]
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(); save(fig, f"fig_diurnal_{setname}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--set", default="core"); a = ap.parse_args()
    V = ["t2", "tmax", "tmin", "ws", "pr24"]
    fig_stations()
    fig_skill(a.set, V, [("pcrps", "PCRPS"), ("tw1", "twPCRPS (tail 1)"), ("tw2", "twPCRPS (tail 2)")], f"fig2_pcrps_skill_{a.set}", None)
    fig_skill(a.set, V, [("oos_idr", "IDR CRPS"), ("oos_tw1", "IDR twCRPS (tail 1)")], f"fig3_oos_skill_{a.set}", None)
    fig_is_vs_oos(a.set, ["t2", "tmax", "tmin", "ws", "pr24"])
    fig_diurnal(a.set)
    fig_zone_heatmap(a.set, V, "pcrps", {"t2": 120, "ws": 120, "tmax": 5, "tmin": 5, "pr24": 5})
    fig_zone_heatmap(a.set, V, "tw1", {"t2": 120, "ws": 120, "tmax": 5, "tmin": 5, "pr24": 5})
    fig_dm_map(a.set, [[("t2", "pcrps", 24), ("t2", "pcrps", 120), ("t2", "pcrps", 240)],
                       [("tmin", "pcrps", 1), ("tmin", "pcrps", 5), ("tmin", "pcrps", 9)],
                       [("pr24", "pcrps", 1), ("pr24", "pcrps", 5), ("pr24", "pcrps", 10)]])
    fig_ic(V)
