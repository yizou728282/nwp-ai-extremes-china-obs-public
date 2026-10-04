import argparse, json, os, pathlib, sys
import pandas as pd, re
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from provenance import method_hash, require_current
ap = argparse.ArgumentParser()
ap.add_argument('--allow-legacy', action='store_true', help='Build a visibly labelled revision draft using historical results')
args = ap.parse_args()
os.chdir(ROOT / 'paper')
stale = []
for comparison in ('core', 'ic'):
    for variable in ('t2', 'tmax', 'tmin', 'ws', 'pr24'):
        try:
            require_current(ROOT / 'results/tables' / f'{comparison}_{variable}_meta.json')
        except RuntimeError as error:
            stale.append(str(error))
if stale and not args.allow_legacy:
    raise SystemExit('Corrected numerical results are unavailable. Run bash run_corrected.sh.\n'
                     'For a labelled historical revision draft only, use --allow-legacy.')
notice = (r'\newcommand{\AbstractRevisionNotice}{\textbf{Revision draft: numerical results and significance claims '
          r'are historical and await the corrected full rerun.} }' if stale else r'\newcommand{\AbstractRevisionNotice}{}')
revision = json.loads((ROOT / 'results/revision_status.json').read_text())
if not stale and revision.get('manuscript_numeric_review_method_sha256') != method_hash():
    notice = (r'\newcommand{\AbstractRevisionNotice}{\textbf{Revision draft: tables and figures have been '
              r'recomputed, but numerical prose and interpretation still require verification against those results.} }')
open('tables/revision_notice.tex', 'w').write(notice + '\n')
stations = pd.read_csv(ROOT / 'data/interim/stations_analysis.csv')
counts = {'Nstations': len(stations), 'NTmax': int(stations.use_tmax.sum()),
          'NTmin': int(stations.use_tmin.sum()), 'Nprecip': int(stations.use_pr24.sum())}
open('tables/counts.tex', 'w').write(''.join(f'\\newcommand{{\\{key}}}{{{value}}}\n' for key, value in counts.items()))
T='../results/tables/'
MN={'graphcast':'GraphCast Op.','pangu':'Pangu-Weather','fourcastnetv2':'FCN~v2-small'}
VN={'t2':'2-m temperature','tmax':'Daily $T_{\\max}$','tmin':'Daily $T_{\\min}$','ws':'10-m wind speed','pr24':'24-h precipitation'}
LEADS={'t2':['L24','L120','L240'],'ws':['L24','L120','L240'],'tmax':['L1','L5','L9'],'tmin':['L1','L5','L9'],'pr24':['L1','L5','L10']}
LL={'L24':'24 h','L120':'120 h','L240':'240 h','L1':'d1','L5':'d5','L9':'d9','L10':'d10'}
SC={'pcrps':'PCRPS','tw1':'twPCRPS$_1$','tw2':'twPCRPS$_2$','oos_idr':'OOS IDR','oos_emos':'OOS EMOS'}
k=pd.read_csv(T+'key_results_core.csv')
def fmt(s,brk):
    m=re.match(r'([+-][\d.]+)(\*?)\s*(\[.*\])?',str(s).strip())
    v,st,b=m.group(1),m.group(2),m.group(3)
    v=v.replace('-','$-$')
    out=(f'\\textbf{{{v}}}' if st else v)
    if brk and b: out+=f' {{\\scriptsize {b.replace("/","/")}}}'
    return out
rows=[]
for v in ['t2','tmax','tmin','ws','pr24']:
    L=LEADS[v]
    rows.append(f'\\midrule\n\\multicolumn{{5}}{{l}}{{\\textit{{{VN[v]}}} (leads: {", ".join(LL[x] for x in L)})}}\\\\')
    for s in ['pcrps','tw1','tw2','oos_idr']:
        for m in ['graphcast','pangu','fourcastnetv2']:
            r=k[(k['var']==v)&(k.score==s)&(k.model==m)]
            if r.empty: continue
            r=r.iloc[0]
            rows.append(f'{SC[s]} & {MN[m]} & '+' & '.join(fmt(r[x],s=='pcrps') for x in L)+r' \\')
txt='\n'.join(rows)+'\n'
i=txt.index('\\midrule\n\\multicolumn{5}{l}{\\textit{10-m')
open('tables/tab_skill_a.tex','w').write(txt[:i]); open('tables/tab_skill_b.tex','w').write(txt[i:])
# IC decomposition
d=pd.read_csv(T+'ic_decomposition_pcrps.csv')
LD={'t2':[24,120,240],'tmax':[1,5,9],'tmin':[1,5,9],'ws':[24,120,240],'pr24':[1,5,10]}
rows=[]
for v in ['t2','tmax','tmin','ws','pr24']:
    rows.append(f'\\midrule\n\\multicolumn{{4}}{{l}}{{\\textit{{{VN[v]}}}}}\\\\')
    for m in ['graphcast','pangu','fourcastnetv2']:
        x=d[(d['var']==v)&(d.arch==m)]
        if x.empty: continue
        cells=[]
        for l in LD[v]:
            y=x[x.lead==l].iloc[0]
            cells.append(f'{100*y.gfs_init:.1f} $\\rightarrow$ {100*y.ifs_init:.1f}'.replace('-','$-$'))
        rows.append(f'{MN[m]} & '+' & '.join(cells)+r' \\')
open('tables/tab_ic.tex','w').write('\n'.join(rows)+'\n')
# zones
rows=[]
Z=['All','Temperate','Arid','Cold','Highland','Tropical']
for v,l in [('t2',120),('tmax',5),('tmin',5),('ws',120),('pr24',5)]:
    b=pd.read_csv(T+f'core_{v}_bootstrap.csv'); b=b[(b.lead==l)&(b.ref=='gfs')]
    lab=f'{l} h' if v in('t2','ws') else f'day {l}'
    rows.append(f'\\midrule\n\\multicolumn{{8}}{{l}}{{\\textit{{{VN[v]}}}, {lab}}}\\\\')
    for s in ['pcrps','tw1']:
        for m in ['graphcast','pangu','fourcastnetv2']:
            x=b[(b.score==s)&(b.model==m)]
            if x.empty: continue
            c=[]
            for z in Z:
                y=x[x.zone==z].iloc[0]; t=f'{100*y.skill:+.1f}'.replace('-','$-$')
                c.append(f'\\textbf{{{t}}}' if (y.ci_lo>0 or y.ci_hi<0) else t)
            rows.append(f'{SC[s]} & {MN[m]} & '+' & '.join(c)+r' \\')
open('tables/tab_zone.tex','w').write('\n'.join(rows)+'\n')

# Full lead-by-lead heavy-precipitation comparison; overall OOS skill does not establish tail skill.
b = pd.read_csv(T + 'core_pr24_bootstrap.csv')
b = b[(b.zone == 'All') & (b.model == 'graphcast') & (b.ref == 'gfs')]
rows = []
for lead in range(1, 11):
    cells = []
    for score in ('tw1', 'oos_tw1'):
        r = b[(b.lead == lead) & (b.score == score)].iloc[0]
        cells.append(f'{100*r.skill:+.2f}')
        cells.append(f'[{100*r.ci_lo:+.2f}, {100*r.ci_hi:+.2f}]')
    rows.append(str(lead) + ' & ' + ' & '.join(cells) + r' \\')
open('tables/tab_oos_precip.tex', 'w').write('\n'.join(rows) + '\n')

# Include full six-hourly diagnostics only after corrected scores have been regenerated.
diurnal = ''
if not stale and pathlib.Path('figures/fig_diurnal_core.pdf').exists():
    diurnal = (r'\begin{figure}[t]\centering' + '\n' +
               r'\includegraphics[width=\textwidth]{fig_diurnal_core}' + '\n' +
               r'\caption{PCRPS skill at every six-hourly lead, separated by valid UTC hour (rows), for 2-m temperature and 10-m wind speed (columns). Models use the line styles and marker shapes in Fig.~\ref{fig:skill}. '
               r'Shading gives 95\% calendar-block bootstrap intervals.}\label{fig:diurnal}' + '\n' +
               r'\end{figure}' + '\n')
open('tables/diurnal_figure.tex', 'w').write(diurnal)

# main.tex = main_src.tex with the generated tables inlined (\input inside tabular is fragile)
if pathlib.Path('main_src.tex').exists():
    src = open('main_src.tex').read()
    open('main.tex', 'w').write(re.sub(r'\\input\{(tables/[^}]+)\}', lambda m: open(m.group(1)).read().rstrip('\n'), src))
else:
    print('Manuscript source is retained privately. Generated LaTeX result tables only.')
print('Built labelled historical draft.' if stale else 'Built paper tables from current corrected results.')
