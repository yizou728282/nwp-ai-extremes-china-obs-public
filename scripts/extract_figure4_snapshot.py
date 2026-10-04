"""Read plotted station positions, significance and RGB fills from the original vector PDF.

This extracts the existing figure's display data. It does not reconstruct numeric
losses, effect magnitudes, p values or q values from colour.
"""
from pathlib import Path
import hashlib, json, re
import numpy as np
import pandas as pd
from pypdf import PdfReader
from pypdf.generic import ContentStream
import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / 'work/accessibility_repo').exists():
    ROOT = ROOT / 'work/accessibility_repo'
SOURCE = ROOT / 'paper/figure_sources/fig4_original.pdf'
if not SOURCE.exists():
    SOURCE = ROOT / 'paper/figures/fig_dm_map_core_graphcast.pdf'

def compose(old, new):
    a,b,c,d,e,f = old
    g,h,i,j,k,l = new
    return (a*g+c*h,b*g+d*h,a*i+c*j,b*i+d*j,a*k+c*l+e,b*k+d*l+f)

reader = PdfReader(SOURCE)
page = reader.pages[0]
height = float(page.mediabox.height)
state = {'matrix': (1.,0.,0.,1.,0.,0.), 'fill': (0.,0.,0.)}
stack = []
points = []
for values, op in ContentStream(page.get_contents(), reader).operations:
    if op == b'q':
        stack.append(state.copy())
    elif op == b'Q':
        state = stack.pop()
    elif op == b'cm':
        state['matrix'] = compose(state['matrix'], tuple(float(v) for v in values))
    elif op == b'rg':
        state['fill'] = tuple(float(v) for v in values)
    elif op == b'g':
        state['fill'] = (float(values[0]),) * 3
    elif op == b'Do' and re.fullmatch(r'/P\d+', str(values[0])):
        marker = int(str(values[0])[2:])
        matrix = state['matrix']
        points.append({'panel': marker // 2, 'significant': marker % 2 == 1,
                       'x': matrix[4], 'top': height-matrix[5], 'rgb': state['fill']})

with pdfplumber.open(SOURCE) as doc:
    lines = [x for x in doc.pages[0].lines if abs(x['linewidth'] - .5) < .01]
    frames = []
    for offset in range(0,36,4):
        group = lines[offset:offset+4]
        frames.append((min(x['x0'] for x in group),min(x['top'] for x in group),
                       max(x['x1'] for x in group),max(x['bottom'] for x in group)))
    old_labels = [(int(a),int(b)) for a,b in re.findall(r'better\s+(\d+)%.*?worse\s+(\d+)%',doc.pages[0].extract_text().replace('\n',' '))]

stations = pd.read_csv(ROOT/'data/interim/stations_analysis.csv',dtype={'station_id':str})
panels = [('t2',24),('t2',120),('t2',240),('tmin',1),('tmin',5),('tmin',9),('pr24',1),('pr24',5),('pr24',10)]
rows = []
maximum_position_error = 0.
for point in points:
    panel = point['panel']
    x0,y0,x1,y1 = frames[panel]
    lon = 72 + 64*(point['x']-x0)/(x1-x0)
    lat = 54.5 - 37.5*(point['top']-y0)/(y1-y0)
    distances = (stations.lon-lon)**2 + (stations.lat-lat)**2
    idx = distances.idxmin()
    error = float(np.sqrt(distances.loc[idx]))
    maximum_position_error = max(maximum_position_error,error)
    assert error < .001, (panel,lon,lat,error)
    station = stations.loc[idx]
    rgb = point['rgb']
    sign = (1 if rgb[2] > rgb[0] else -1) if point['significant'] else 0
    variable, lead = panels[panel]
    rows.append({'panel':panel,'var':variable,'lead':lead,'station_id':station.station_id,
                 'lon':station.lon,'lat':station.lat,'plot_significant':point['significant'],
                 'plot_sign':sign,'fill_r':rgb[0],'fill_g':rgb[1],'fill_b':rgb[2]})

df = pd.DataFrame(rows)
assert len(df)==2289
comparisons=[]
for panel,(variable,lead) in enumerate(panels):
    d=df[df.panel==panel]
    expected=334 if variable=='t2' else (180 if variable=='tmin' else 249)
    assert len(d)==expected and d.station_id.nunique()==expected
    summary=pd.read_csv(ROOT/f'results/tables/core_{variable}_dm_fdr.csv')
    record=summary[(summary.model=='graphcast') & (summary.ref=='gfs') & (summary.score=='pcrps') & (summary.lead==lead)].iloc[0]
    for sign,key in [(1,'frac_better_fdr05'),(-1,'frac_worse_fdr05')]:
        fraction=float((d.plot_sign==sign).mean())
        assert np.isclose(round(fraction,4),record[key],rtol=0,atol=1e-12),(panel,key,fraction,record[key])
    comparisons.append({'panel':panel,'stations':expected,'better':int((d.plot_sign==1).sum()),'worse':int((d.plot_sign==-1).sum())})

dest=ROOT/'data/interim/fig4_plot_snapshot.csv'
df.to_csv(dest,index=False,lineterminator='\n')
provenance={'source':'existing vector figure, with station identities matched to the unchanged station list',
            'source_pdf_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'snapshot_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
            'maximum_position_match_error_degrees':maximum_position_error,
            'all_18_station_fraction_checks_passed':True,'summary_fraction_precision':4,'panels':comparisons,
            'scope':'positions, plotted significance/sign, and original RGB fills only; no numeric effect sizes or p/q values reconstructed'}
(ROOT/'data/interim/fig4_plot_snapshot_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
print(json.dumps(provenance,indent=2))
