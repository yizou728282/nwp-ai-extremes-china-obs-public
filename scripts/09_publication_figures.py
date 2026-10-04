"""Redraw the seven manuscript figures without changing verification results.

Layouts are designed for a 5.15-inch text column, with 9.5-point minimum
source text, redundant line/marker coding and explicit signed heatmap labels.
Figure 4 uses the documented display-only snapshot of the original vector PDF.
"""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from common import ROOT, RESULTS
from provenance import require_current

TB = RESULTS / 'tables'
FIG = RESULTS / 'figures'
FIG.mkdir(parents=True, exist_ok=True)
PAPER = ROOT / 'paper/figures'
PAPER.mkdir(parents=True, exist_ok=True)
W = 5.15
FS = 9.5
plt.rcParams.update({'font.size': FS, 'axes.titlesize': FS,
    'axes.labelsize': FS, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
    'legend.fontsize': FS, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .6, 'savefig.bbox': None})
MODELS = ['graphcast', 'pangu', 'fourcastnetv2']
COL = dict(zip(MODELS, ['#087f5b', '#c44a00', '#6653a2']))
MARK = dict(zip(MODELS, ['o', 's', '^']))
LS = dict(zip(MODELS, ['-', '--', ':']))
LAB = dict(zip(MODELS, ['GraphCast Operational', 'Pangu-Weather', 'FourCastNet v2-small']))
V = ['t2', 'tmax', 'tmin', 'ws', 'pr24']
VL = dict(zip(V, ['2 m temperature', 'Daily Tmax', 'Daily Tmin', '10 m wind speed', '24 h precipitation']))
TAIL = {'t2': ['> q95', '< q05'], 'tmax': ['> q90', '> q95'],
    'tmin': ['< q10', '< q05'], 'ws': ['> q90', '> q95'],
    'pr24': ['> wet-day q90', '>= 25 mm']}
ZONES = ['Tropical', 'Temperate', 'Arid', 'Cold', 'Highland']
ZCOL = ['#c23a3a', '#32823b', '#d37513', '#306db5', '#8552a1']
checks = []

def save(fig, name):
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    omitted = set()
    for ax in fig.axes:
        for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            lo,hi = sorted(limits)
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lo-1e-9 <= tick.get_loc() <= hi+1e-9:
                    omitted.update([id(tick.label1),id(tick.label2)])
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip() and id(t) not in omitted]
    for t in texts:
        b = t.get_window_extent(renderer)
        assert b.x0 >= -1 and b.y0 >= -1 and b.x1 <= fig.bbox.width+1 and b.y1 <= fig.bbox.height+1, (name,t.get_text(),tuple(b.bounds))
    assert all(t.get_fontsize() >= FS for t in texts), name
    for ext in ('pdf', 'png'):
        fig.savefig(FIG / f'{name}.{ext}', dpi=300)
    shutil.copyfile(FIG / f'{name}.pdf', PAPER / f'{name}.pdf')
    checks.append({'figure': name, 'width_inches': W,
        'height_inches': fig.get_figheight(), 'minimum_source_font_pt': min(t.get_fontsize() for t in texts),
        'maximum_source_font_pt': max(t.get_fontsize() for t in texts),
        'pdf_sha256': hashlib.sha256((PAPER / f'{name}.pdf').read_bytes()).hexdigest()})
    plt.close(fig)
    print('saved', name)

def data(kind, var):
    return pd.read_csv(TB / f'{kind}_{var}_bootstrap.csv')

def select(b, model, score, zone='All', ref='gfs'):
    return b[(b.zone == zone) & (b.model == model) & (b.ref == ref) & (b.score == score)].sort_values('lead')

def xlead(var, leads):
    return np.asarray(leads) / 24 if var in ('t2', 'ws') else np.asarray(leads)

def decorate(ax, xlabel=False, ylabel=False):
    ax.axhline(0, lw=.7, color='0.15', zorder=1)
    ax.set_xlim(0, 10.2)
    ax.set_xticks([0, 5, 10])
    ax.grid(axis='y', color='0.85', lw=.4)
    if xlabel: ax.set_xlabel('Lead (days)')
    if ylabel: ax.set_ylabel('Skill (%)')

def draw(ax, d, var, m, ls=None, band=True, open_marker=False):
    if d.empty: return
    x = xlead(var, d.lead)
    step = max(1, len(d) // 5)
    ax.plot(x, 100*d.skill, color=COL[m], ls=ls or LS[m], lw=1.2,
        marker=MARK[m], markersize=3.2, markevery=step,
        markerfacecolor='white' if open_marker else COL[m], markeredgewidth=.6)
    if band:
        ax.fill_between(x, 100*d.ci_lo, 100*d.ci_hi, color=COL[m], alpha=.12, lw=0)

def model_handles(styles=True):
    return [Line2D([], [], color=COL[m], ls=LS[m] if styles else '-',
        marker=MARK[m], markersize=4, lw=1.2, label=LAB[m]) for m in MODELS]

def coast(ax):
    gj = json.loads((ROOT/'data/aux/ne_50m_coastline.geojson').read_text())
    for f in gj['features']:
        g = f['geometry']
        for p in g['coordinates'] if g['type'] == 'MultiLineString' else [g['coordinates']]:
            p = np.asarray(p)
            if ((p[:,0]>65) & (p[:,0]<140) & (p[:,1]>10) & (p[:,1]<58)).any():
                ax.plot(p[:,0],p[:,1],color='0.5',lw=.4,zorder=0)
    ax.set_xlim(72,136); ax.set_ylim(17,54.5); ax.set_aspect(1.25)

def stations():
    st = pd.read_csv(ROOT/'data/interim/stations_analysis.csv', dtype={'station_id': str})
    fig, ax = plt.subplots(figsize=(W,4.7))
    fig.subplots_adjust(left=.12,right=.98,bottom=.31,top=.98)
    coast(ax)
    for z,c,m in zip(ZONES,ZCOL,['o','s','^','D','v']):
        d = st[st.zone == z]
        ax.scatter(d.lon,d.lat,s=17,c=c,marker=m,edgecolors='0.15',linewidths=.35,
            label=f'{z} (n={len(d)})')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,frameon=False,
        title='Koppen-Geiger zone',title_fontsize=FS,
        borderpad=.2,labelspacing=.4,handletextpad=.3,columnspacing=.8)
    ax.set_xlabel('Longitude (degrees E)'); ax.set_ylabel('Latitude (degrees N)')
    save(fig,'fig1_stations')

def skill():
    fig, axes = plt.subplots(5,3,figsize=(W,6.45))
    fig.subplots_adjust(left=.135,right=.98,bottom=.075,top=.83,wspace=.40,hspace=.67)
    fig.legend(handles=model_handles(),loc='upper center',ncol=1,frameon=False,
        bbox_to_anchor=(.52,1),labelspacing=.15,handlelength=2.7)
    for i,v in enumerate(V):
        b = data('core',v)
        for j,score in enumerate(['pcrps','tw1','tw2']):
            ax = axes[i,j]
            for m in MODELS: draw(ax,select(b,m,score),v,m)
            decorate(ax,xlabel=i==4,ylabel=j==0)
            title = VL[v] if j==0 else TAIL[v][j-1]
            if i==0: title = ['PCRPS','twPCRPS: tail 1','twPCRPS: tail 2'][j]+'\n'+title
            ax.set_title(title,pad=3)
    save(fig,'fig2_pcrps_skill_core')

def diurnal():
    fig,axes = plt.subplots(4,2,figsize=(W,5.5),sharex=True,sharey='col')
    fig.subplots_adjust(left=.135,right=.98,bottom=.085,top=.81,wspace=.30,hspace=.60)
    fig.legend(handles=model_handles(),loc='upper center',ncol=1,frameon=False,
        bbox_to_anchor=(.52,1),labelspacing=.15)
    for i,h in enumerate([0,6,12,18]):
        for j,v in enumerate(['t2','ws']):
            ax=axes[i,j]; b=data('core',v)
            for m in MODELS:
                d=select(b,m,'pcrps'); draw(ax,d[d.lead%24==h],v,m)
            decorate(ax,xlabel=i==3,ylabel=j==0)
            ax.set_title((VL[v]+'\n' if i==0 else '')+f'Valid hour {h:02d} UTC',pad=3)
            ax.tick_params(labelbottom=True)
    save(fig,'fig_diurnal_core')

def dm_map():
    d = pd.read_csv(ROOT/'data/interim/fig4_plot_snapshot.csv',dtype={'station_id':str})
    prov = json.loads((ROOT/'data/interim/fig4_plot_snapshot_provenance.json').read_text())
    assert prov['snapshot_sha256']==hashlib.sha256((ROOT/'data/interim/fig4_plot_snapshot.csv').read_bytes()).hexdigest()
    fig,axes = plt.subplots(3,3,figsize=(W,6.1))
    fig.subplots_adjust(left=.025,right=.98,bottom=.30,top=.93,wspace=.06,hspace=1.05)
    for p,ax in enumerate(axes.flat):
        s = d[d.panel==p]; coast(ax)
        n = s[s.plot_sign==0]
        ax.scatter(n.lon,n.lat,s=5,facecolors='none',edgecolors='0.55',linewidths=.35,zorder=2)
        for sign,marker in [(1,'^'),(-1,'v')]:
            n=s[s.plot_sign==sign]
            ax.scatter(n.lon,n.lat,s=12,c=n[['fill_r','fill_g','fill_b']].values,
                marker=marker,edgecolors='0.05',linewidths=.35,zorder=3)
        v=s['var'].iloc[0]; L=s.lead.iloc[0]; day=L/24 if v=='t2' else L
        ax.set_title(f'({chr(97+p)}) {VL[v]}\nDay {day:g}',loc='left',pad=3)
        ax.text(.02,-.17,f'Better {(s.plot_sign==1).mean()*100:.0f}%\nWorse {(s.plot_sign==-1).mean()*100:.0f}%',
            transform=ax.transAxes,va='top',ha='left')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(.5)
    handles = [Line2D([],[],ls='none',marker=m,color='0.15',markerfacecolor=fc,
        markersize=5,label=label) for m,fc,label in [('^','0.5','Significant improvement'),
        ('v','0.5','Significant deterioration'),('o','white','Not significant')]]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.10,.13),ncol=1,
        frameon=False,labelspacing=.15)
    cax=fig.add_axes([.20,.055,.62,.020])
    cb=fig.colorbar(ScalarMappable(norm=Normalize(-30,30),cmap='RdBu'),cax=cax,orientation='horizontal')
    cb.set_ticks([-30,0,30],labels=['-30','0','+30'])
    cax.text(.5,2.0,'GraphCast PCRPS improvement over GFS (%)',transform=cax.transAxes,ha='center')
    save(fig,'fig_dm_map_core_graphcast')

def zones():
    st=pd.read_csv(ROOT/'data/interim/stations_analysis.csv')
    fig,axes=plt.subplots(3,2,figsize=(W,6.2))
    fig.subplots_adjust(left=.27,right=.985,bottom=.075,top=.945,wspace=.90,hspace=.67)
    cmap=plt.get_cmap('RdBu')(np.linspace(0,1,256))
    cmap[:,:3]=.50+.50*cmap[:,:3]
    for k,(ax,v) in enumerate(zip(axes.flat,V)):
        models=MODELS if v!='pr24' else MODELS[:1]
        b=data('core',v); lead=120 if v in ('t2','ws') else 5
        mat=np.full((6,len(models)),np.nan); sig=np.zeros_like(mat,dtype=bool)
        for i,z in enumerate(ZONES+['All']):
            for j,m in enumerate(models):
                s=select(b,m,'pcrps',zone=z); s=s[s.lead==lead]
                if len(s): mat[i,j]=100*s.skill.iloc[0]; sig[i,j]=s.ci_lo.iloc[0]>0 or s.ci_hi.iloc[0]<0
        ax.imshow(mat,cmap=ListedColormap(cmap),vmin=-20,vmax=20,aspect='auto')
        for (i,j),val in np.ndenumerate(mat):
            if np.isfinite(val):
                label=f'{val:+.0f}' if val!=0 else '0'
                ax.text(j,i,label+('*' if sig[i,j] else ''),ha='center',va='center',color='black')
        use=st if v in ('t2','ws') else st[st[f'use_{v}']]
        counts=use.zone.value_counts()
        labels=([f'{z} ({counts.get(z,0)})' for z in ZONES]+[f'All ({len(use)})'] if k%2==0
            else [f'n={counts.get(z,0)}' for z in ZONES]+[f'n={len(use)}'])
        ax.set_yticks(range(6),labels=labels)
        ax.set_xticks(range(len(models)),labels=['GC','PW','FCN2'][:len(models)])
        ax.set_title(VL[v]+'\nDay 5',pad=4,y=1.0)
        ax.tick_params(length=0)
        for sp in ax.spines.values(): sp.set_visible(False)
    axes[2,1].axis('off')
    axes[2,1].text(-.6,.98,'PCRPS skill vs GFS (%)\n\nGC: GraphCast\nPW: Pangu-Weather\nFCN2: FourCastNet v2-small\n\n+ / - : skill sign\n* : 95% interval excludes 0\nCounts: eligible stations',
        transform=axes[2,1].transAxes,va='top',linespacing=1.35)
    save(fig,'fig_zone_core_pcrps')

def paired(kind):
    fig,axes=plt.subplots(3,2,figsize=(W,5.7))
    fig.subplots_adjust(left=.135,right=.98,bottom=.10,top=.94,wspace=.33,hspace=.64)
    for ax,v in zip(axes.flat,V):
        b=data('core' if kind=='oos' else 'ic',v)
        for m in MODELS:
            if kind=='oos':
                for sc,ls in [('pcrps','-'),('oos_idr','--'),('oos_emos',':')]:
                    draw(ax,select(b,m,sc),v,m,ls=ls,band=False)
            else:
                for suffix,ls,op in [('', '-',False),('_ifs','--',True)]:
                    draw(ax,select(b,m+suffix,'pcrps'),v,m,ls=ls,open_marker=op)
        decorate(ax,xlabel=True,ylabel=True)
        ax.set_title(VL[v],pad=4,y=1.0)
    ax=axes[2,1]; ax.axis('off')
    h=model_handles(styles=False)
    h += [Line2D([],[],color='black',ls=ls,lw=1.2,label=label) for ls,label in
        ([('-','PCRPS (in-sample IDR)'),('--','IDR out-of-sample'),(':','EMOS out-of-sample')]
        if kind=='oos' else [('-','GFS initial conditions'),('--','IFS initial conditions')])]
    if kind=='ic':
        h[-2].set_marker('o'); h[-2].set_markerfacecolor('black')
        h[-1].set_marker('o'); h[-1].set_markerfacecolor('white')
        h += [Patch(facecolor='0.9',edgecolor='none',label='95% intervals')]
    ax.legend(handles=h,loc='upper left',bbox_to_anchor=(-.06,.98),frameon=False,
        labelspacing=.55,borderpad=0,handlelength=2.2,handletextpad=.5)
    save(fig,'fig_insample_vs_oos_core' if kind=='oos' else 'fig_ic_decomposition_pcrps')

if __name__=='__main__':
    for kind in ('core','ic'):
        for v in V: require_current(TB/f'{kind}_{v}_meta.json')
    stations(); skill(); diurnal(); dm_map(); zones(); paired('oos'); paired('ic')
    check_path = ROOT / 'paper/submission/apjas/Figure_source_checks.json'
    check_path.parent.mkdir(parents=True, exist_ok=True)
    check_path.write_text(
        json.dumps({'source_font_minimum_pt':FS,'designed_column_width_inches':W,
        'numerical_tables_modified':False,'figures':checks},indent=2)+'\n')
