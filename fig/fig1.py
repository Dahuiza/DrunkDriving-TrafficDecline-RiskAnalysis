"""
FIG1: Spatial distribution of drunk driving indicators under the pre-pandemic baseline
and national monthly trends across three scenarios.
Author: Hui Liu
Date: 2026-04-15
Description: This script processes Excel and text-based drunk driving datasets,
             performs Albers projection on geographic shapefiles, and generates
             a publication-quality composite figure (Maps + Time-series)
             highlighting trends across three distinct time scenarios.
"""

import struct, os, math
import numpy as np
import pandas as pd
import matplotlib
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon, Rectangle
import warnings, time
from config import __version, __date

start_ = time.time()
warnings.filterwarnings('ignore')
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# ═══════════════════════════════════════════════════════════════════════
# LAYOUT
# ═══════════════════════════════════════════════════════════════════════
FIG_WIDTH  = 7.0
FIG_HEIGHT = 6.1

MAP_COL_LEFT  = [0.01, 0.345, 0.68]
MAP_COL_WIDTH = 0.315
MAP_ROW_HEIGHT = 2.5 / FIG_HEIGHT
TOP_MARGIN     = 0.025
MAP_ROW_BOTTOM_I = 1 - TOP_MARGIN - MAP_ROW_HEIGHT

CBAR_LEFT   = 0.08
CBAR_BOTTOM = -0.05
CBAR_WIDTH  = 0.58
CBAR_HEIGHT = 0.045

CBAR_FIG_BOTTOM = MAP_ROW_BOTTOM_I + CBAR_BOTTOM * MAP_ROW_HEIGHT
TL_LEFT   = 0.07
TL_WIDTH  = 0.86
TL_HEIGHT = 2.2 / FIG_HEIGHT - 0.02
TL_TOP    = CBAR_FIG_BOTTOM - 0.025
TL_BOTTOM = TL_TOP - TL_HEIGHT - 0.04

LEGEND_Y  = TL_BOTTOM - 0.15
# ═══════════════════════════════════════════════════════════════════════


FIG1_PATH = r'../data/fig1'
EXCEL_FP = os.path.join(FIG1_PATH, "fig1_data_without35_th15.xlsx")
TXT_FP   = os.path.join(FIG1_PATH, 'pop_1year_all_fig0_1211.txt')
CITY_SHP = os.path.join(FIG1_PATH, "shp2020", 'City.shp')
CITY_DBF = os.path.join(FIG1_PATH, "shp2020", 'City.dbf')
PROV_SHP = os.path.join(FIG1_PATH, "shp2020", 'Province.shp')
NINE_SHP = os.path.join(FIG1_PATH, "shp2020", 'Nineline.shp')
NINE_DBF = os.path.join(FIG1_PATH, "shp2020", 'Nineline.dbf')
SAVE_DIR = os.path.join(r'../output', f"{__date}-{__version}")
save_fig_fn = f'fig1.jpg'
save_fig_pdf = f'fig1.pdf'
os.makedirs(SAVE_DIR, exist_ok=True)


# ── Albers Projection ──────────────────────────────────────────────────
def albers(lon, lat, lon0=105, lat0=0, lat1=25, lat2=47):
    lon_r=math.radians(lon); lat_r=math.radians(lat)
    lon0_r=math.radians(lon0)
    lat1_r=math.radians(lat1); lat2_r=math.radians(lat2)
    n=(math.sin(lat1_r)+math.sin(lat2_r))/2
    C=math.cos(lat1_r)**2+2*n*math.sin(lat1_r)
    rho0=math.sqrt(C-2*n*math.sin(math.radians(lat0)))/n
    theta=n*(lon_r-lon0_r)
    rho=math.sqrt(max(0,C-2*n*math.sin(lat_r)))/n
    return rho*math.sin(theta), rho0-rho*math.cos(theta)

def project_ring(ring):
    out=[]
    for x,y in ring:
        try: out.append(albers(x,y))
        except: pass
    return out

def read_dbf_records(dbf_path):
    records=[]
    with open(dbf_path,'rb') as f:
        hdr=f.read(32)
        num_records=struct.unpack('<I',hdr[4:8])[0]
        header_size=struct.unpack('<H',hdr[8:10])[0]
        record_size=struct.unpack('<H',hdr[10:12])[0]
        fields=[]; pos=32
        while pos<header_size-1:
            fd=f.read(32)
            if not fd or fd[0]==0x0D: break
            name=fd[0:11].replace(b'\x00',b'').decode('gbk',errors='replace')
            flen=fd[16]; fields.append((name,flen)); pos+=32
        f.seek(header_size)
        for _ in range(num_records):
            raw=f.read(record_size)
            if not raw: break
            rec={}; offset=1
            for name,flen in fields:
                rec[name]=raw[offset:offset+flen].decode('gbk',errors='replace').strip()
                offset+=flen
            records.append(rec)
    return records

def read_shp(shp_path, dbf_path=None):
    records=read_dbf_records(dbf_path) if dbf_path else []
    shapes=[]
    with open(shp_path,'rb') as f:
        f.read(100); idx=0
        while True:
            rec_hdr=f.read(8)
            if len(rec_hdr)<8: break
            content_len=struct.unpack('>i',rec_hdr[4:8])[0]*2
            content=f.read(content_len)
            if len(content)<4: break
            shp_type=struct.unpack('<i',content[0:4])[0]
            attr=records[idx] if idx<len(records) else {}
            if shp_type in (3,5):
                nparts=struct.unpack('<i',content[36:40])[0]
                npoints=struct.unpack('<i',content[40:44])[0]
                parts=[
                    struct.unpack('<i',
                                  content[44+i*4:48+i*4])[0] for i in range(nparts)
                ]
                ps=44+nparts*4
                points=[(struct.unpack('<d',content[ps+i*16:ps+i*16+8])[0],
                         struct.unpack('<d',content[ps+i*16+8:ps+i*16+16])[0])
                        for i in range(npoints)]
                rings=[]
                for pi,p0 in enumerate(parts):
                    p1=parts[pi+1] if pi+1<len(parts) else len(points)
                    rings.append(project_ring(points[p0:p1]))
                shapes.append((attr,rings))
            else:
                shapes.append((attr,[]))
            idx+=1
    return shapes

print("Loading shapefiles...")
city_shapes=read_shp(CITY_SHP,CITY_DBF)
prov_shapes=read_shp(PROV_SHP)
nine_shapes=read_shp(NINE_SHP,NINE_DBF)

all_px,all_py=[],[]
for _,rings in city_shapes:
    for ring in rings:
        all_px+=[p[0] for p in ring]; all_py+=[p[1] for p in ring]
XMIN,XMAX=min(all_px),max(all_px)
YMIN,YMAX=min(all_py),max(all_py)

# ── Data Loading & Processing ──────────────────────────────────────────
def assign_sc(ym):
    if '2016-10'<=ym<='2019-09': return 'I'
    if '2020-01'<=ym<='2020-03': return 'II'
    if '2020-04'<=ym<='2020-09': return 'III'
    return 'other'

print("Loading data...")
df_city = pd.read_excel(EXCEL_FP)
df_city['scenario'] = df_city['year_and_month'].apply(assign_sc)
df_city = df_city[df_city['drink_num'] > 0]

sc1_all = df_city[df_city['scenario']=='I'].groupby('city').agg(
    drunk=('drink_num','mean'), acc=('acc_num','mean'), rate=('acc_ratio','mean')
).reset_index()
city_vals_abc = {str(int(r['city'])): {'drunk':r['drunk'],'acc':r['acc'],'rate':r['rate']}
                 for _,r in sc1_all.iterrows()}

rows=[]
with open(TXT_FP,'r') as f:
    for line in f:
        p=line.strip().split(',')
        if len(p)<4: continue
        rows.append({'ym':p[0],'drunk':int(p[1]),'acc':int(p[2]),'rate':float(p[3])})
df_time=pd.DataFrame(rows)
df_time['date']     = pd.to_datetime(df_time['ym'],format='%Y-%m')
df_time['scenario'] = df_time['ym'].apply(assign_sc)
df_time=(df_time[df_time['scenario']!='other'].sort_values('date').reset_index(drop=True))

# ── North arrow ───────────────────────────────────────────────────────
def add_north_arrow_data(ax,anchor_x,anchor_y,dc_w,dc_h):
    hw=dc_w*0.028; h=dc_h*0.10; gap=dc_h*0.012
    cx=anchor_x; cy_tip=anchor_y+h; notch_y=anchor_y+h*0.32
    tip=(cx,cy_tip); bl=(cx-hw,anchor_y); br=(cx+hw,anchor_y); notch=(cx,notch_y)
    ax.add_patch(
        Polygon([tip,bl,notch],
                closed=True,facecolor='grey',edgecolor='black',linewidth=0.8,zorder=25
                ))
    ax.add_patch(
        Polygon([tip,br,notch],
                closed=True,facecolor='white',edgecolor='black',hatch='////',linewidth=0.8,zorder=25
                ))
    ax.add_patch(
        Polygon([tip,bl,notch,br],
                closed=True,facecolor='none',edgecolor='black',linewidth=1.0,zorder=26
                ))
    ax.text(cx,cy_tip+gap,'N',ha='center',va='bottom',fontsize=7,fontweight='bold',zorder=27)

def add_scale_bar_data(ax,xmin,xmax,ymin,ymax,x_left_frac):
    scale_proj=(500*1.60934)/6371.0
    W=xmax-xmin; H=ymax-ymin
    bar_x0=xmin+W*x_left_frac; bar_y0=ymin+H*0.04
    bar_x1=bar_x0+scale_proj; bar_xm=(bar_x0+bar_x1)/2; bar_h=H*0.011
    ax.add_patch(
        Rectangle((bar_x0,bar_y0),
                  scale_proj/2,bar_h,facecolor='black',edgecolor='black',linewidth=0.5,zorder=10
                  ))
    ax.add_patch(
        Rectangle((bar_xm,bar_y0),
                  scale_proj/2,bar_h,facecolor='white',edgecolor='black',linewidth=0.5,zorder=10
                  ))
    label_y=bar_y0+bar_h+H*0.012
    for val,bx in [(0,bar_x0),(250,bar_xm),(500,bar_x1)]:
        ax.text(bx,label_y,str(val),ha='center',va='bottom',fontsize=6,zorder=11)
    ax.text(bar_xm,bar_y0-H*0.013,'mi',ha='center',va='top',fontsize=6,zorder=11)
    return bar_x0,bar_x1,label_y+H*0.055,bar_xm

def draw_map(ax,value_key,cmap_name,panel_letter,city_data,norm,show_cbar=False):
    cmap=plt.cm.get_cmap(cmap_name)
    vals=[city_data.get(s[0].get('市代码','').strip(),{}).get(value_key,np.nan) for s in city_shapes]
    for (attr,rings),val in zip(city_shapes,vals):
        color=cmap(norm(val)) if not np.isnan(val) else '#d0d0d0'
        for ring in rings:
            if len(ring)<3: continue
            xs=[p[0] for p in ring]; ys=[p[1] for p in ring]
            ax.fill(xs,ys,color=color,linewidth=0,zorder=1)
            ax.plot(xs+[xs[0]],ys+[ys[0]],color='#999999',linewidth=0.07,zorder=2)
    for _,rings in prov_shapes:
        for ring in rings:
            if len(ring)<2: continue
            xs=[p[0] for p in ring]; ys=[p[1] for p in ring]
            ax.plot(xs+[xs[0]],ys+[ys[0]],color='#2a2a2a',linewidth=0.55,zorder=3)
    for _,rings in nine_shapes:
        for ring in rings:
            if len(ring)<2: continue
            xs=[p[0] for p in ring]; ys=[p[1] for p in ring]
            ax.plot(xs,ys,color='#111111',linewidth=0.7,zorder=4)
    pad_x=(XMAX-XMIN)*0.02
    x0=XMIN-pad_x; x1=XMAX+pad_x; x_range=x1-x0
    y_range=x_range*(MAP_ROW_HEIGHT*FIG_HEIGHT/(MAP_COL_WIDTH*FIG_WIDTH))
    y_center=(YMIN+YMAX)/2
    ax.set_xlim(x0,x1); ax.set_ylim(y_center-y_range/2,y_center+y_range/2)
    ax.axis('off')
    ax.text(0.02,0.99,panel_letter,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top',ha='left')
    if show_cbar:
        cbar_ax=ax.inset_axes([CBAR_LEFT,CBAR_BOTTOM,CBAR_WIDTH,CBAR_HEIGHT])
        sm=ScalarMappable(cmap=cmap,norm=norm); sm.set_array([])
        cbar=plt.colorbar(sm,cax=cbar_ax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=6,pad=1.5,length=2)
        cbar.outline.set_linewidth(0.4)
        if value_key=='rate':
            cbar.ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        else:
            fmt=ticker.ScalarFormatter(useMathText=True)
            fmt.set_scientific(True); fmt.set_powerlimits((0,0))
            cbar.ax.xaxis.set_major_formatter(fmt)
            cbar.ax.xaxis.offsetText.set_visible(False)
            cbar.ax.figure.canvas.draw()
            offset_str=cbar.ax.xaxis.get_major_formatter().get_offset()
            if offset_str:
                cbar_ax.text(1.04,0.5,offset_str,transform=cbar_ax.transAxes,ha='left',va='center',fontsize=6)
    xlim=ax.get_xlim(); ylim=ax.get_ylim()
    dc_w=xlim[1]-xlim[0]; dc_h=ylim[1]-ylim[0]
    _,_,bar_y_top,bar_xm=add_scale_bar_data(ax,xlim[0],xlim[1],ylim[0],ylim[1],CBAR_LEFT)
    add_north_arrow_data(ax,bar_xm,bar_y_top,dc_w,dc_h)

# ═══════════════════════════════════════════════════════════════════════
# PLOT FIGURE A
# ═══════════════════════════════════════════════════════════════════════
fig=plt.figure(figsize=(FIG_WIDTH,FIG_HEIGHT))
SC_COLORS={'I':'#6BAED6','II':'#E8696B','III':'#74C476'}

abc_norms={}
for vkey,field in [('drunk','drink_num'),('acc','acc_num'),('rate','acc_ratio')]:
    vals=df_city[df_city['scenario']=='I'][field].values
    abc_norms[vkey]=mcolors.Normalize(vmin=np.percentile(vals,2),vmax=np.percentile(vals,98))

ABC_CMAPS={'drunk':'YlOrRd','acc':'BuGn','rate':'PuRd'}

for col,(vkey,letter) in enumerate(zip(['drunk','acc','rate'],['a','b','c'])):
    ax=fig.add_axes([MAP_COL_LEFT[col], MAP_ROW_BOTTOM_I, MAP_COL_WIDTH, MAP_ROW_HEIGHT])
    draw_map(ax, vkey, ABC_CMAPS[vkey], letter, city_vals_abc, abc_norms[vkey], show_cbar=True)

# ── Time-series Plot ───────────────────────────────────────────────────
ax_bar =fig.add_axes([TL_LEFT, TL_BOTTOM, TL_WIDTH, TL_HEIGHT])
ax_line=ax_bar.twinx()
sc_arr=df_time['scenario'].values
x_idx=np.arange(len(df_time))
x_lbls=[d.strftime('%Y.%m') for d in df_time['date']]
bar_w=0.35
ax_bar.bar(x_idx-bar_w/2,df_time['drunk'],width=bar_w,color=[SC_COLORS[s] for s in sc_arr],alpha=0.90,zorder=2)
ax_bar.bar(x_idx+bar_w/2,df_time['acc'],
           width=bar_w,color=[SC_COLORS[s] for s in sc_arr],
           alpha=0.45,hatch='///',edgecolor='white',linewidth=0.15,zorder=2)
ax_line.plot(x_idx,df_time['rate'],color='#C0392B',marker='o',markersize=2.8,linewidth=1.4,zorder=5)
ax_line.set_ylabel('Crash Incidence Rate',fontsize=8,color='black')
ax_line.tick_params(axis='y',labelsize=7,labelcolor='black')
ax_line.set_ylim(0,0.75)
ax_line.spines['top'].set_visible(False)
ax_line.spines['right'].set_color('black')

cur_sc=sc_arr[0]; seg_start=0
for j in range(1,len(sc_arr)):
    if sc_arr[j]!=cur_sc:
        ax_bar.axvspan(seg_start-.5,j-.5,color=SC_COLORS[cur_sc],alpha=0.10,zorder=0)
        ax_bar.axvline(j-.5,color='#E8A015',linestyle='--',linewidth=0.9,alpha=0.8,zorder=3)
        cur_sc=sc_arr[j]; seg_start=j
ax_bar.axvspan(seg_start-.5,len(sc_arr)-.5,color=SC_COLORS[cur_sc],alpha=0.10,zorder=0)

# ── Scenario Labels: Moved to legend, removed text from plot area ──────
ylim_top=df_time['drunk'].max()*1.10; ax_bar.set_ylim(0,ylim_top)

ax_bar.set_xticks(x_idx[::3]); ax_bar.set_xticklabels(x_lbls[::3],rotation=38,ha='right',fontsize=7)
ax_bar.set_ylabel('Count',fontsize=8); ax_bar.spines['top'].set_visible(False)
ax_bar.tick_params(axis='y',labelsize=7)
fmt=ticker.ScalarFormatter(useMathText=True); fmt.set_scientific(True); fmt.set_powerlimits((0,0))
ax_bar.yaxis.set_major_formatter(fmt); ax_bar.yaxis.offsetText.set_fontsize(7)
ax_bar.text(-0.065,1.08,'d',transform=ax_bar.transAxes,fontsize=11,fontweight='bold',va='top')

# ── Legend ─────────────────────────────────────────────────────────────
# Row 1: Data types (Bars + Line)
data_handles = [
    mpatches.Patch(facecolor='#888888', edgecolor='#555555', alpha=0.90,
                   label='Drunk driving cases'),
    mpatches.Patch(facecolor='#aaaaaa', edgecolor='#555555', alpha=0.60,
                   hatch='///', label='Drunk driving related crashes'),
    plt.Line2D([0],[0], color='#C0392B', marker='o', markersize=4,
               label='Crash incidence rate'),
]

sc_handles = [
    mpatches.Patch(facecolor=SC_COLORS['I'],  alpha=0.35, edgecolor=SC_COLORS['I'],
                   label='Scenario I'),
    mpatches.Patch(facecolor=SC_COLORS['II'], alpha=0.35, edgecolor=SC_COLORS['II'],
                   label='Scenario II'),
    mpatches.Patch(facecolor=SC_COLORS['III'],alpha=0.35, edgecolor=SC_COLORS['III'],
                   label='Scenario III'),
]

# Stacked legends: Row 1 directly below the plot, Row 2 below Row 1
LEG1_Y = LEGEND_Y + 0.03     # Row 1 (Data types)
LEG2_Y = LEGEND_Y            # Row 2 (Scenarios)

fig.legend(handles=data_handles, loc='lower center',
           bbox_to_anchor=(0.5, LEG1_Y),
           fontsize=7.5, ncol=3, frameon=False,
           columnspacing=1.2, handlelength=1.8)

fig.legend(handles=sc_handles, loc='lower center',
           bbox_to_anchor=(0.5, LEG2_Y),
           fontsize=7.5, ncol=3, frameon=False,
           columnspacing=1.0, handlelength=1.8)

out=os.path.join(SAVE_DIR, save_fig_fn)
out_pdf=os.path.join(SAVE_DIR, save_fig_pdf)
fig.savefig(out,dpi=500,bbox_inches='tight')
fig.savefig(out_pdf,dpi=150,bbox_inches='tight')
print(f"Time spent: {time.time()-start_:.1f}s  →  {out}")
plt.show()
plt.close(fig)
