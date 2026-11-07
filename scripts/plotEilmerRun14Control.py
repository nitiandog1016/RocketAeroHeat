"""生成可导出的科研对照图；明确未收敛身份，不把基准压力/热流混作求解结果。"""
from pathlib import Path
import argparse
import csv
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(result,case):
    start=time.monotonic();result=Path(result);case=Path(case)
    report=json.loads((result/'report.json').read_text());wall=np.genfromtxt(result/'wall.csv',delimiter=',',names=True)
    stations=list(csv.DictReader((result/'stations.csv').open()))
    hist=np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
    fig,ax=plt.subplots(2,3,figsize=(14,7),layout='constrained')
    blue='#21618c';red='#a93226'
    for axis,quantity,key,scale,label in [(ax[0,0],'p','pressure_Pa',1000,'Pressure (kPa)'),
                                        (ax[0,1],'q','heat_into_wall_W_m2',1000,'Heat flux into wall (kW/m²)')]:
        axis.plot(wall['x_m']/.1017016,wall[key]/scale,color=blue,label='Eilmer research field')
        data=[x for x in stations if x['quantity']==quantity and x['measurement']]
        inside=[x for x in data if x['prediction']];outside=[x for x in data if not x['prediction']]
        axis.scatter([float(x['x_m'])/.1017016 for x in inside],[float(x['measurement'])/scale for x in inside],
                     s=15,facecolors='none',edgecolors='black',label='CUBRC Run 14')
        if outside:axis.scatter([float(x['x_m'])/.1017016 for x in outside],[float(x['measurement'])/scale for x in outside],
                                s=30,marker='x',color=red,label='Outside computed wall')
        axis.axvline(.1017/.1017016,color='.6',linewidth=.7,linestyle='--')
        axis.set(xlabel='Axial x / L (L = 0.1017016 m)',ylabel=label);axis.legend(fontsize=8)
    for key,label in [('global-residual-rel','Global'),('mass-rel','Mass'),('total-energy-rel','Energy')]:
        ax[1,0].semilogy(hist['step'],hist[key],label=label,linewidth=1)
    ax[1,0].axhline(1e-8,color='.5',linestyle=':',label='Original residual target')
    ax[1,0].axvline(50,color='.6',linestyle='--',linewidth=.7)
    ax[1,0].set(xlabel='Native step (second order after step 50)',ylabel='Relative residual');ax[1,0].legend(fontsize=8)
    ax[1,1].plot(wall['x_m']/.1017016,wall['signed_tangent_stress_Pa'],color=blue)
    ax[1,1].axhline(0,color='.5',linewidth=.7);ax[1,1].set(xlabel='Axial x / L',ylabel='Signed tangential viscous traction (Pa)')
    # 全局热峰在尖锐前缘；另列扩张段，避免前缘峰压缩再附附近的可读范围。
    take=wall['block']==1
    ax[0,2].plot(wall['x_m'][take]/.1017016,wall['heat_into_wall_W_m2'][take]/1000,color=blue)
    data=[s for s in stations if s['quantity']=='q' and float(s['x_m'])>=.1017 and s['measurement']]
    ax[0,2].scatter([float(s['x_m'])/.1017016 for s in data],[float(s['measurement'])/1000 for s in data],
                    s=15,facecolors='none',edgecolors='black')
    ax[0,2].set(xlabel='Axial x / L',ylabel='Flare heat flux (kW/m²)')
    ax[0,2].set_title('Flare only — excludes leading-edge peak',fontsize=10)
    ax[1,2].plot(wall['x_m']/.1017016,wall['independent_y_plus'],color=blue)
    ax[1,2].set(xlabel='Axial x / L',ylabel='y+ (normal distance / tangential shear)',title='No mesh-independence qualification')
    for axis in ax.ravel():
        axis.grid(alpha=.15);axis.spines[['top','right']].set_visible(False)
    fig.suptitle('CUBRC Run 14 / Eilmer — research comparison, not physically qualified\n'
                 'Native step %d; ideal N₂, Sutherland transport; Tw = 295.2 K'%report['lastSavedLoadsStep'],fontsize=12)
    for ext in ['png','pdf']:
        dest=result/('comparison.'+ext);assert not dest.exists();fig.savefig(dest,dpi=170)
    plt.close(fig)
    (result/'plot_report.json').write_text(json.dumps(dict(elapsedSeconds=time.monotonic()-start,unconvergedIdentityShown=True))+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--case',required=True)
    a=p.parse_args();plot(a.result,a.case)
