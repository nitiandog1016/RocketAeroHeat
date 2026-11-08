"""将原19752面责任分区画成科研说明图；这不是新的压力或热流结果。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.monotonic();base=ROOT/'outputs/v212_validation'
    r=json.loads((base/'design01/regions.json').read_text())
    d=loadmat(str(ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'),variable_names=['sourceVertices','sourceFaces'])
    v=d['sourceVertices'];f=d['sourceFaces'].astype(int)-1
    colors={'nose':'#236797','contraction':'#167d70','afterbody':'#be7838','remainder':'#bec4cc'}
    facecolors=np.tile(matplotlib.colors.to_rgba(colors['remainder']),(len(f),1))
    for name,reg in r['regions'].items():facecolors[np.array(reg['ownedFaceIDs'])-1]=matplotlib.colors.to_rgba(colors[name])
    fig=plt.figure(figsize=(12,6.7));fig.patch.set_facecolor('white')
    # 模块一：完整原表面、单一责任着色；重叠只在表格中列出，避免面积双算。
    ax=fig.add_axes([.02,.20,.67,.70],projection='3d')
    ax.add_collection3d(Poly3DCollection(v[f],facecolors=facecolors,edgecolors='none'))
    ax.set_xlim(v[:,0].min(),v[:,0].max());ax.set_ylim(v[:,1].min(),v[:,1].max());ax.set_zlim(v[:,2].min(),v[:,2].max())
    if hasattr(ax,'set_box_aspect'):ax.set_box_aspect(np.ptp(v,axis=0))
    ax.view_init(elev=18,azim=-68);ax.set_xlabel('')
    # 细长箭体横向轴刻度易重叠；几何比例保持，只显示轴向定量刻度。
    ax.set_yticks([]);ax.set_zticks([])
    fig.text(.37,.34,'X [m]',size=11,color='#364255')
    ax.grid(False)
    for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.pane.fill=False
    fig.text(.045,.945,'Sequential overlapping regions on the original rocket',size=17,weight='bold',color='#243347')
    fig.text(.045,.908,'Geometry / ownership design only. No new pressure or heat-flux result.',size=11,color='#545e6d')
    ax.legend(handles=[Patch(facecolor=colors[k],label=t) for k,t in [('nose','Nose'),('contraction','Contraction'),
        ('afterbody','Short afterbody'),('remainder','Remaining full-body responsibility')]],loc='upper left',fontsize=8,frameon=False)
    # 模块二：必要载荷与待解面积分开显示，不把内存下界写成实际峰值。
    fig.text(.72,.80,'Required layer payload',size=12,weight='bold',color='#243347')
    y=.74
    for k,label in [('nose','Nose'),('contraction','Contraction'),('afterbody','Short afterbody')]:
        reg=r['regions'][k]
        fig.text(.72,y,label+'   %.2f GiB'%reg['layerPayloadLowerBoundGiB'],color=colors[k],size=11,weight='bold')
        fig.text(.72,y-.035,'%d solve faces / %d owned'%(reg['faceCount'],reg['ownedFaceCount']),size=9,color='#545e6d');y-=.10
    fig.text(.72,.40,'Exterior mesh and runtime peak\nare not included. Fit is unproven.',size=10,color='#7d4a2c',linespacing=1.5)
    fig.text(.72,.29,'Overlap: 736 / 1,104 original faces\n336 missing-heat roots inside core\n32 layers / first height 2 micrometres',size=9,color='#364255',linespacing=1.6)
    fig.text(.045,.115,'All 19,752 original faces retained  |  Total area 1.9646292566 m2  |  New qualified heat area: 0',size=10,color='#364255')
    fig.text(.045,.07,'Internal coupling, terminal boundary, convergence and heat validation remain open.',size=10,color='#7d4a2c')
    fig.savefig(str(base/'region_design.png'),dpi=150,facecolor='white');plt.close(fig)
    print(json.dumps(dict(seconds=time.monotonic()-start,faces=len(f),figure='outputs/v212_validation/region_design.png')))

if __name__=='__main__':main()
