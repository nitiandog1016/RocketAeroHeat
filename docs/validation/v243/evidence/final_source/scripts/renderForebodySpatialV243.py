"""模块职责：用已保存原父面数据画粗细区域对照，明确未收敛与非网格无关身份。"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,Normalize
import runForebodySpatialV243 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;OLD=ROOT/'outputs/v242_validation';CLOSED=True


def render():
    out=BASE/'visual01';out.mkdir(exist_ok=False);a=read(BASE/'flow_audit02/report.json');assert a['passed']
    data=np.loadtxt(str(BASE/'flow_audit02/coarse_fine_original_faces.csv'),delimiter=',',skiprows=1)
    g=np.load(str(ROOT/'outputs/v230_validation/audit02/geometry.npz'));ids=data[:,0].astype(int)-1;tri=g['vertices'][g['faces'][ids]]
    right=np.array([1.,0.,.08]);right/=np.linalg.norm(right);up=np.array([0.,.25,1.]);up-=right*np.dot(right,up);up/=np.linalg.norm(up)
    points=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(3,2,figsize=(12,9),gridspec_kw={'height_ratios':[1,1,1.25]})
    fig.subplots_adjust(left=.085,right=.87,top=.88,bottom=.12,hspace=.45,wspace=.28)
    for row,label,values in [(0,'Parent wall pressure [kPa]',data[:,2:4]/1000),
      (1,'Parent H deviation from freestream [%]',100*(data[:,4:6]/a['freestreamHJkg']-1))]:
        lo=float(values.min());hi=float(values.max())
        if row==1:hi=max(abs(lo),abs(hi));lo=-hi
        norm=LogNorm(lo,hi) if row==0 else Normalize(lo,hi);cmap=plt.get_cmap('viridis' if row==0 else 'coolwarm')
        for col in (0,1):
            ax=axs[row,col];ax.add_collection(PolyCollection(points[order],facecolors=cmap(norm(values[order,col])),edgecolors='none'))
            ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off()
            ax.set_title(('Coarse / update 400' if col==0 else 'Refined / update 800')+' | '+label,loc='left',fontsize=9)
        box=axs[row,1].get_position();cax=fig.add_axes([.90,box.y0+.05*box.height,.015,.9*box.height])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    for folder,label,color,style in [(OLD,'coarse','#7b858c','--'),(BASE,'refined','#216d8b','-')]:
        h=np.loadtxt(str(folder/'run01/physical_history.csv'),delimiter=',',skiprows=1)
        axs[2,0].plot(h[:,0],np.log10(h[:,8]),style,color=color,label=label,lw=1.3)
    axs[2,0].axhline(-8,c='black',ls=':',lw=.8,label='original gate');axs[2,0].legend(fontsize=8)
    axs[2,0].set_xlabel('Steady branch updates');axs[2,0].set_ylabel('Physical energy log10 RMS');axs[2,0].grid(alpha=.2)
    axs[2,1].scatter(tri.mean(axis=1)[:,0],data[:,6],s=3,color='#216d8b',alpha=.25,rasterized=True)
    axs[2,1].axhline(0,c='black',lw=.6);axs[2,1].set_xlabel('Original parent centroid X [m]');axs[2,1].set_ylabel('Refined / coarse pressure change [%]');axs[2,1].grid(alpha=.2)
    fig.suptitle('RocketAeroHeat | same physical forebody / parent-face comparison\nBoth branches unqualified: spatial and iterative effects are NOT separated',fontsize=12,y=.97)
    fig.text(.085,.035,'4434 original parents, 0.321489 m²; refined 17736 children area-aggregated before comparison.\nPressure difference (area-weighted): %.3f%%. No new wall heat or accepted boundary-layer edge.'%(100*a['spatialComparison']['pressureAreaRelativeDifference']),fontsize=9,color='#444')
    path=out/'spatial_comparison.png';fig.savefig(str(path),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,commonPressureAndHColorScales=True,originalParents=4434,notGridIndependence=True,
      imageSHA256=sha(path),sourceComparisonSHA256=sha(BASE/'flow_audit02/coarse_fine_original_faces.csv'),physicalQualification=False))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 spatial render closed; use NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':render()
    else:
        run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC;size=run.guard.size_bytes
        run.guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        locks=[Path(__file__).resolve(),ROOT/'scripts/runForebodySpatialV243.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'CFD_PLAN.md',
          BASE/'flow_audit02/report.json',BASE/'flow_audit02/coarse_fine_original_faces.csv',ROOT/'outputs/v230_validation/audit02/geometry.npz',BASE/'run01/physical_history.csv',OLD/'run01/physical_history.csv']
        raise SystemExit(0 if run.guard.supervise('render',[sys.executable,str(Path(__file__).resolve()),'task'],30,384,locks,BASE,'audit') else 1)
