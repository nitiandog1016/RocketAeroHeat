"""模块职责：用已保存原面字段画科研图；完整箭体未知面为灰色，前体局部单独放大。"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,Normalize
import runForebodyEulerV240 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC


def render():
    out=BASE/'visual02';out.mkdir(exist_ok=False)
    m=read(BASE/'workbench01/forebody_research.json');a=read(BASE/'audit02/report.json');assert read(BASE/'workbench01/report.json')['passed']
    xyz=np.asarray(m['vertices']);tri=xyz[np.asarray(m['faces'])]
    right=np.array([1.,0.,.08]);right/=np.linalg.norm(right);up=np.array([0.,.25,1.]);up-=right*np.dot(right,up);up/=np.linalg.norm(up)
    points=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    selected=np.zeros(19752,bool);selected[np.asarray(m['researchIdentity']['originalFaceIDs'])-1]=True
    fig,axs=plt.subplots(3,2,figsize=(13,9),gridspec_kw={'height_ratios':[1,1,1.25]},facecolor='white')
    fig.subplots_adjust(left=.075,right=.87,top=.88,bottom=.115,hspace=.5,wspace=.25)
    for row,key,label in [(0,'pressure','Wall pressure [kPa]'),(1,'Hdeviation','Total H deviation from freestream [%]')]:
        f=next(f for f in m['fields'] if f['key']==key+'_200')
        values=np.array([float('nan') if v is None else v for v in f['values']]);scale=.001 if key=='pressure' else 1.
        lo,hi=np.asarray(m['comparisonColorRanges'][key+'_200'])*scale
        # 此图只比较同一末帧的全模型/放大图；总焓用末帧对称色标，避免早期瞬态掩盖局部误差。
        if key=='Hdeviation':hi=float(np.max(abs(values[selected])));lo=-hi
        norm=LogNorm(lo,hi) if key=='pressure' else Normalize(lo,hi);cmap=plt.get_cmap('viridis' if key=='pressure' else 'coolwarm')
        colors=np.tile([.78,.80,.82,1.],(len(tri),1));colors[selected]=cmap(norm(values[selected]*scale))
        for col in (0,1):
            ax=axs[row,col];which=order if col==0 else order[selected[order]]
            ax.add_collection(PolyCollection(points[which],facecolors=colors[which],linewidths=0));ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.margins(.03);ax.set_axis_off()
            ax.set_title(label+' | '+('full original surface' if col==0 else 'forebody detail'),loc='left',fontsize=9,pad=2)
        box=axs[row,1].get_position();cax=fig.add_axes([.90,box.y0+.05*box.height,.014,.9*box.height])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    h=np.loadtxt(str(run.CASE/'physical_history.csv'),delimiter=',',skiprows=1)
    for j,name in enumerate(('mass','x momentum','y momentum','z momentum','physical energy')):
        axs[2,0].plot(h[:,0],np.log10(h[:,4+j]),label=name,lw=1)
    axs[2,0].axhline(-8,c='black',ls='--',lw=.8);axs[2,0].set_ylabel('Physical log10 RMS');axs[2,0].set_xlabel('Steady branch updates');axs[2,0].grid(alpha=.2);axs[2,0].legend(fontsize=7,ncol=2)
    oh=np.loadtxt(str(run.CASE/'outlet_history.csv'),delimiter=',',skiprows=1)
    axs[2,1].plot(oh[:,0],oh[:,2],label='minimum native outlet normal Mach',color='#275979');axs[2,1].axhline(1,c='black',ls='--',lw=.8,label='incoming-characteristic threshold')
    axs[2,1].set_xlabel('Steady branch updates');axs[2,1].set_ylabel('Normal Mach');axs[2,1].grid(alpha=.2);axs[2,1].legend(fontsize=7)
    fig.suptitle('RocketAeroHeat | M5 / alpha 20 deg | planar-outlet forebody / update 200\nResearch only: pressure stable; residuals and net balances NOT qualified',fontsize=12,y=.97)
    fig.text(.075,.035,'4434 / 19752 original faces | regional area 0.321489 m² | grey = outside this regional calculation\nCorrected global-ID face mapping. No viscous heat flux or solid wall temperature is solved. Iteration is not physical time.',fontsize=9,color='#444444')
    target=out/'forebody_research.png';fig.savefig(str(target),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,savedFieldsOnly=True,completeOriginalSurfaceShown=True,unknownGray=True,
      imageSHA256=sha(target),sourcePayloadSHA256=sha(BASE/'workbench01/forebody_research.json'),totalHUsesCurrentFrameSymmetricRange=True,physicalQualification=False))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':render()
    else:
        run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC
        size=run.guard.size_bytes;run.guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        inputs=[Path(__file__).resolve(),ROOT/'scripts/runForebodyEulerV240.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
          BASE/'workbench01/forebody_research.json',BASE/'workbench01/report.json',BASE/'audit02/report.json',run.CASE/'physical_history.csv',run.CASE/'outlet_history.csv']
        raise SystemExit(0 if run.guard.supervise('render2',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,inputs,BASE,'audit') else 1)
