"""模块职责：用已审核保存结果生成共同色标科研图；不重算CFD、热流或原面映射。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize,SymLogNorm,LogNorm
import runMUSCLRampV239 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC


def main():
    out=BASE/'visual01';out.mkdir(exist_ok=False);m=read(BASE/'workbench01/rocket_restart_comparison.json')
    assert read(BASE/'workbench01/report.json')['passed'];g=np.load(str(run.GEOMETRY));updates=m['researchIdentity']['successfulUpdates']
    frames={s:np.loadtxt(str(BASE/'real_audit01'/('original_faces_%06d.csv'%s)),delimiter=',',skiprows=1) for s in updates};last=max(updates)
    tri=g['originalVertices'][g['originalFaces']];right=np.array([1.,0.,.10]);right/=np.linalg.norm(right)
    up=np.array([-.025,1.,.25]);up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    projected=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(4,2,figsize=(14,9.5),facecolor='white',gridspec_kw=dict(height_ratios=[1,1,1,1.4]))
    fig.subplots_adjust(left=.07,right=.87,top=.89,bottom=.11,hspace=.4,wspace=.2)
    specs=(('pressure',2,1e-3,'Wall pressure [kPa], log scale','viridis'),
        ('heatFlux',3,1e-6,'Heat flux into wall [MW/m²]','coolwarm'),('wallDeviation',None,1,'Wall fluid T minus given Tw [K]','coolwarm'))
    for row,(key,col,scale,label,cmap) in enumerate(specs):
        lo,hi=np.array(m['comparisonColorRanges'][key+'_0'])*scale
        if key=='pressure':assert lo>0;norm=LogNorm(lo,hi)
        elif key=='heatFlux':limit=max(abs(lo),abs(hi));norm=SymLogNorm(linthresh=.01,vmin=-limit,vmax=limit,base=10)
        else:norm=Normalize(lo,hi)
        for column,step in enumerate((0,last)):
            ax=axs[row,column];a=frames[step];values=(a[:,col] if col is not None else a[:,4]-a[:,5])*scale
            ax.add_collection(PolyCollection(projected[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(values[order]))))
            ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.margins(.025);ax.set_axis_off()
            ax.set_title('%s\nBranch update %d'%(label,step),loc='left',fontsize=9,pad=1)
        box=axs[row,1].get_position();cax=fig.add_axes([.9,box.y0+.1*box.height,.012,.8*box.height])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    completed=read(BASE/'real_audit01/report.json')['successfulUpdates']
    with (run.CASE/'stages.csv').open() as f:
        stages=[r for r in csv.DictReader(f) if r['stage']=='after_time' and int(r['update'])<=completed]
    for j,name in enumerate(('mass','x momentum','y momentum','z momentum','physical energy','SA')):
        selected=[r for r in stages if int(r['system'])==(0 if j<5 else 2)]
        axs[3,0].plot([int(r['update']) for r in selected],[np.log10(max(float(r['rms%d'%j]),1e-300)) for r in selected],label=name,lw=1)
    axs[3,0].axhline(-8,c='black',ls='--',lw=.7);axs[3,0].set_xlabel('Branch updates');axs[3,0].set_ylabel('Physical component log10 RMS');axs[3,0].legend(fontsize=7,ncol=3);axs[3,0].grid(alpha=.2)
    # 阶次渐进期间改变空间算子，用淡灰底区分，避免把跨阶次变化当完整二阶收敛。
    axs[3,0].axvspan(0,min(101,completed),color='#777777',alpha=.09,zorder=0)
    if completed>=102:axs[3,0].axvline(101.5,color='#777777',ls=':',lw=.8)
    axs[3,1].plot(updates,[g['faceArea']@frames[s][:,3]/1e6 for s in updates],marker='o',ms=3,color='#265b79')
    axs[3,1].set_xlabel('Branch updates');axs[3,1].set_ylabel('Original-surface heat power [MW]');axs[3,1].grid(alpha=.2)
    fig.suptitle('RocketAeroHeat | M5 / alpha 20 deg / given Tw 300 K | original 19,752 faces\nNative MUSCL-ramp HLLC branch: current-state wall diagnostics; NOT physically qualified',fontsize=12,y=.98)
    fig.text(.07,.015,'SU2 8.6 | phase-aware Air | SA-neg WITHFT2 | full MUSCL order from update 102 (startup shaded)\nIteration index is not physical time. Solid wall temperature is unsolved. Signed heat flux retained.',fontsize=9,color='#444444')
    p=out/'rocket_restart_comparison.png';fig.savefig(str(p),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,savedDataRendering=True,sourceModelSHA256=sha(BASE/'workbench01/rocket_restart_comparison.json'),imageSHA256=sha(p),
        pressureUsesExplicitLogScale=True,physicalQualification=False))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runMUSCLRampV239.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
          run.GEOMETRY,run.CASE/'stages.csv',BASE/'workbench01/report.json',BASE/'workbench01/rocket_restart_comparison.json']+list((BASE/'real_audit01').glob('original_faces_*.csv'))
        raise SystemExit(0 if run.run.guard.supervise('render',[sys.executable,str(Path(__file__).resolve()),'task'],30,1024,locks,BASE,'audit') else 1)
