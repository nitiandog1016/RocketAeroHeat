"""模块职责：将已审核原面研究结果画回真实箭体，并显式显示未收敛状态。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,SymLogNorm
from matplotlib.collections import PolyCollection
from runSU2RocketPhaseV222 import ROOT,BASE,CASE,DOC,read,save,sha,supervise

def main():
    report=read(BASE/'audit01/report.json');assert report['passed']
    last=report['records'][-1];step=last['iteration'];folder=BASE/'visual02';folder.mkdir(exist_ok=False)
    source=BASE/'audit01'/('original_faces_%06d.csv'%step)
    a=np.loadtxt(str(source),delimiter=',',skiprows=1);g=np.load(str(BASE/'geometry01/cache.npz'))
    v=g['originalVertices'];tri=v[g['originalFaces']]
    # 正交投影保留长度比；按深度绘制真实三角面，避免3D方形坐标框把细长箭体挤到角落。
    right=np.array([1.,0.,.10]);right/=np.linalg.norm(right)
    up=np.array([-.025,1.,.25]);up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    depth=np.cross(right,up);projected=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@depth)
    fig=plt.figure(figsize=(13,7.5),facecolor='white');gs=fig.add_gridspec(3,2,height_ratios=[.65,.65,1.05])
    for row,(field,scale,label,cmap) in enumerate(((2,1e-3,'Wall pressure [kPa]','viridis'),(3,1e-6,'Heat flux into wall [MW/m2]','coolwarm'))):
        ax=fig.add_subplot(gs[row,:]);values=a[:,field]*scale
        if field==2:norm=Normalize(vmin=float(values.min()),vmax=float(values.max()))
        else:
            limit=float(np.max(abs(values)));norm=SymLogNorm(linthresh=.01,vmin=-limit,vmax=limit,base=10)
        poly=PolyCollection(projected[order],linewidths=0,rasterized=True,facecolors=plt.get_cmap(cmap)(norm(values[order])))
        ax.add_collection(poly);ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.margins(.03);ax.set_axis_off();ax.set_title(label,fontsize=11,loc='left',pad=0)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.018,pad=.01,shrink=.9)
    with (CASE/'history.csv').open() as stream:rows=list(csv.DictReader(stream,skipinitialspace=True))
    h={k.strip():np.array([float(r[k]) for r in rows]) for k in rows[0]};x=h['Inner_Iter']+1
    ax=fig.add_subplot(gs[2,0])
    for name in ('Rho','RhoU','RhoV','RhoW','RhoE','nu'):ax.plot(x,h['rms['+name+']'],label=name,lw=1)
    ax.axhline(-8,color='black',ls='--',lw=.8);ax.set_xlabel('Successful updates');ax.set_ylabel('Native log10 RMS');ax.legend(ncol=3,fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[2,1]);ax.plot(x,h['HF']/1e6,lw=1.2,color='#245e81');ax.set_xlabel('Successful updates');ax.set_ylabel('Native wall heat power [MW]');ax.grid(alpha=.2)
    fig.suptitle('RocketAeroHeat | original 19,752 faces | M5 / alpha 20 deg / prescribed Tw 300 K\nResearch field, update %d: NOT converged or physically qualified'%(step+1),fontsize=12,y=.985)
    fig.text(.06,.02,'SU2 8.6 | phase-aware standard Air | general Roe / implicit Euler | SA-neg WITHFT2\nSigned heat flux retained. Solid wall temperature is unsolved. Geometry aspect ratio is preserved.',fontsize=9,color='#444444')
    fig.subplots_adjust(top=.91,bottom=.13,left=.07,right=.94,hspace=.26,wspace=.28)
    fig.savefig(str(folder/'rocket_research.png'),dpi=170);plt.close(fig)
    save(folder/'report.json',dict(passed=True,iteration=step,sourceSHA256=sha(source),imageSHA256=sha(folder/'rocket_research.png'),physicalQualification=False))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        locks=[Path(__file__).resolve(),BASE/'audit01/report.json',BASE/'geometry01/cache.npz',CASE/'history.csv',DOC/'PLAN.md']
        raise SystemExit(0 if supervise('visual2',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE,'audit') else 1)
