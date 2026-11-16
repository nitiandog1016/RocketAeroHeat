"""模块职责：展示控制壳壁法向抵消和同阶段热壁量，禁止冒充真实箭体云图。"""
from pathlib import Path
import csv,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import runFASCompatibilityV245 as run
from runSU2RocketPhaseV222 import ROOT,read,save
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def render():
    out=BASE/'visual01';out.mkdir(exist_ok=False);r=read(BASE/'audit01/report.json')
    wall=np.loadtxt(str(BASE/'observe01/walls.csv'),delimiter=',',skiprows=1);f=wall[wall[:,0]==0];c=wall[wall[:,0]==1]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(12,5));a=fig.add_subplot(121,projection='3d');p=f[:,6:9]
    a.scatter(*p.T,color='#30668c',s=42,label='6 fine wall nodes')
    for i in range(6):
        for j in range(i+1,6):
            if abs(p[i]@p[j])<.5:a.plot(*np.array([p[i],p[j]]).T,color='#a8b3bd',linewidth=.8)
    a.quiver(*p.T,*f[:,3:6].T,length=.45,color='#30668c')
    a.scatter(*c[:,6:9].T,color='#a84335',marker='x',s=80,label='1 coarse wall node')
    a.set(xlabel='x [m]',ylabel='y [m]',zlabel='z [m]',title='Existing diagnostic shell: 24 -> 4 volume nodes')
    a.set_box_aspect((1,1,1));a.legend(loc='lower left',fontsize=8)
    b=fig.add_subplot(122);heat=r['wallHeatStages']
    for l,color,label in ((0,'#30668c','Fine-layer native wall term'),(1,'#a84335','Coarse-layer native wall term')):
        selected=[v for v in heat if v['level']==l]
        b.plot(np.arange(len(selected)),np.abs([v['heatIntoWallW'] for v in selected]),'o-',ms=4,color=color,label=label)
    b.set(yscale='log',xlabel='Recorded wall evaluation (within 5 updates)',ylabel='Absolute wall energy RHS contribution [W]',title='Algebra passes; coarse wall direction is lost')
    b.grid(axis='y',alpha=.2);b.legend(fontsize=8,loc='center right')
    ratio=r['wallCancellation'][0]['ratio']
    fig.suptitle('v245 FAS compatibility control — thermal compatibility REJECTED',fontsize=13)
    fig.text(.5,.03,'Coarse/fine normal-magnitude ratio = %.3e | Diagnostic control only; no rocket heat qualification'%ratio,ha='center',fontsize=9)
    fig.subplots_adjust(left=.04,right=.98,bottom=.17,top=.84,wspace=.3)
    path=out/'fas_wall_cancellation.png';fig.savefig(str(path),dpi=160);plt.close(fig)
    with (out/'wall_heat_stages.csv').open('w',encoding='utf8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['inner','level','heatIntoWallW']);writer.writeheader();writer.writerows(heat)
    save(out/'report.json',dict(passed=True,controlOnly=True,rocketResults=False,scientificImage=str(path),numericRows=len(heat),browserTested=False))


def main():
    if CLOSED:raise SystemExit('v245 rendering closed; follow NEXT_DECISION')
    g=run.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',BASE/'audit01/report.json',BASE/'observe01/walls.csv']
    return g.supervise('render',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v245 rendering closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':render()
    else:raise SystemExit(0 if main() else 1)
