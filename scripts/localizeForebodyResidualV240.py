"""模块职责：仅用已存200场定位新区域残差贡献；不冻结重算、不证明唯一数值原因。"""
from pathlib import Path
import sys
import numpy as np
import runForebodyEulerV240 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh
BASE=run.BASE


def inspect():
    x,_,markers=read_tetra_mesh(run.CASE/'mesh.su2')
    a=np.loadtxt(str(run.CASE/'state_000200.csv'),delimiter=',',skiprows=1);a=a[np.argsort(a[:,0])]
    assert np.array_equal(a[:,1:4],x)
    r=-a[:,10:15].copy();r[:,4]-=300000*r[:,0];square=r*r;total=square.sum(axis=0)
    wall=np.zeros(len(x),bool);wall[np.unique(markers['wall'])]=True
    outer=np.zeros(len(x),bool);outer[np.unique(markers['farfield'])]=True
    outlet=np.zeros(len(x),bool);outlet[np.unique(markers['outlet'])]=True
    groups={'wall_outlet_seam':wall&outlet,'wall_other':wall&~outlet,'outlet_nonwall':outlet&~wall,
      'farfield_only':outer&~wall&~outlet,'interior':~wall&~outer&~outlet}
    assert np.all(sum(mask.astype(int) for mask in groups.values())==1)
    rows=[]
    for name,mask in groups.items():rows.append(dict(group=name,nodes=int(sum(mask)),squaredPhysicalRHSFractions=(square[mask].sum(axis=0)/total).tolist()))
    idx=np.argsort(square[:,4])[-10:][::-1]
    hot=[dict(globalNode=int(i),xyz=x[i].tolist(),group=next(k for k,m in groups.items() if m[i]),physicalEnergyRHS=float(r[i,4]),energySquareFraction=float(square[i,4]/total[4])) for i in idx]
    save(BASE/'diagnostic01/report.json',dict(passed=True,update=200,groups=rows,largestPhysicalEnergyRHSNodes=hot,
      currentStateAndRHSPhaseDifferent=True,uniqueCauseProved=False,flowReevaluations=0,nativeQueries=0,physicalQualification=False,
      sourceSHA256=sha(run.CASE/'state_000200.csv')))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':inspect()
    else:
        run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=run.DOC
        size=run.guard.size_bytes;run.guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        inputs=[Path(__file__).resolve(),ROOT/'scripts/runForebodyEulerV240.py',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/phaseEulerRegionalVolume.py',run.DOC/'PLAN.md',run.CASE/'mesh.su2',run.CASE/'state_000200.csv']
        raise SystemExit(0 if run.guard.supervise('diagnostic',[sys.executable,str(Path(__file__).resolve()),'task'],15,256,inputs,BASE,'audit') else 1)
