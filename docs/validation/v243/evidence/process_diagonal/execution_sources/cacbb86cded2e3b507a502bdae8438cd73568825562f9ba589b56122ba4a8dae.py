"""模块职责：只核对失败母单元的三种中心八面体对角线质量上界，不改网格。"""
from pathlib import Path
import sys
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import runSU2RocketPhaseV222 as guard
from reviewSU2ExplicitQuality import quality
BASE=ROOT/'outputs/v243_validation';DOC=ROOT/'docs/validation/v243'
CLOSED=False


def review():
    z=np.load(str(BASE/'refine_audit01/mesh.npz'));x=z['vertices'];t=z['tetrahedra'];support=z['nodeSourceSupport'];parent=z['parentCells']
    bad=np.unique(parent[z['SICN']<.1]);records=[]
    for father in bad:
        group=np.flatnonzero(parent==father);nodes=np.unique(t[group]);mid=[int(i) for i in nodes if support[i,0]!=support[i,1]];assert len(mid)==6
        pairs=[(a,b) for j,a in enumerate(mid) for b in mid[:j] if not set(support[a])&set(support[b])];assert len(pairs)==3
        variants=[]
        corners=t[group][np.any(support[t[group],0]==support[t[group],1],axis=1)];assert len(corners)==4
        for a,b in pairs:
            rest=sorted(set(mid)-{a,b});cycle=[rest[0]]
            while len(cycle)<4:
                cycle.append(next(v for v in rest if v not in cycle and set(support[v])&set(support[cycle[-1]])))
            assert set(support[cycle[-1]])&set(support[cycle[0]])
            cells=np.array([[a,b,cycle[j],cycle[(j+1)%4]] for j in range(4)],dtype=np.int64)
            delta=x[cells[:,1:]]-x[cells[:,0,None]];negative=np.linalg.det(delta)<0
            cells[negative]=cells[negative][:,[1,0,2,3]]
            q=quality(x[np.r_[corners,cells]])
            variants.append(dict(diagonalNodes=[a,b],minimumSICN=float(min(q)),cells=cells.tolist()))
        records.append(dict(parent=int(father),originalMinimum=float(min(z['SICN'][group])),bestPossibleMinimum=max(r['minimumSICN'] for r in variants),variants=variants))
    save(BASE/'diagonal01/report.json',dict(passed=True,badParents=records,canMeetOriginalQuality=all(r['bestPossibleMinimum']>=.1 for r in records),meshChanged=False,physicalQualification=False))
    print(records)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 diagonal review closed; use NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':review()
    else:
        guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
        guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',DOC/'PLAN.md',DOC/'METHOD_DECISION.md',BASE/'refine_audit01/mesh.npz',BASE/'refine_audit01/report.json']
        raise SystemExit(0 if guard.supervise('diagonal',[sys.executable,str(Path(__file__).resolve()),'task'],20,256,locks,BASE,'offline') else 1)
