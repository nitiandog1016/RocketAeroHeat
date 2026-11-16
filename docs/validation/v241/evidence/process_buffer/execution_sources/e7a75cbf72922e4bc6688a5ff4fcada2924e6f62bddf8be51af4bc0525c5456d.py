"""模块职责：用保存场定位出口缓冲区的固定面约束，输出一次成熟网格采用依据。"""
from pathlib import Path
import sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from phaseEulerRegionalVolume import read_tetra_mesh, topology, extract_region
from reviewSU2ExplicitQuality import quality, bound

ROOT=guard.ROOT; BASE=ROOT/'outputs/v241_validation'; DOC=ROOT/'docs/validation/v241'
OLD=ROOT/'outputs/v240_validation/run01'


def inspect():
    """模块一：只对具名候选和热点求质量，不重新评价CFD算子。"""
    x,t,m=read_tetra_mesh(OLD/'mesh.su2'); topo=topology(x,t,m)
    wallnodes=np.unique(m['wall']); radius=float(np.linalg.norm(x[wallnodes,1:],axis=1).max())
    cut=float(x[np.unique(m['outlet']),0].min())
    near=x[:,0]<=cut+radius
    safe=np.linalg.norm(x[:,1:],axis=1)>2*radius
    mask=np.any(near[t],axis=1)&np.all(safe[t],axis=1)
    reg=extract_region(topo,mask); b=reg['boundary']; outlet=topo['markers']['outlet']
    cap=np.intersect1d(b,outlet); fixed=np.setdiff1d(b,cap)
    q=quality(x[t]); upper=bound(x[topo['faces'][fixed]])[0]
    capupper=bound(x[topo['faces'][cap]])[0]
    diagnostics=guard.read(ROOT/'outputs/v240_validation/diagnostic01/report.json')
    rows=[]
    for row in diagnostics['largestPhysicalEnergyRHSNodes']:
        node=row['globalNode']; adj=np.flatnonzero(np.any(t==node,axis=1)); owned=np.isin(adj,reg['cells'])
        row=dict(row,incidentCells=len(adj),incidentSelected=int(sum(owned)),
                 minimumIncidentSICN=float(q[adj].min()),
                 selectedAllIncident=bool(owned.all()))
        rows.append(row)
    out=BASE/'buffer01';out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'cavity.npz',vertices=x,tetrahedra=t,selected=reg['cells'],
                        boundaryFaces=topo['faces'][b],fixedFaces=topo['faces'][fixed],
                        outletFaces=topo['faces'][cap],wall=m['wall'],farfield=m['farfield'],outlet=m['outlet'])
    report=dict(passed=True,scope='named outlet buffer geometric feasibility only',
      radiusM=radius,cutXM=cut,protectedRadiusM=2*radius,selected=reg['summary'],
      selectedMinimumSICN=float(q[mask].min()),selectedBelow01=int(sum(q[mask]<.1)),
      unselectedMinimumSICN=float(q[~mask].min()),unselectedBelow01=int(sum(q[~mask]<.1)),
      outletTriangles=len(cap),fixedTriangles=len(fixed),
      fixedMinimumSingleFaceBound=float(upper.min()),fixedBoundsBelow01=int(sum(upper<.1)),
      outletMinimumSingleFaceBound=float(capupper.min()),outletBoundsBelow01=int(sum(capupper<.1)),
      hotNodes=rows,geometryOnly=True,uniqueResidualCauseProved=False,newQualifiedAreaM2=0,
      sourceSHA256={str(p):guard.sha(p) for p in [OLD/'mesh.su2',ROOT/'outputs/v240_validation/diagnostic01/report.json']})
    guard.save(out/'report.json',report);print(report)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':inspect()
    else:
        BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
        size=guard.size_bytes;guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        locks=[Path(__file__).resolve(),DOC/'PLAN.md',OLD/'mesh.su2',
          ROOT/'outputs/v240_validation/diagnostic01/report.json',
          ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/phaseEulerRegionalVolume.py',
          ROOT/'scripts/reviewSU2ExplicitQuality.py']
        raise SystemExit(0 if guard.supervise('buffer',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE) else 1)
