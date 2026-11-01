"""只读审核v136完整外域和事前两层局部选区的CGAL导入条件，不生成新网格。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from cgalTetrahedralContract import condition
from reviewSU2ExplicitQuality import quality
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    start=time.monotonic();source=ROOT/'outputs/v136_validation/explicit';out=ROOT/'docs/validation/v137/real_admissibility'
    out.mkdir(parents=True,exist_ok=False)
    paths=[source/'hybrid_volume.mat',source/'mesh_process.json',source/'geometry_audit.json',Path(__file__),
        ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'scripts/reviewSU2ExplicitQuality.py']
    locked={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    receipt=json.loads((source/'mesh_process.json').read_text());assert receipt['exitCode']==1 and receipt['sourceStable']
    assert sha(source/'hybrid_volume.mat')==receipt['outputSHA256']['hybrid_volume.mat']
    d=loadmat(str(source/'hybrid_volume.mat'));P=d['vertices'];T=d['tetrahedra'].astype(int)-1
    full,graph,_=condition(P,T);bad=np.flatnonzero(quality(P[T])<.1);assert len(bad)==205
    selected=np.zeros(len(T),bool);selected[bad]=True
    for _ in range(2):selected|=np.asarray(graph@selected.astype(int)>0).ravel()
    rows=np.flatnonzero(selected);local,_,labels=condition(P,T[rows]);parts=[]
    for c in range(local['components']):
        ids=rows[labels==c];report,_,_=condition(P,T[ids]);report.update(component=c,sourceTetRowsOneBased=(ids+1).tolist())
        parts.append(report)
    edges=np.unique(np.sort(T[rows][:,[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]].reshape(-1,2),axis=1),axis=0)
    target=float(np.median(np.linalg.norm(P[edges[:,0]]-P[edges[:,1]],axis=1)))
    report=dict(sourceSHA256=locked,sourceStableAfter={n:sha(ROOT/n)==h for n,h in locked.items()},
        fullExterior=full,selectedUnion=local,components=parts,targetEdgeLengthM=target,badTetCount=len(bad),
        permittedDirectImport=bool(full['directImportNecessaryConditions'] or all(p['directImportNecessaryConditions'] for p in parts)),
        realRemeshingExecuted=False,meshUnchanged=True,physicsQualified=False,newQualifiedArea=0,seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['components','sourceSHA256','sourceStableAfter']},indent=2))
    print('components:',[(p['tets'],p['relativeConvexHullGap']) for p in parts])
if __name__=='__main__':main()
