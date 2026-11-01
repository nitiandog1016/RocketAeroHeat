"""核对剩余坏单元的可分裂内部边及完整边星；仅制定下一候选，不增点。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from prepareCGALRemeshingRuntime import ROOT,sha
from reviewMMGLocalCandidate import EDGES,code_edges,key_cells
from cgalTetrahedralContract import faces_and_neighbors

def main():
    folder=ROOT/'outputs/v138_validation/actual';d=loadmat(str(folder/'hybrid_volume.mat'))
    old=loadmat(str(ROOT/'outputs/v136_validation/explicit/hybrid_volume.mat'));P=d['vertices'];T=d['tetrahedra'].astype(int)-1
    A=old['tetrahedra'].astype(int)-1;bad=np.array(json.loads((folder/'candidate_report.json').read_text())['badTetRowsOneBased'])-1
    sel=json.loads((ROOT/'docs/validation/v137/real_admissibility/report.json').read_text());mask=np.zeros(len(A),bool)
    for part in sel['components']:mask[np.array(part['sourceTetRowsOneBased'])-1]=True
    frozen=np.isin(key_cells(T),key_cells(A[~mask]));assert sum(frozen)==170560
    F,_=faces_and_neighbors(T);edges=np.sort(F[:,[[0,1],[1,2],[2,0]]].reshape(-1,2).astype(np.int64),axis=1)
    boundary=np.unique(edges[:,0]*len(P)+edges[:,1]);protected=np.union1d(boundary,code_edges(T[frozen],len(P)))
    eligible=np.setdiff1d(code_edges(T[bad],len(P)),protected);allEdges=code_edges(T,len(P)).reshape(-1,6)
    stars=[]
    for e in eligible:
        rows=np.flatnonzero(np.any(allEdges==e,axis=1));assert len(rows)>=3 and not np.any(frozen[rows])
        a,b=divmod(int(e),len(P));link=T[rows];opposite=[]
        for cell in link:opposite.append([int(v) for v in cell if v not in [a,b]])
        _,counts=np.unique(opposite,return_counts=True);assert np.all(counts==2),'edge link is not closed'
        stars.append(dict(vertexIDsOneBased=[a+1,b+1],incidentTetRowsOneBased=(rows+1).tolist(),
            lengthM=float(np.linalg.norm(P[a]-P[b]))))
    assert np.all(np.any(np.isin(allEdges[bad],eligible),axis=1))
    report=dict(eligibleInteriorEdges=len(eligible),allRemainingBadTetsHaveEligibleEdge=True,fullEdgeStarsAvoidFrozenTets=True,
        fullEdgeStarsHaveClosedLinks=True,minimumStarTets=min(len(r['incidentTetRowsOneBased']) for r in stars),
        maximumStarTets=max(len(r['incidentTetRowsOneBased']) for r in stars),edges=stars,
        proposedOnly=True,noNewVerticesCreated=True,noNewNativeExecution=True,
        limitation='edge eligibility does not prove that seeded cells or later optimization meet SICN 0.1',
        sourceSHA256={str(p.resolve()):sha(p) for p in [Path(__file__).resolve(),folder/'hybrid_volume.mat',
            ROOT/'scripts/reviewMMGLocalCandidate.py']})
    target=folder/'interior_seed_review.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='edges'},indent=2))
if __name__=='__main__':main()
