"""独立复核MMG候选质量与局部冻结约束；不重复重网格或改变拒绝身份。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from prepareCGALRemeshingRuntime import ROOT,sha
from cgalTetrahedralContract import faces_and_neighbors
from optimizeSU2FanExterior import sicn

EDGES=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
def code_edges(T,n):
    e=np.sort(T[:,EDGES].reshape(-1,2).astype(np.int64),axis=1);return e[:,0]*n+e[:,1]
def key_cells(T):
    a=np.ascontiguousarray(np.sort(T,axis=1));return a.view(np.dtype((np.void,a.dtype.itemsize*4))).ravel()

def main():
    source=ROOT/'outputs/v136_validation/explicit';folder=ROOT/'outputs/v138_validation/actual'
    p=folder/'mesh_process.json';receipt=json.loads(p.read_text());assert receipt['exitCode']==1 and receipt['sourceStable'] and not receipt['timedOut']
    for n,h in receipt['outputSHA256'].items():assert sha(folder/n)==h
    old=loadmat(str(source/'hybrid_volume.mat'));new=loadmat(str(folder/'hybrid_volume.mat'))
    P=old['vertices'];A=old['tetrahedra'].astype(int)-1;B=new['tetrahedra'].astype(int)-1
    assert np.array_equal(P,new['vertices'])
    for n in ['wedges','pyramids','farfield','sourceVertices','sourceFaces']:assert np.array_equal(old[n],new[n])
    q=sicn(P,B);report=json.loads((folder/'candidate_report.json').read_text());bad=np.flatnonzero(q<.1)
    assert np.array_equal(bad+1,report['badTetRowsOneBased']);err=float(max(abs(q[bad]-report['badTetSICN'])));assert err<1e-12
    sel=json.loads((ROOT/'docs/validation/v137/real_admissibility/report.json').read_text());selected=np.zeros(len(A),bool)
    for c in sel['components']:selected[np.array(c['sourceTetRowsOneBased'])-1]=True
    boundary,_=faces_and_neighbors(A);edges=np.sort(boundary[:,[[0,1],[1,2],[2,0]]].reshape(-1,2).astype(np.int64),axis=1)
    phys=np.unique(edges[:,0]*len(P)+edges[:,1]);frozen=np.unique(code_edges(A[~selected],len(P)))
    ee=code_edges(B[bad],len(P));physcount=np.isin(ee,phys).reshape(-1,6).sum(axis=1)
    frozenCount=np.isin(ee,frozen).reshape(-1,6).sum(axis=1);either=np.isin(ee,np.union1d(phys,frozen)).reshape(-1,6).sum(axis=1)
    boundaryNodes=np.unique(boundary);nodecount=np.isin(B[bad],boundaryNodes).sum(axis=1)
    oldkeys=key_cells(A);newkeys=key_cells(B);removed=~np.isin(oldkeys,newkeys);created=~np.isin(newkeys,oldkeys)
    centers=P[B[bad]].mean(axis=1)
    data=dict(passed=True,originalNodesExactlyEqual=True,allPrismPyramidOriginalWallFarfieldArraysEqual=True,
        independentQualityMaximumDifference=err,badTets=len(bad),minimumSICN=float(min(q)),
        removedOriginalTets=int(sum(removed)),createdNewTets=int(sum(created)),
        retainedOriginalBadTets=int(sum(np.isin(newkeys[bad],oldkeys))),
        badTetPhysicalBoundaryVertexHistogram=np.bincount(nodecount,minlength=5).tolist(),
        badTetPhysicalBoundaryEdgeHistogram=np.bincount(physcount,minlength=7).tolist(),
        badTetFrozenNeighborEdgeHistogram=np.bincount(frozenCount,minlength=7).tolist(),
        badTetUnionProtectedEdgeHistogram=np.bincount(either,minlength=7).tolist(),
        badTetsWithAllSixEdgesProtected=int(sum(either==6)),noseBadTets=int(sum(centers[:,0]>0)),
        physicalBoundaryNodes=len(boundaryNodes),allNativeNodesFrozen=True,
        scope='constraints of this fixed selection; not mathematical impossibility of boundary-preserving remeshing',
        physicsQualified=False,newQualifiedArea=0,sourceSHA256={str(p.resolve()):sha(p) for p in [Path(__file__).resolve(),
            ROOT/'scripts/optimizeSU2FanExterior.py',source/'hybrid_volume.mat',folder/'hybrid_volume.mat',folder/'mesh_process.json']})
    target=folder/'independent_review02.json';assert not target.exists();target.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))
if __name__=='__main__':main()
