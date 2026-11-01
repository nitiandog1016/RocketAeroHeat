"""独立复核连续包络的层高、原保护点、原生连接与剩余坏四面体边界约束。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from prepareCGALRemeshingRuntime import ROOT,sha
from optimizeSU2FanExterior import sicn
from reviewSU2ExplicitQuality import bound
from cgalTetrahedralContract import faces_and_neighbors
from mmgLocalInterface import sorted_cells

def main():
    start=time.monotonic();folder=ROOT/'outputs/v141_validation/actual';source=ROOT/'outputs/v133_validation/fan_workbench'
    paths=[Path(__file__).resolve(),ROOT/'scripts/optimizeSU2FanExterior.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',
        ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'scripts/mmgLocalInterface.py',source/'hybrid_volume.mat',
        folder/'prepared_envelope.mat',folder/'hybrid_volume.mat',folder/'gmsh_raw.mat',folder/'mesh_process.json',
        folder/'candidate_report.json',folder/'boundary_audit/audit_input.mat',folder/'boundary_audit/process.json']
    locked={str(p):sha(p) for p in paths};receipt=json.loads((folder/'mesh_process.json').read_text());assert receipt['exitCode']==1 and receipt['sourceStable']
    for n,h in receipt['outputSHA256'].items():assert sha(folder/n)==h
    old=loadmat(str(source/'hybrid_volume.mat'));prepared=loadmat(str(folder/'prepared_envelope.mat'));new=loadmat(str(folder/'hybrid_volume.mat'));raw=loadmat(str(folder/'gmsh_raw.mat'))
    p=old['vertices'];P=prepared['vertices'];V=old['sourceVertices'];nv=len(V);F=old['sourceFaces'].astype(np.int64)-1
    h=prepared['originalNormalHeightsM'].ravel();H=prepared['normalHeightsByColumnM'];meta=json.loads((folder/'preparation_report.json').read_text())
    u=np.maximum(0,np.minimum(1,(V[:,0]-meta['joinSourceXM'])/(meta['cutSourceXM']-meta['joinSourceXM'])))
    weight=6*u**5-15*u**4+10*u**3;expectedH=np.tile(h[:,None],(1,nv))
    for j in range(20,33):expectedH[j]=h[j]+weight*(h[19]+(j-19)*(h[19]-h[18])-h[j])
    heightError=float(np.max(abs(H-expectedH)));assert heightError<2e-16
    chosen=np.flatnonzero(weight>0);normal=(p[32*nv+chosen]-V[chosen])/h[32];expected=p.copy()
    for j in range(20,33):expected[j*nv+chosen]=p[j*nv+chosen]+(expectedH[j,chosen]-h[j])[:,None]*normal
    pointError=float(np.max(abs(P-expected)));assert pointError<2e-13
    assert np.array_equal(P[:20*nv],p[:20*nv]) and np.array_equal(prepared['wedges'],old['wedges'])
    # 逐列外层增厚比不低于1，不以把13层挤成极薄首层获得小包络。
    spacing=np.diff(H,axis=0);growth=spacing[19:]/spacing[18:-1];assert min(growth.ravel())>1-2e-12
    ids=new['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1;original=ids<len(P)
    assert np.array_equal(new['vertices'][original],P[ids[original]])
    W=prepared['wedges'].astype(np.int64)-1;U=new['wedges'].astype(np.int64)-1
    assert np.array_equal(ids[U],W) and np.array_equal(P[W[:19*len(F)]],p[W[:19*len(F)]]) and np.array_equal(P[W[-2852:]],p[W[-2852:]])
    assert np.array_equal(new['sourceFaces'],old['sourceFaces']) and np.array_equal(new['sourceVertices'],V)
    assert np.array_equal(ids[new['farfield']-1]+1,old['farfield'])
    fixed=np.unique(np.r_[prepared['farfield'].ravel()-1,prepared['exteriorFarfield'].ravel()-1]);tags=raw['nodeTags'].ravel().astype(np.int64)
    adds=~np.isin(tags,fixed+1);assert np.array_equal(new['vertices'][~original],raw['vertices'][adds])
    mapRaw=np.full(int(max(tags))+1,-1,np.int64);mapRaw[fixed+1]=fixed;mapRaw[tags[adds]]=np.arange(len(P),len(P)+sum(adds))
    C=new['tetrahedra'].astype(np.int64)-1;assert np.array_equal(ids[C],mapRaw[raw['elements4'].reshape(-1,4)])
    values=sicn(new['vertices'],C);r=json.loads((folder/'candidate_report.json').read_text());bad=np.flatnonzero(values<.1)
    assert np.array_equal(bad+1,r['badTetRowsOneBased']);err=float(max(abs(values[bad]-r['badTetSICN'])));assert err<1e-12
    boundary,_=faces_and_neighbors(C);bnodes=np.unique(boundary);counts=np.isin(C[bad],bnodes).sum(axis=1)
    boundaryKeys={tuple(v) for v in np.sort(boundary,axis=1)}
    faces=C[bad][:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]];faceCounts=np.array([sum(tuple(sorted(f)) in boundaryKeys for f in row) for row in faces])
    upper,_,_,_=bound(new['vertices'][boundary]);center=new['vertices'][C[bad]].mean(axis=1)
    auditInput=loadmat(str(folder/'boundary_audit/audit_input.mat'));top=prepared['farfield'].astype(np.int64)-1;used,inverse=np.unique(top,return_inverse=True)
    assert np.array_equal(auditInput['V'],P[used]) and np.array_equal(auditInput['F']-1,inverse.reshape(top.shape))
    stable={n:sha(Path(n))==v for n,v in locked.items()};assert all(stable.values())
    report=dict(passed=True,seconds=time.monotonic()-start,sourceSHA256=locked,sourceStableAfter=stable,
        independentHeightMaximumDifferenceM=heightError,independentCoordinateMaximumDifferenceM=pointError,
        minimumOuterLayerGrowth=float(min(growth.ravel())),maximumOuterLayerGrowth=float(max(growth.ravel())),
        protected19LayerPointsExactlyPreserved=True,tailFanAndOriginalWallFarfieldExactlyPreserved=True,
        allPrismConnectivityExactlyPreserved=True,allNativeCellsAndCoordinatesExactlyMatched=True,cgalActualNewEnvelopeVerified=True,
        independentSICNMaximumDifference=err,badTets=len(bad),minimumSICN=float(min(values)),
        badTetBoundaryVertexHistogram=np.bincount(counts,minlength=5).tolist(),badTetBoundaryFaceHistogram=np.bincount(faceCounts,minlength=5).tolist(),
        badTetCenterXRangeM=[float(min(center[:,0])),float(max(center[:,0]))],
        minimumBoundarySingleFaceOptimalSICN=float(min(upper)),boundaryFacesWithOptimalSICNBelow01=int(sum(upper<.1)),
        singleFaceBoundScope='necessary single-fixed-triangle constraint only; not a global tetrahedralization existence proof',
        candidateAccepted=False,physicsQualified=False,newQualifiedArea=0)
    target=folder/'independent_review.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','sourceStableAfter']},indent=2))

if __name__=='__main__':main()
