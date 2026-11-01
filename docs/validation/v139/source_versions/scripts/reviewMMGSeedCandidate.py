"""只读独立复核v139细分谱系、固定域与逆条件数；为下一几何路线定位失败。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareCGALRemeshingRuntime import ROOT,sha
from optimizeSU2FanExterior import sicn
from mmgLocalInterface import sorted_cells,read_output
from cgalTetrahedralContract import faces_and_neighbors

def volume(P,T):
    q=P[T];return np.linalg.det(np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0],q[:,3]-q[:,0]],axis=2))/6

def keys(cells):
    a=np.ascontiguousarray(np.sort(cells.astype(np.int64),axis=1));return a.view(np.dtype((np.void,8*cells.shape[1]))).ravel()

def main():
    start=time.monotonic();folder=ROOT/'outputs/v139_validation/actual';source=ROOT/'outputs/v138_validation/actual'
    paths=[Path(__file__).resolve(),ROOT/'scripts/optimizeSU2FanExterior.py',ROOT/'scripts/mmgLocalInterface.py',
        ROOT/'scripts/cgalTetrahedralContract.py',source/'hybrid_volume.mat',folder/'hybrid_volume.mat',folder/'seeded_volume.mat',
        folder/'mesh_process.json',folder/'candidate_report.json',folder/'split_report.json',folder/'native/native_result.dat',
        ROOT/'docs/validation/v139/seed_test/fixture.npz',ROOT/'docs/validation/v139/seed_test/report.json']
    locked={str(p):sha(p) for p in paths};receipt=json.loads((folder/'mesh_process.json').read_text())
    assert receipt['exitCode']==1 and receipt['sourceStable'] and not receipt['timedOut']
    for name,h in receipt['outputSHA256'].items():assert sha(folder/name)==h
    old=loadmat(str(source/'hybrid_volume.mat'));seed=loadmat(str(folder/'seeded_volume.mat'));new=loadmat(str(folder/'hybrid_volume.mat'))
    P=old['vertices'];A=old['tetrahedra'].astype(np.int64)-1;S=seed['vertices'];B=seed['tetrahedra'].astype(np.int64)-1
    Q=new['vertices'];C=new['tetrahedra'].astype(np.int64)-1;parent=seed['parentTetRowsOneBased'].ravel().astype(np.int64)-1
    mask=seed['selectedTets'].ravel().astype(bool);events=json.loads((folder/'split_report.json').read_text())['events']
    assert np.array_equal(S[:len(P)],P) and np.array_equal(Q[:len(P)],P)
    for n in ['sourceVertices','sourceFaces','wedges','pyramids','farfield']:
        assert np.array_equal(old[n],seed[n]) and np.array_equal(old[n],new[n])
    # 独立行列式积分检查每个原父单元，不能只看全域相互抵消的体积误差。
    va=volume(P,A);vb=volume(S,B);vc=volume(Q,C)
    aggregate=np.bincount(parent,weights=vb,minlength=len(A));err=float(max(abs(aggregate-va)/va));assert err<2e-10
    assert min(vb)>0 and min(vc)>0
    for event in events:
        a,b=[v-1 for v in event['edgeIDsOneBased']];node=event['newVertexIDOneBased']-1
        assert np.array_equal(S[node],(P[a]+P[b])/2)
    frozen=B[~mask];assert len(frozen)==170560 and np.all(np.isin(keys(frozen),keys(C)))
    boundary,_=faces_and_neighbors(A);newboundary,_=faces_and_neighbors(C)
    assert np.array_equal(sorted_cells(boundary),sorted_cells(newboundary))
    q=sicn(Q,C);r=json.loads((folder/'candidate_report.json').read_text());bad=np.flatnonzero(q<.1)
    assert np.array_equal(bad+1,r['badTetRowsOneBased']);qerror=float(max(abs(q[bad]-r['badTetSICN'])));assert qerror<1e-12
    seedq=sicn(S,B);badC=C[bad];physical=np.unique(boundary);counts=np.isin(badC,physical).sum(axis=1)
    sides=new['pyramids'].astype(np.int64)-1;sides=sides[:,[[0,1,4],[1,2,4],[2,3,4],[3,0,4]]].reshape(-1,3)
    bf=badC[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
    touchesPyramid=np.isin(keys(bf),keys(sides)).reshape(-1,4).any(axis=1)
    touchesBoundary=np.isin(keys(bf),keys(boundary)).reshape(-1,4).sum(axis=1)
    cut=json.loads((source/'mesh_report.json').read_text())['noseCutSourceXM'];center=Q[badC].mean(axis=1)
    raw=read_output(folder/'native/native_result.dat');free=raw[0][raw[3]==0];assert len(free)==247 and len(Q)-len(P)==247
    distance,_=cKDTree(S[len(P):]).query(free)
    # 小例用同一独立逆矩阵公式核对质量和解析体积；不重跑原生程序。
    small=ROOT/'docs/validation/v139/seed_test';fix=np.load(str(small/'fixture.npz'));smallraw=read_output(small/'native/native_result.dat')
    smallq=sicn(smallraw[0],smallraw[1]);sm=json.loads((small/'report.json').read_text())
    assert abs(min(smallq)-sm['minimumSICNAfter'])<1e-12 and abs(sum(volume(smallraw[0],smallraw[1]))-sm['knownVolume'])<2e-12
    stable={p:sha(Path(p))==h for p,h in locked.items()};assert all(stable.values())
    report=dict(passed=True,sourceSHA256=locked,sourceStableAfter=stable,seconds=time.monotonic()-start,
        allOriginalArraysExactlyPreserved=True,frozenTetsExactlyPreserved=len(frozen),seedMidpointsExactlyVerified=len(events),
        originalParentVolumeMaxRelative=err,boundaryExactlyPreserved=True,independentSICNMaxDifference=qerror,
        initialBadTets=int(sum(sicn(P,A)<.1)),seededBadTets=int(sum(seedq<.1)),finalBadTets=len(bad),
        seededMinimumSICN=float(min(seedq)),finalMinimumSICN=float(min(q)),
        badTetPhysicalBoundaryVertexHistogram=np.bincount(counts,minlength=5).tolist(),
        badTetPhysicalBoundaryFaceHistogram=np.bincount(touchesBoundary,minlength=5).tolist(),
        badTetsSharingPyramidFace=int(sum(touchesPyramid)),badTetsContainingInsertedNode=int(sum(np.any(badC>=len(P),axis=1))),
        badTetCenterXRangeM=[float(min(center[:,0])),float(max(center[:,0]))],noseCutSourceXM=cut,
        badTetCentersAtOrAboveNoseCut=int(sum(center[:,0]>=cut)),badTetCentersBelowNoseCut=int(sum(center[:,0]<cut)),
        freeOutputNodes=len(free),freePointsNotAtAnySeedCoordinate=int(sum(distance>1e-12)),maximumNearestSeedDistanceM=float(max(distance)),
        displacementScope='unlabelled free-point set difference, not per-node trajectories',
        smallIndependentQualityAndVolumePassed=True,candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,
        limitation='failure of this bounded candidate only, not proof that any fixed-interface mesh is impossible')
    target=folder/'independent_review.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','sourceStableAfter']},indent=2))

if __name__=='__main__':main()
