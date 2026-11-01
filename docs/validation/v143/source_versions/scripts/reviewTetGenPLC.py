"""独立只读复核PLC原生面/节点映射、逆矩阵质量及剩余坏单元位置。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from prepareTetGenRuntime import ROOT,sha
from reviewTetGenFixedRepair import read_raw,keys,independent_topology
from optimizeSU2FanExterior import sicn
from cgalTetrahedralContract import faces_and_neighbors

def main():
    start=time.monotonic();folder=ROOT/'outputs/v143_validation/actual';src=ROOT/'outputs/v141_validation/actual'
    paths=[Path(__file__).resolve(),ROOT/'scripts/reviewTetGenFixedRepair.py',ROOT/'scripts/optimizeSU2FanExterior.py',
        ROOT/'scripts/cgalTetrahedralContract.py',folder/'hybrid_volume.mat',folder/'candidate_report.json',folder/'mesh_process.json',
        folder/'native/input.dat',folder/'native/native_import.dat',folder/'native/native_result.dat',src/'hybrid_volume.mat',
        ROOT/'docs/validation/v143/hollow_test/report.json']
    locked={str(p):sha(p) for p in paths};r=json.loads((folder/'mesh_process.json').read_text())
    assert r['exitCode']==1 and r['sourceStable']
    for n,h in r['outputSHA256'].items():assert sha(folder/n)==h
    for n,h in r['sourceSHA256'].items():assert sha(Path(n))==h
    small=ROOT/'docs/validation/v143/hollow_test';sp,st,sf=read_raw(small/'native/native_result.dat')
    sb,sv=independent_topology(sp[:,:3],st);sq=sicn(sp[:,:3],st)
    assert abs(sv-26)<1e-12 and min(sq)>=.1 and sorted(sb)==keys(sf[:,:3])
    assert (small/'native/input.dat').read_bytes()==(small/'native/native_import.dat').read_bytes()
    # 输入与导入文本采用相同数字值，但不依赖空白相同作真实坐标对应。
    lines=(folder/'native/input.dat').read_text().splitlines();imported=(folder/'native/native_import.dat').read_text().splitlines()
    assert [s.split() for s in lines]==[s.split() for s in imported]
    old=loadmat(str(src/'hybrid_volume.mat'));d=loadmat(str(folder/'hybrid_volume.mat'));P=old['vertices'];p=d['vertices']
    Q,T,F=read_raw(folder/'native/native_result.dat');refs=Q[:,3].astype(np.int64);fixed=refs>0
    srcB,_=faces_and_neighbors(old['tetrahedra'].astype(np.int64)-1);boundaryIDs=np.unique(srcB)
    assert np.array_equal(d['sourcePLCBoundaryNodeIDsOneBased'].ravel()-1,boundaryIDs)
    assert np.array_equal(Q[fixed,:3],P[boundaryIDs[refs[fixed]-1]])
    nativeToPLC=np.empty(len(Q),np.int64);nativeToPLC[fixed]=refs[fixed]-1
    nativeToPLC[~fixed]=np.arange(len(boundaryIDs),len(boundaryIDs)+sum(~fixed),dtype=np.int64)
    plcToHybrid=d['plcToHybridNodeIDsOneBased'].ravel().astype(np.int64)-1;nativeToHybrid=plcToHybrid[nativeToPLC]
    assert np.array_equal(p[nativeToHybrid],Q[:,:3])
    C=d['tetrahedra'].astype(np.int64)-1;assert np.array_equal(C,nativeToHybrid[T])
    origin=d['sourceV141NodeIDsOneBased'].ravel().astype(np.int64)-1;retained=origin>=0
    assert np.array_equal(p[retained],P[origin[retained]])
    assert np.array_equal(origin[d['wedges']-1]+1,old['wedges']) and np.array_equal(origin[d['farfield']-1]+1,old['farfield'])
    assert np.array_equal(d['sourceVertices'],old['sourceVertices']) and np.array_equal(d['sourceFaces'],old['sourceFaces'])
    assert np.array_equal(d['normalHeightsByColumnM'],old['normalHeightsByColumnM'])
    boundary,volume=independent_topology(Q[:,:3],T);assert sorted(boundary)==keys(F[:,:3])
    assert keys(origin[nativeToHybrid[F[:,:3]]])==keys(srcB)
    val=sicn(p,C);bad=np.flatnonzero(val<.1);candidate=json.loads((folder/'candidate_report.json').read_text())
    assert np.array_equal(bad+1,candidate['badTetRowsOneBased']);error=float(max(abs(val[bad]-candidate['badTetSICN'])));assert error<1e-12
    assert abs(volume/candidate['expectedFluidVolume']-1)<2e-10
    # 模块2：仅定位剩余失败，不运行优化或以单面必要界证明全域可行。
    globalBoundary=np.unique(nativeToHybrid[F[:,:3]]);counts=np.isin(C[bad],globalBoundary).sum(axis=1)
    facekeys=set(keys(nativeToHybrid[F[:,:3]]));faceCounts=[];minAngle=[];maxAngle=[];radiusRatio=[]
    for t in p[C[bad]]:
        faces=t[np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])];normal=np.cross(faces[:,1]-faces[:,0],faces[:,2]-faces[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        angle=np.degrees(np.arccos(np.clip(-(normal@normal.T)[np.triu_indices(4,1)],-1,1)))
        minAngle.append(float(min(angle)));maxAngle.append(float(max(angle)))
        e=t[1:]-t[0];center=np.linalg.solve(2*e,np.sum(e*e,axis=1))
        edges=np.array([np.linalg.norm(t[i]-t[j]) for i in range(4) for j in range(i)])
        radiusRatio.append(float(np.linalg.norm(center)/min(edges)))
    for row in C[bad]:faceCounts.append(sum(tuple(sorted(row[list(idx)])) in facekeys for idx in [(0,2,1),(0,1,3),(1,2,3),(2,0,3)]))
    centers=p[C[bad]].mean(axis=1);farIDs=np.unique(d['farfield']-1)
    result=dict(passed=True,seconds=time.monotonic()-start,smallIndependentVolume=sv,smallIndependentMinimumSICN=float(min(sq)),
        allRawNativeCoordinatesAndConnectivityExactlyMapped=True,allOriginalPrismsWallFarfieldAndHeightsPreserved=True,
        explicitHolePreserved=True,independentTetVolume=volume,minimumSICN=float(min(val)),independentBadSICNDifference=error,
        badTets=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetNodeIDsOneBased=(C[bad]+1).tolist(),badTetCentersM=centers.tolist(),
        badTetBoundaryVertexCounts=counts.tolist(),badTetBoundaryFaceCounts=faceCounts,badTetFarfieldVertexCounts=np.isin(C[bad],farIDs).sum(axis=1).tolist(),
        badTetMinDihedralDeg=minAngle,badTetMaxDihedralDeg=maxAngle,badTetRadiusEdgeRatio=radiusRatio,
        candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,sourceSHA256=locked,
        sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()})
    assert all(result['sourceStableAfter'].values());dst=folder/'independent_review.json';assert not dst.exists()
    dst.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['sourceSHA256','sourceStableAfter','badTetNodeIDsOneBased']},indent=2))

if __name__=='__main__':main()
