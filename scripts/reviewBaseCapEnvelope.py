"""独立复核底面外层切向重分布、保护区域、完整原生映射和全部四面体质量。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from prepareTetGenRuntime import ROOT,sha
from reviewTetGenFixedRepair import read_raw,keys,independent_topology
from optimizeSU2FanExterior import sicn
from reviewSU2ExplicitQuality import quality

def inverse_ids(ids):
    good=ids>=0;out=np.full(int(max(ids))+1,-1,np.int64);out[ids[good]]=np.flatnonzero(good);return out

def main():
    start=time.monotonic();folder=ROOT/'outputs/v144_validation/actual';source=ROOT/'outputs/v143_validation/actual'
    fanPath=ROOT/'outputs/v133_validation/fan_workbench/mesh_report.json'
    paths=[Path(__file__).resolve(),ROOT/'scripts/reviewTetGenFixedRepair.py',ROOT/'scripts/optimizeSU2FanExterior.py',
        ROOT/'scripts/reviewSU2ExplicitQuality.py',fanPath,source/'hybrid_volume.mat']
    paths += [folder/n for n in ['mesh_process.json','prepared_cap.mat','hybrid_volume.mat','candidate_report.json','preparation_report.json',
        'boundary_audit/audit_input.mat','boundary_audit/process.json','boundary_audit/intersections.json',
        'native/input.dat','native/native_import.dat','native/native_result.dat','native/process/process.json']]
    locked={str(p):sha(p) for p in paths};receipt=json.loads((folder/'mesh_process.json').read_text())
    assert receipt['exitCode']==0 and receipt['sourceStable']
    for n,h in receipt['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in receipt['outputSHA256'].items():assert sha(folder/n)==h
    nativeReceipt=json.loads((folder/'native/process/process.json').read_text());assert nativeReceipt['exitCode']==0
    for n,h in nativeReceipt['sourceSHA256'].items():assert sha(Path(n))==h
    old=loadmat(str(source/'hybrid_volume.mat'));prepared=loadmat(str(folder/'prepared_cap.mat'));new=loadmat(str(folder/'hybrid_volume.mat'))
    P=old['vertices'];Q=prepared['vertices'];V=old['sourceVertices'];F=old['sourceFaces'].astype(np.int64)-1;nv=len(V)
    oldPrepared=old['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1;inverse=inverse_ids(oldPrepared)
    col=int(prepared['baseColumnSourceIDOneBased'].ravel()[0])-1;assert col==9876
    ids=inverse[np.arange(33)*nv+col];fan=json.loads(fanPath.read_text());edge=np.array(fan['rimEdges'],np.int64)
    outeredge=inverse[edge+32*nv];a=P[outeredge[:,0],1:];b=P[outeredge[:,1],1:]
    # 独立分扇面积积分求环质心，与生成器的有向边shoelace公式分开。
    origin=P[inverse[np.array(fan['rimNodeIDs'])+32*nv],1:].mean(axis=0)
    aa=a-origin;bb=b-origin;areas=abs(aa[:,0]*bb[:,1]-aa[:,1]*bb[:,0])/2
    centroid=np.sum(areas[:,None]*(origin+a+b)/3,axis=0)/sum(areas)
    h=old['normalHeightsByColumnM'][:,col];N=(P[ids[-1]]-P[ids[0]])/h[-1]
    tangent=np.r_[0.,centroid-P[ids[-1],1:]];tangent-=N*np.dot(N,tangent)
    u=np.clip((h-h[19])/(h[-1]-h[19]),0,1);weights=6*u**5-15*u**4+10*u**3
    offset=weights[:,None]*tangent;expected=P.copy();expected[ids[20:]]+=offset[20:]
    coordinateError=float(np.max(abs(Q-expected)));assert coordinateError<2e-14
    offsetError=float(np.max(abs(prepared['baseColumnTangentialOffsetsM']-offset)));assert offsetError<2e-14
    moved=np.flatnonzero(np.any(Q!=P,axis=1));assert np.array_equal(moved,np.sort(ids[20:]))
    protected=inverse[np.arange(20*nv)];assert min(protected)>=0 and np.array_equal(Q[protected],P[protected])
    fanSource=np.unique(np.r_[np.array(fan['sideNodeIDs']).ravel(),np.array(fan['baseNodeIDs']).ravel()]);fanIDs=inverse[fanSource]
    assert not np.any(np.isin(ids,fanIDs)) and np.array_equal(Q[fanIDs],P[fanIDs])
    assert np.array_equal(prepared['wedges'],old['wedges']) and np.array_equal(Q[old['wedges'][-2852:]-1],P[old['wedges'][-2852:]-1])
    assert np.array_equal(Q[old['farfield']-1],P[old['farfield']-1]) and np.array_equal(Q[:nv],V)
    assert np.array_equal(prepared['normalHeightsByColumnM'],old['normalHeightsByColumnM'])
    # 完整CGAL证据必须对应此次实际新外壳，不能继承旧表面结论。
    cg=json.loads((folder/'boundary_audit/process.json').read_text());assert cg['exitCode']==0 and cg['runtimeStableAfter']
    for n,hv in cg['outputSHA256'].items():assert sha(folder/'boundary_audit'/n)==hv
    assert json.loads((folder/'boundary_audit/intersections.json').read_text())['selfIntersectionFree']
    ci=loadmat(str(folder/'boundary_audit/audit_input.mat'));top=prepared['farfield'].astype(np.int64)-1;used,back=np.unique(top,return_inverse=True)
    assert np.array_equal(ci['V'],Q[used]) and np.array_equal(ci['F']-1,back.reshape(top.shape))
    # 原生全部点/连接按标记映射到最终网格；旧点中13点改变已单独解释。
    lines=(folder/'native/input.dat').read_text().splitlines();imported=(folder/'native/native_import.dat').read_text().splitlines()
    assert [s.split() for s in lines]==[s.split() for s in imported]
    raw,T,G=read_raw(folder/'native/native_result.dat');ref=raw[:,3].astype(np.int64);fixed=ref>0
    boundaryIDs=new['sourcePLCBoundaryNodeIDsOneBased'].ravel().astype(np.int64)-1
    assert np.array_equal(raw[fixed,:3],Q[boundaryIDs[ref[fixed]-1]])
    rawToLocal=np.empty(len(raw),np.int64);rawToLocal[fixed]=ref[fixed]-1;rawToLocal[~fixed]=np.arange(len(boundaryIDs),len(boundaryIDs)+sum(~fixed))
    localToHybrid=new['plcToHybridNodeIDsOneBased'].ravel().astype(np.int64)-1;nativeToHybrid=localToHybrid[rawToLocal]
    p=new['vertices'];C=new['tetrahedra'].astype(np.int64)-1
    assert np.array_equal(p[nativeToHybrid],raw[:,:3]) and np.array_equal(C,nativeToHybrid[T])
    sourceIDs=new['sourceV143NodeIDsOneBased'].ravel().astype(np.int64)-1;retained=sourceIDs>=0
    assert np.array_equal(p[retained],Q[sourceIDs[retained]])
    assert np.array_equal(sourceIDs[new['wedges']-1]+1,old['wedges']) and np.array_equal(sourceIDs[new['farfield']-1]+1,old['farfield'])
    assert np.array_equal(new['sourceVertices'],V) and np.array_equal(new['sourceFaces'],F+1)
    assert np.array_equal(new['baseColumnTangentialOffsetsM'],prepared['baseColumnTangentialOffsetsM'])
    B,volume=independent_topology(raw[:,:3],T);assert sorted(B)==keys(G[:,:3])
    expectedBoundary=np.r_[top[:,[0,2,1]],old['farfield']-1]
    assert keys(sourceIDs[nativeToHybrid[G[:,:3]]])==keys(expectedBoundary)
    q=sicn(p,C);svd=quality(p[C]);error=float(max(abs(q-svd)));assert error<1e-12 and min(q)>=.1
    report=json.loads((folder/'candidate_report.json').read_text());assert report['candidateAccepted'] and report['tetsBelowSICN01']==0
    assert abs(min(q)-report['minimumTetSICN'])<1e-12 and abs(volume/report['expectedFluidVolume']-1)<2e-10
    face=V[F];area=float(np.linalg.norm(np.cross(face[:,1]-face[:,0],face[:,2]-face[:,0]),axis=1).sum()/2)
    assert area==report['originalWallAreaM2']
    stable={p:sha(Path(p))==hv for p,hv in locked.items()};assert all(stable.values())
    result=dict(passed=True,seconds=time.monotonic()-start,independentCoordinateDifferenceM=coordinateError,
        independentOffsetDifferenceM=offsetError,movedOriginalNodes=len(moved),
        originalWallAreaM2=area,protected19LayersAndEntireTailFanPreserved=True,allOriginalPrismConnectivityPreserved=True,
        actualNewCGALSurfaceVerified=True,allRawCoordinatesAndConnectivityExactlyMapped=True,
        independentTetVolume=volume,minimumSICN=float(min(q)),independentFullSICNMaximumDifference=error,tetsBelowSICN01=int(sum(q<.1)),
        maximumOuterColumnTangentialSlope=float(max(np.linalg.norm(np.diff(offset,axis=0),axis=1)/np.diff(h))),
        geometryCandidateAccepted=True,flowAccuracyQualified=False,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256=locked,sourceStableAfter=stable)
    dst=folder/'independent_review.json';assert not dst.exists();dst.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['sourceSHA256','sourceStableAfter']},indent=2))

if __name__=='__main__':main()
