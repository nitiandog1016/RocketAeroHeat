"""只读v135失败原始网格：计数真实配对缺陷、全域体积和质量，不修复或重新生成。"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from su2PyramidGeometry import pyramid_geometry
from meshSU2HybridVolume import wedge_jacobian_min
from auditSU2HybridTopology import polyvolume
from prepareSU2CheckpointContinuation import ROOT,read,sha

def face_review(faces,boundary,P):
    keys,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
    j=faces.argmin(axis=1);rows=np.arange(len(faces));n=faces.shape[1]
    direction=np.bincount(inv,weights=np.sign(faces[rows,(j+1)%n]-faces[rows,(j-1)%n]))
    bad=(count>2)|((count==2)&(direction!=0));actual={tuple(f) for f in keys[count==1]};expected={tuple(f) for f in np.sort(boundary,axis=1)}
    return dict(nonManifoldFaces=int(sum(count>2)),sameDirectionSharedFaces=int(sum((count==2)&(direction!=0))),
        boundaryFaces=int(sum(count==1)),missingBoundaryFaces=len(expected-actual),unexpectedBoundaryFaces=len(actual-expected),
        badFaceNodesZeroBased=keys[bad].tolist(),badFaceCentersM=P[keys[bad]].mean(axis=1).tolist())

def main():
    start=time.monotonic();folder=ROOT/'outputs/v135_validation/actual';source=ROOT/'outputs/v133_validation/fan_workbench'
    r=read(folder/'mesh_process.json');assert r['status']=='terminal' and r['exitCode']==1 and r['sourceStable']
    for name,h in r['outputSHA256'].items():assert sha(folder/name)==h
    d=loadmat(str(folder/'hybrid_volume.mat'));old=loadmat(str(source/'hybrid_volume.mat'));raw=loadmat(str(folder/'gmsh_raw.mat'))
    P=d['vertices'];W=d['wedges'].astype(int)-1;Q=d['pyramids'].astype(int)-1;T=d['tetrahedra'].astype(int)-1
    V=d['sourceVertices'];F=d['sourceFaces'].astype(int)-1;far=d['farfield'].astype(int)-1
    prior=read(source/'mesh_report.json');h=np.array(prior['normalHeightsM']);cut=min(V[V[:,0]>=.03,0]);cap=np.all(V[F,0]>=cut,axis=1)
    retained=int(np.flatnonzero(h<=.001)[-1]);assert retained==19
    originalW=old['wedges'].astype(int)-1;mask=np.r_[~np.concatenate([cap & (j>=retained) for j in range(len(h)-1)]),np.ones(len(originalW)-(len(h)-1)*len(F),bool)]
    expectedW=originalW[mask];assert np.array_equal(P[W],old['vertices'][expectedW])
    assert np.array_equal(V,old['sourceVertices']) and np.array_equal(d['sourceFaces'],old['sourceFaces'])
    assert np.array_equal(P[:len(V)],V) and np.array_equal(P[far],old['vertices'][old['farfield'].astype(int)-1])
    pv,pi=pyramid_geometry(P,Q);prismJac=float(min(wedge_jacobian_min(P[W])));assert prismJac>0
    p=P[T];tv=np.einsum('ij,ij->i',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))/6;assert min(tv)>0
    prism=np.zeros(len(W));ref=P[W].mean(axis=1)
    for face in [[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]:
        coords=P[W[:,face]]-ref[:,None];fc=coords.mean(axis=1)
        for a in range(len(face)):
            x=coords[:,a];y=coords[:,(a+1)%len(face)]
            prism+=np.einsum('ij,ij->i',x,np.cross(y-x,fc-x))/6
    assert min(prism)>0
    tri=np.concatenate([W[:,[0,2,1]],W[:,[3,4,5]],T[:,[0,2,1]],T[:,[0,1,3]],T[:,[1,2,3]],T[:,[2,0,3]],Q[:,[0,1,4]],Q[:,[1,2,4]],Q[:,[2,3,4]],Q[:,[3,0,4]]])
    quad=np.concatenate([W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]],Q[:,[0,3,2,1]]])
    tr=face_review(tri,np.r_[F[:,[0,2,1]],far],P);qr=face_review(quad,np.empty((0,4),int),P)
    reference=np.array([[1,.5,.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2/3)]])
    J=np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=2)@np.linalg.inv(reference)
    singular=np.linalg.svd(J,compute_uv=False);quality=3/np.sqrt(np.sum(singular**2,axis=1)*np.sum(singular**-2,axis=1))
    bad=np.flatnonzero(quality<.1)
    _,counts=np.unique(np.sort(T,axis=1),axis=0,return_counts=True)
    rawT=raw['elements4'].reshape(-1,4);_,rawCounts=np.unique(np.sort(rawT,axis=1),axis=0,return_counts=True)
    assert sum(counts-1)==sum(rawCounts-1)
    # 旧17坏单元只按原顶点和连接判定；新编号不构成相同单元的证据。
    oldT=old['tetrahedra'].astype(int)-1;op=old['vertices'][oldT]
    oldJ=np.stack([op[:,1]-op[:,0],op[:,2]-op[:,0],op[:,3]-op[:,0]],axis=2)@np.linalg.inv(reference)
    s=np.linalg.svd(oldJ,compute_uv=False);oldBad=np.flatnonzero(3/np.sqrt(np.sum(s*s,axis=1)*np.sum(s**-2,axis=1))<.1);assert len(oldBad)==17
    dist,idx=cKDTree(P).query(op[oldBad].reshape(-1,3));idx=idx.reshape(-1,4);exact=np.all(dist.reshape(-1,4)==0,axis=1)
    current={tuple(row) for row in np.sort(T,axis=1)};retainedOld=sum(tuple(np.sort(row)) in current for row in idx[exact])
    total=float(sum(prism)+sum(pv)+sum(tv));expected=polyvolume(P,far)-polyvolume(V,F)
    report=dict(exitCode=r['exitCode'],generationSeconds=r['elapsedSeconds'],nodes=len(P),prisms=len(W),pyramids=len(Q),tetrahedra=len(T),
        cells=len(W)+len(Q)+len(T),originalWallAndFacesUnchanged=True,retainedPrismsExactlyUnchanged=True,farfieldExactlyUnchanged=True,
        retainedLayers=retained,retainedHeightM=float(h[retained]),capFaces=int(sum(cap)),noseCutXM=float(cut),removedPrisms=int(sum(~mask)),
        positiveCellVolumes=True,minimumPrismJacobian=prismJac,pyramidGeometry=pi,triangular=tr,quadrilateral=qr,
        summedCellVolumeM3=total,expectedFluidVolumeM3=expected,volumeClosureRelative=abs(total-expected)/expected,
        duplicateTetrahedra=int(sum(counts-1)),duplicateTetrahedraAlsoInRawGmsh=int(sum(rawCounts-1)),
        minimumTetSICN=float(min(quality)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=quality[bad].tolist(),
        badTetCenterBoundsM=[p[bad].mean(axis=1).min(axis=0).tolist(),p[bad].mean(axis=1).max(axis=0).tolist()],
        oldBadTetrahedra=17,exactlyRetainedOldBadTetrahedra=retainedOld,geometryAccepted=False,physicsQualified=False,newQualifiedArea=0,
        elapsedSeconds=time.monotonic()-start,inputSHA256={str(f.relative_to(ROOT)):sha(f) for f in [folder/'hybrid_volume.mat',folder/'gmsh_raw.mat',folder/'mesh_process.json',source/'hybrid_volume.mat',Path(__file__).resolve()]})
    out=folder/'independent_diagnosis.json';assert not out.exists();out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['triangular','quadrilateral','badTetRowsOneBased','badTetSICN','inputSHA256']},indent=2));print(json.dumps({k:{x:y for x,y in v.items() if not x.startswith('badFace')} for k,v in [('triangular',tr),('quadrilateral',qr)]}))

if __name__=='__main__':main()
