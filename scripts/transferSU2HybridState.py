"""在真实旧四面体中定位新外层点并插值守恒初猜；重合节点严格保留，不作物理资格判断。"""
from pathlib import Path
import argparse,hashlib,json,os
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def transfer(source_mesh,target_mesh,restart,output):
    assert not output.exists();a=loadmat(str(source_mesh));b=loadmat(str(target_mesh))
    A=a['vertices'];B=b['vertices'];T=a['tetrahedra'].astype(int)-1
    state=np.loadtxt(str(restart),delimiter=',',skiprows=1)
    assert np.array_equal(state[:,0],np.arange(len(A))) and np.max(abs(state[:,1:4]-A))<1e-14
    assert np.array_equal(a['sourceVertices'],b['sourceVertices']) and np.array_equal(a['sourceFaces'],b['sourceFaces'])
    distance,old=cKDTree(A).query(B);same=distance<1e-14
    # 数值格式造成的1e-14以内匹配只用于精确坐标同一节点，不把近点配对当成重合。
    same &= np.all(A[old]==B,axis=1)
    q=state[:,4:9];result=np.empty((len(B),5));result[same]=q[old[same]]
    index=np.flatnonzero(~same);remaining=index.copy();p=A[T];centroid=p.mean(axis=1)
    inv=np.linalg.inv(np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=2));tree=cKDTree(centroid)
    found=np.full(len(B),-1,dtype=int);weight=np.zeros((len(B),4));maximumCandidates=0
    for k in [16,64,256,1024]:
        if not len(remaining):break
        maximumCandidates=k;unfound=[]
        for start in range(0,len(remaining),128):
            ids=remaining[start:start+128];_,near=tree.query(B[ids],k=min(k,len(T)))
            local=np.einsum('nkij,nkj->nki',inv[near],B[ids,None,:]-p[near,0])
            w=np.concatenate([1-local.sum(axis=2,keepdims=True),local],axis=2)
            inside=(w.min(axis=2)>=-2e-11)&(w.max(axis=2)<=1+2e-11)
            hit=inside.any(axis=1);pick=np.argmax(inside,axis=1);good=ids[hit];selected=near[np.arange(len(ids)),pick][hit]
            found[good]=selected;weight[good]=w[np.arange(len(ids)),pick][hit];unfound.extend(ids[~hit])
        remaining=np.array(unfound,dtype=int)
    assert not len(remaining),'新点未被旧真实外域覆盖；不使用最近点或静默外推'
    w=weight[index];ids=T[found[index]]
    assert np.max(abs(w.sum(axis=1)-1))<1e-14
    coordinateError=np.max(abs(np.einsum('ni,nij->nj',w,A[ids])-B[index]))
    assert coordinateError<1e-12
    result[index]=np.einsum('ni,nij->nj',w,q[ids])
    rho=result[:,0];internal=result[:,4]/rho-.5*np.sum((result[:,1:4]/rho[:,None])**2,axis=1)
    assert np.isfinite(result).all() and (rho>0).all() and (internal>0).all(),'非正状态不能靠裁剪通过'
    assert np.array_equal(result[same],q[old[same]])
    nv=len(a['sourceVertices']);assert same[:25*nv].all() and np.array_equal(old[:25*nv],np.arange(25*nv))
    values=np.column_stack([np.arange(len(B)),B,result]);header='"PointID","x","y","z","Density","Momentum_x","Momentum_y","Momentum_z","Energy"'
    np.savetxt(str(output),values,delimiter=',',header=header,comments='',fmt=['%d']+['%.17e']*8)
    report=dict(originalPoints=len(A),newPoints=len(B),exactRetainedNodes=int(same.sum()),interpolatedNodes=len(index),
        first24LayersUnchanged=True,maximumCandidates=maximumCandidates,maximumLinearCoordinateError=float(coordinateError),
        minimumBarycentricWeight=float(w.min()),minimumDensity=float(rho.min()),minimumInternalEnergy=float(internal.min()),
        strictlyConservativeRemap=False,initialGuessOnly=True,physicsQualified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),source_mesh,target_mesh,restart]},outputSHA256=sha(output))
    output.with_suffix('.transfer.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-mesh',type=Path,required=True);p.add_argument('--target-mesh',type=Path,required=True)
    p.add_argument('--restart',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    transfer(a.source_mesh.resolve(),a.target_mesh.resolve(),a.restart.resolve(),a.output.resolve())
