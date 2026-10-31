"""独立检查近壁延长候选和旧单元插值；不调用迁移器，不重跑网格或求解器。"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v116_validation/study01'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'preflight_audit.json';assert not target.exists()
    meta=json.loads((OUT/'preparation.json').read_text(encoding='utf-8'))
    oldpath=ROOT/meta['oldMeshFolder']/'hybrid_volume.mat';newpath=OUT/'mesh32/hybrid_volume.mat'
    a=loadmat(str(oldpath));b=loadmat(str(newpath));A=a['vertices'];B=b['vertices']
    oldfile=ROOT/meta['oldEngineFolder']/'restart_flow_002999.csv'
    initial=OUT/'candidate/initial.csv';old=np.loadtxt(str(oldfile),delimiter=',',skiprows=1)
    new=np.loadtxt(str(initial),delimiter=',',skiprows=1)
    assert sha(OUT/'control/initial.csv')==sha(oldfile)
    assert (OUT/'control/case.cfg').read_bytes()==(OUT/'candidate/case.cfg').read_bytes()
    assert np.array_equal(a['sourceVertices'],b['sourceVertices'])
    assert np.array_equal(a['sourceFaces'],b['sourceFaces'])
    nv=len(a['sourceVertices']);assert np.array_equal(A[:25*nv],B[:25*nv])
    assert np.array_equal(new[:25*nv,4:],old[:25*nv,4:])
    assert np.max(abs(old[:,1:4]-A))<1e-14
    assert np.array_equal(new[:,0],np.arange(len(B))) and np.array_equal(new[:,1:4],B)
    # 按真实旧四面体，用行列式体积比重新定位，独立于迁移器的逆矩阵计算。
    lookup={tuple(p):i for i,p in enumerate(A)};unmatched=[]
    for i,p in enumerate(B):
        j=lookup.get(tuple(p))
        if j is None:unmatched.append(i)
        else:assert np.array_equal(new[i,4:],old[j,4:])
    T=a['tetrahedra'].astype(int)-1;P=A[T];d=P[:,1:]-P[:,0,None,:]
    determinant=np.einsum('ij,ij->i',d[:,0],np.cross(d[:,1],d[:,2]))
    tree=cKDTree(P.mean(axis=1));remaining=np.array(unmatched);error=0.;coord=0.;minimum=1.
    for k in (32,128,512,2048):
        missing=[]
        for start in range(0,len(remaining),64):
            ids=remaining[start:start+64];_,near=tree.query(B[ids],k=min(k,len(T)))
            delta=B[ids,None,:]-P[near,0];edges=d[near]
            w1=np.einsum('nki,nki->nk',delta,np.cross(edges[:,:,1],edges[:,:,2]))/determinant[near]
            w2=np.einsum('nki,nki->nk',edges[:,:,0],np.cross(delta,edges[:,:,2]))/determinant[near]
            w3=np.einsum('nki,nki->nk',edges[:,:,0],np.cross(edges[:,:,1],delta))/determinant[near]
            w=np.stack([1-w1-w2-w3,w1,w2,w3],axis=2)
            inside=(w.min(axis=2)>=-2e-11)&(w.max(axis=2)<=1+2e-11)
            hit=inside.any(axis=1);pick=inside.argmax(axis=1);good=ids[hit]
            weights=w[np.arange(len(ids)),pick][hit];cells=T[near[np.arange(len(ids)),pick][hit]]
            if len(good):
                expected=np.einsum('ni,nij->nj',weights,old[cells,4:])
                error=max(error,float(np.max(abs(expected-new[good,4:])/np.maximum(1,abs(new[good,4:])))))
                coord=max(coord,float(np.max(abs(np.einsum('ni,nij->nj',weights,A[cells])-B[good]))))
                minimum=min(minimum,float(weights.min()))
            missing.extend(ids[~hit])
        remaining=np.array(missing,dtype=int)
        if not len(remaining):break
    assert not len(remaining) and error<2e-11 and coord<1e-12
    # 完整体域的三角/四边形共享与方向，确保延长的层和四面体没有洞或叠面。
    W=b['wedges'].astype(int)-1;T=b['tetrahedra'].astype(int)-1;F=b['sourceFaces'].astype(int)-1
    far=b['farfield'].astype(int)-1
    triangles=np.concatenate([W[:,[0,2,1]],W[:,[3,4,5]],T[:,[1,2,3]],T[:,[0,3,2]],T[:,[0,1,3]],T[:,[0,2,1]]])
    quads=np.concatenate([W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]]])
    topology={}
    for name,faces,boundary in [('triangles',triangles,np.concatenate([F[:,[0,2,1]],far])),('quads',quads,np.zeros((0,4),int))]:
        keys,group,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
        least=faces.argmin(axis=1);rows=np.arange(len(faces));n=faces.shape[1]
        signs=np.sign(faces[rows,(least+1)%n]-faces[rows,(least-1)%n])
        direction=np.bincount(group,weights=signs)
        assert (count<=2).all() and (direction[count==2]==0).all()
        assert np.array_equal(keys[count==1],np.unique(np.sort(boundary,axis=1),axis=0))
        topology[name]=dict(internal=int(sum(count==2)),boundary=int(sum(count==1)))
    rho=new[:,4];e=new[:,8]/rho-.5*np.sum((new[:,5:8]/rho[:,None])**2,axis=1)
    assert np.isfinite(new).all() and min(rho)>0 and min(e)>0
    transfer=json.loads((OUT/'candidate/initial.transfer.json').read_text())
    for file,h in transfer['sourceSHA256'].items():assert sha(Path(file))==h
    assert sha(initial)==transfer['outputSHA256']
    receipt=json.loads((ROOT/'outputs/v116_validation/prepare01_process.json').read_text(encoding='utf-8'))
    assert receipt['status']=='terminal' and receipt['exitCode']==0
    paths=[Path(__file__).resolve(),oldpath,newpath,oldfile,initial,OUT/'control/case.cfg',OUT/'candidate/case.cfg',OUT/'candidate/initial.transfer.json']
    result=dict(exactRetainedNodes=len(B)-len(unmatched),interpolatedNodes=len(unmatched),
        first24LayersExact=True,normalizedStateDifference=error,linearCoordinateError=coord,
        minimumWeight=minimum,topology=topology,minimumInternalEnergy=float(min(e)),
        minimumTemperatureK=float(min(e)/(287.05287/.4)),
        stateWarning='Positive but near-zero inherited temperature; unconverged numerical diagnostic only.',
        physicsQualified=False,preflightPassed=True,
        SHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='SHA256'}))
if __name__=='__main__':main()
