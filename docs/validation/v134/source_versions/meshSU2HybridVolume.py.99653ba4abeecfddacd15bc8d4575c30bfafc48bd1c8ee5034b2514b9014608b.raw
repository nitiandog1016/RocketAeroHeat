"""原三角面法向棱柱层与Gmsh四面体外域；有限时子进程，不计算流场。"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import numpy as np
from scipy.io import loadmat,savemat

ROOT=Path(__file__).resolve().parents[1]
RUNTIME=ROOT/'tools/runtime/gmsh/4.15.2'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def wedge_jacobian_min(p):
    # 线性棱柱Jacobian对三角坐标为仿射、对高度为二次；检查三个角点的精确极小值。
    low=p[:,:3];delta=p[:,3:]-low;mins=[]
    for corner in range(3):
        values=[]
        for z in (0.,.5,1.):
            a=low[:,1]-low[:,0]+z*(delta[:,1]-delta[:,0])
            b=low[:,2]-low[:,0]+z*(delta[:,2]-delta[:,0])
            values.append(np.einsum('ij,ij->i',np.cross(a,b),delta[:,corner]))
        c=values[0];aa=2*(values[2]+c-2*values[1]);bb=values[2]-c-aa
        m=np.minimum(values[0],values[2]);active=aa>0
        z=np.zeros(len(m));z[active]=-bb[active]/(2*aa[active]);active&=(z>0)&(z<1)
        m[active]=np.minimum(m[active],aa[active]*z[active]**2+bb[active]*z[active]+c[active]);mins.append(m)
    return np.min(mins,axis=0)

def prepare_layers(V,F,opt):
    assert V.ndim==2 and V.shape[1]==3 and np.isfinite(V).all()
    assert F.ndim==2 and F.shape[1]==3 and F.min()>=0 and F.max()<len(V)
    E=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);keys=np.sort(E,axis=1)
    _,inv,counts=np.unique(keys,axis=0,return_inverse=True,return_counts=True)
    assert np.all(counts==2) and np.all(np.bincount(inv,weights=np.sign(E[:,1]-E[:,0]))==0),'Surface must be closed and consistently oriented'
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    A=coo_matrix((np.ones(len(E)),(E[:,0],E[:,1])),shape=(len(V),len(V))).tocsr()
    assert connected_components(A,directed=False,return_labels=False)==1,'One connected shell required'
    p=V[F];av=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])/2;area=np.linalg.norm(av,axis=1)
    assert np.all(area>0);normal=av/area[:,None]
    body=np.einsum('ij,ij->',p[:,0],av)/3;assert body>0,'Surface must point out of solid'
    N=np.zeros_like(V)
    for c in range(3):
        a=p[:,(c+1)%3]-p[:,c];b=p[:,(c+2)%3]-p[:,c]
        angles=np.arctan2(np.linalg.norm(np.cross(a,b),axis=1),np.einsum('ij,ij->i',a,b))
        np.add.at(N,F[:,c],angles[:,None]*normal)
    N/=np.linalg.norm(N,axis=1)[:,None]
    alignment=np.einsum('fcj,fj->fc',N[F],normal)
    assert np.isfinite(N).all() and np.all(alignment>0),'Normal points into adjacent surface; no geometry alteration allowed'
    h=np.asarray(opt['layerHeights'],float)
    assert len(h)>=2 and len(h)<=65 and h[0]==0 and np.all(np.diff(h)>0) and np.isfinite(h).all()
    layers=np.concatenate([V+d*N for d in h]);nv=len(V)
    wedges=np.concatenate([np.concatenate([F+j*nv,F+(j+1)*nv],axis=1) for j in range(len(h)-1)])
    jac=wedge_jacobian_min(layers[wedges]);assert np.all(jac>0),'Inverted prism; reduce geometry complexity with explicit new input, not hidden clipping'
    return layers,wedges,area,body,alignment,jac

def build(V,F,opt):
    sys.path.insert(0,str(RUNTIME));import gmsh
    assert gmsh.__version__=='4.15.2'
    layers,W,areas,body,alignment,jac=prepare_layers(V,F,opt)
    nv=len(V);nl=len(opt['layerHeights']);shell=layers[-nv:]
    L=float(np.ptp(V[:,0]));assert L>0
    lo=np.array([V[:,0].min()-opt['downstreamLengths']*L,-opt['lateralLengths']*L,-opt['lateralLengths']*L])
    hi=np.array([V[:,0].max()+opt['upstreamLengths']*L,opt['lateralLengths']*L,opt['lateralLengths']*L])
    assert np.all(shell>lo) and np.all(shell<hi),'Entire body and collar must be inside farfield'
    gmsh.initialize();gmsh.option.setNumber('General.NumThreads',2)
    try:
        gmsh.model.add('RocketAeroHeat-original-faces')
        box=gmsh.model.occ.addBox(*lo,*(hi-lo));gmsh.model.occ.synchronize()
        outside=[t for d,t in gmsh.model.getBoundary([(3,box)],oriented=False)]
        gmsh.model.occ.remove([(3,box)],recursive=False);gmsh.model.occ.synchronize()
        # 先创建CAD外域，再分配不冲突离散面tag，避免两种内核共享tag覆盖壳面。
        inner=gmsh.model.addDiscreteEntity(2,max(outside)+1)
        gmsh.model.mesh.addNodes(2,inner,np.arange(1,nv+1),shell.ravel())
        gmsh.model.mesh.addElementsByType(inner,2,np.arange(1,len(F)+1),(F+1).ravel())
        outloop=gmsh.model.geo.addSurfaceLoop(outside);inloop=gmsh.model.geo.addSurfaceLoop([inner])
        region=gmsh.model.geo.addVolume([outloop,inloop]);gmsh.model.geo.synchronize()
        gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1);gmsh.option.setNumber('Mesh.Algorithm3D',1)
        gmsh.option.setNumber('Mesh.Renumber',0)
        gmsh.option.setNumber('Mesh.MeshSizeMin',opt['outerMinSize'])
        gmsh.option.setNumber('Mesh.MeshSizeMax',opt['farSize'])
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',1)
        gmsh.option.setNumber('Mesh.Optimize',1);gmsh.option.setNumber('Mesh.SaveAll',1)
        gmsh.model.mesh.generate(3)
        types,tags,nodeLists=gmsh.model.mesh.getElements(2,inner)
        assert list(types)==[2] and len(tags[0])==len(F)
        # 外壳作为空腔边界可被引擎反向或重排；只允许同一组原三角连接，不允许细分/换对角线。
        assert np.array_equal(np.unique(np.sort(nodeLists[0].reshape(-1,3),axis=1),axis=0),
                              np.unique(np.sort(F+1,axis=1),axis=0)),'Interface triangles changed'
        nodeTags,coordinates,_=gmsh.model.mesh.getNodes();coordinates=coordinates.reshape(-1,3)
        order=np.argsort(nodeTags);nodeTags=nodeTags[order];coordinates=coordinates[order]
        original=np.searchsorted(nodeTags,np.arange(1,nv+1))
        assert np.array_equal(coordinates[original],shell),'Interface nodes moved'
        newMask=~np.isin(nodeTags,np.arange(1,nv+1));newTags=nodeTags[newMask]
        P=np.concatenate([layers,coordinates[newMask]])
        mapping=np.full(int(nodeTags.max())+1,-1,dtype=np.int64)
        mapping[1:nv+1]=np.arange((nl-1)*nv,nl*nv)
        mapping[newTags]=np.arange(len(layers),len(P))
        types,tags,nodeLists=gmsh.model.mesh.getElements(3,region);assert list(types)==[4],'Only linear tetrahedra supported in outer domain'
        quality=gmsh.model.mesh.getElementQualities(tags[0],'minSICN')
        T=mapping[nodeLists[0].reshape(-1,4)];assert np.all(T>=0)
        far=[]
        for face in outside:
            types,_,ids=gmsh.model.mesh.getElements(2,face);assert list(types)==[2]
            f=mapping[ids[0].reshape(-1,3)];q=P[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
            flip=np.einsum('ij,ij->i',q.mean(axis=1)-(lo+hi)/2,n)<0
            f[flip]=f[flip][:,[0,2,1]];far.append(f)
        far=np.concatenate(far)
        q=P[T];det=np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))
        assert np.all(det>0),'Nonpositive Gmsh tetrahedron'
        assert len(P)<=250000 and len(T)+len(W)<=600000,'Geometry-only research cell budget exceeded'
        report=dict(gmshVersion=gmsh.__version__,originalFaces=len(F),originalVertices=nv,
                    nodes=len(P),prisms=len(W),tetrahedra=len(T),outerBoundaryFaces=len(far),
                    minimumPrismJacobian=float(jac.min()),minimumTetVolume=float(det.min()/6),
                    minimumTetSICN=float(np.min(quality)),tetsBelowSICN01=int(np.sum(quality<.1)),
                    minimumNormalAlignment=float(alignment.min()),bodyArea=float(areas.sum()),bodyVolume=float(body),
                    farfieldBox=[lo.tolist(),hi.tolist()],farfieldVolume=float(np.prod(hi-lo)),
                    originalFacesUnchanged=True,interfaceExactlyPreserved=True,physicsQualified=False,
                    newQualifiedArea=0,nearWallResolutionValidated=False,farfieldAdequacyValidated=False)
        return dict(vertices=P,wedges=W+1,tetrahedra=T+1,farfield=far+1,sourceVertices=V,sourceFaces=F+1),report
    finally:gmsh.finalize()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--folder',required=True)
    ap.add_argument('--worker',action='store_true');ap.add_argument('--maximum-seconds',type=float,default=720)
    a=ap.parse_args();src=Path(a.input).resolve();folder=Path(a.folder).resolve()
    if a.worker:
        data=loadmat(str(src));V=data['V'].astype(float);F=data['F'].astype(np.int64)-1
        opt=json.loads(str(data['optionsJSON'][0]));result,report=build(V,F,opt)
        savemat(str(folder/'hybrid_volume.mat'),result,do_compression=True)
        (folder/'mesh_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');return
    assert 0<a.maximum_seconds<=720 and not folder.exists();folder.mkdir(parents=True)
    assert (RUNTIME/'runtime.json').is_file(),'Run scripts/installGmshRuntime.py first'
    runtime=json.loads((RUNTIME/'runtime.json').read_text())
    for entry in runtime['files']:assert sha(RUNTIME/entry['path'])==entry['sha256']
    record=dict(status='running',input=str(src),inputSHA256=sha(src),scriptSHA256=sha(__file__),
                runtime=runtime,maximumSeconds=a.maximum_seconds,threads=2,supervisorPID=os.getpid())
    receipt=folder/'mesh_process.json'
    def save():receipt.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    start=time.monotonic();save();timed=False
    with (folder/'gmsh.log').open('wb') as log:
        process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker','--input',str(src),'--folder',str(folder)],
             stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0,
             env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
        record['processID']=process.pid;save()
        try:
            try:code=process.wait(timeout=a.maximum_seconds)
            except subprocess.TimeoutExpired:timed=True;process.kill();code=process.wait()
        finally:
            if process.poll() is None:process.kill();process.wait()
    record.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start,
                  inputStable=sha(src)==record['inputSHA256'],sourceStable=sha(__file__)==record['scriptSHA256'])
    record['outputSHA256']={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p!=receipt};save()
    print(json.dumps({k:record[k] for k in ('status','exitCode','timedOut','elapsedSeconds')}))
    raise SystemExit(124 if timed else code)

if __name__=='__main__':main()
