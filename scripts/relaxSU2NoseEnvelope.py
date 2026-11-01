"""Gmsh球面重定位仅修整体网格外层，原壁面不动；先保存失败候选再作门槛判断。"""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
import numpy as np
from scipy.io import loadmat,savemat
from meshSU2HybridVolume import wedge_jacobian_min,RUNTIME
from buildSU2NoseEnvelope import edge_stencils
from prepareSU2CheckpointContinuation import ROOT,read,sha

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')

def build(source,folder):
    old=loadmat(str(source/'hybrid_volume.mat'));report=read(source/'mesh_report.json');proc=read(source/'mesh_process.json')
    assert proc['exitCode']==0 and proc['inputStable'] and proc['sourceStable']
    for name in ['hybrid_volume.mat','mesh_report.json']:assert sha(source/name)==proc['outputSHA256'][name]
    V=old['sourceVertices'];F=old['sourceFaces'].astype(int)-1;nv=len(V);P=old['vertices'].copy()
    h=np.array(report['normalHeightsM']);nl=len(h)-1;S=P[nl*nv:(nl+1)*nv].copy();mask=V[:,0]>.03
    A=np.column_stack([2*V[mask],np.ones(sum(mask))]);fit=np.linalg.lstsq(A,np.sum(V[mask]**2,axis=1),rcond=None)[0]
    c=fit[:3];R=np.sqrt(fit[3]+c@c)+h[-1];target=c+R*(V-c)/np.linalg.norm(V-c,axis=1)[:,None]
    u=np.clip((V[:,0]-.03)/.01,0,1);weight=u*u*(3-2*u);shell=S+weight[:,None]*(target-S)
    cut=float(min(V[V[:,0]>=.04,0]));rim=np.flatnonzero(V[:,0]==cut);cap=np.flatnonzero(V[:,0]>=cut)
    apex=int(np.argmax(V[:,0]));interior=np.setdiff1d(cap,np.r_[rim,apex]);cutX=float(np.mean(target[rim,0]))
    circleR=np.sqrt(R**2-(cutX-c[0])**2);yz=target[rim,1:]-c[1:]
    shell[rim,0]=cutX;shell[rim,1:]=c[1:]+circleR*yz/np.linalg.norm(yz,axis=1)[:,None];shell[apex]=c+[R,0,0]
    faceMask=np.all(V[F,0]>=cut,axis=1);faces=F[faceMask]
    edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    _,keep,count=np.unique(np.sort(edges,axis=1),axis=0,return_index=True,return_counts=True)
    boundary=edges[keep[count==1]];assert np.array_equal(np.unique(boundary),rim)
    detail=dict(method='gmsh-spherical-cad-relocate2d',optimizationStarted=False,optimizationCompleted=False,
        capNodes=len(cap),capFaces=len(faces),cutSourceXM=cut,cutShellXM=cutX,sphereCenterM=c.tolist(),
        sphereOuterRadiusM=float(R),displacementBoundM=.002,optimizer='Relocate2D',iterations=5,physicsQualified=False,newQualifiedArea=0)
    write(folder/'surface_report.json',detail)
    sys.path.insert(0,str(RUNTIME));import gmsh
    assert gmsh.__version__=='4.15.2';gmsh.initialize();gmsh.option.setNumber('General.NumThreads',2)
    try:
        gmsh.model.add('RocketAeroHeat-nose-envelope')
        vol=gmsh.model.occ.addSphere(*c,R,angle1=float(np.arcsin((cutX-c[0])/R)),angle2=np.pi/2)
        gmsh.model.occ.rotate([(3,vol)],*c,0,1,0,np.pi/2);gmsh.model.occ.synchronize()
        spheres=[t for dim,t in gmsh.model.getEntities(2) if gmsh.model.getType(dim,t)=='Sphere'];assert len(spheres)==1
        sphere=spheres[0];curves=[t for dim,t in gmsh.model.getBoundary([(2,sphere)],oriented=False) if dim==1]
        circles=[t for t in curves if gmsh.model.getType(1,t)=='Circle'];assert len(circles)==1;curve=circles[0]
        vertex=min(gmsh.model.getEntities(0),key=lambda dt:np.linalg.norm(np.array(gmsh.model.getValue(*dt,[]))-shell[apex]))
        assert np.linalg.norm(np.array(gmsh.model.getValue(*vertex,[]))-shell[apex])<1e-10
        for dim,tag,ids in [(2,sphere,interior),(1,curve,rim),(0,vertex[1],np.array([apex]))]:
            xyz=shell[ids].ravel();param=gmsh.model.getParametrization(dim,tag,xyz) if dim else []
            gmsh.model.mesh.addNodes(dim,tag,ids+1,xyz,param)
        gmsh.model.mesh.addElementsByType(curve,1,[],(boundary+1).ravel())
        gmsh.model.mesh.addElementsByType(sphere,2,np.flatnonzero(faceMask)+1,(faces+1).ravel())
        gmsh.option.setNumber('Mesh.Renumber',0);detail['optimizationStarted']=True;write(folder/'surface_report.json',detail)
        gmsh.model.mesh.optimize('Relocate2D',force=True,niter=5,dimTags=[(2,sphere)])
        tags,xyz,_=gmsh.model.mesh.getNodes();xyz=xyz.reshape(-1,3)
        assert len(tags)==len(cap) and np.array_equal(np.sort(tags),cap+1)
        result=shell.copy();result[tags.astype(int)-1]=xyz
        types,_,conn=gmsh.model.mesh.getElements(2,sphere);assert list(types)==[2]
        assert np.array_equal(np.unique(np.sort(conn[0].reshape(-1,3),axis=1),axis=0),np.unique(np.sort(faces+1,axis=1),axis=0))
        detail['optimizationCompleted']=True
        savemat(str(folder/'surface_raw.mat'),dict(oldShell=S,newShell=result,sourceVertices=V,sourceFaces=F+1),do_compression=True)
        assert np.array_equal(result[rim],shell[rim]) and np.array_equal(result[apex],shell[apex])
        st=edge_stencils(F);st=st[np.any(mask[st],axis=1)];q=result[st]
        normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        distance=np.sum((q[:,3]-q[:,0])*normal,axis=1)
        detail.update(maximumDisplacementM=float(max(np.linalg.norm(result-S,axis=1))),
            concaveEdges=int(sum(distance>1e-12)),maximumConcavityM=float(max(distance)))
        write(folder/'surface_report.json',detail)
        assert detail['maximumDisplacementM']<=.002+1e-12,'外壳位移超2mm，保留原始结果'
        assert np.all(normal[:,0]>0) and max(distance)<=1e-10,'保留连接仍有凹折，不启动外域重建'
    finally:gmsh.finalize()
    delta=result-S
    for j in range(nl+1):
        u=max(0.,(h[j]-.001)/(h[-1]-.001));w=u*u*(3-2*u)
        if w>0:P[j*nv:(j+1)*nv]+=w*delta
    jac=wedge_jacobian_min(P[old['wedges'].astype(int)-1]);detail['minimumPrismJacobian']=float(min(jac))
    write(folder/'surface_report.json',detail);assert min(jac)>0
    assert np.array_equal(P[:sum(h<=.001)*nv],old['vertices'][:sum(h<=.001)*nv])
    savemat(str(folder/'collar_prepared.mat'),dict(vertices=P,wedges=old['wedges'],farfield=old['farfield'],
        sourceVertices=V,sourceFaces=F+1,oldOuterShell=S,newOuterShell=result),do_compression=True)
    report.update(noseEnvelope=detail,noseEnvelopePrepared=True,exteriorRequiresRebuild=True,physicsQualified=False,newQualifiedArea=0)
    write(folder/'collar_report.json',report);savemat(str(folder/'mesh_input.mat'),loadmat(str(source/'mesh_input.mat')),do_compression=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--target',type=Path,required=True)
    p.add_argument('--worker',action='store_true');a=p.parse_args();source=a.source.resolve();folder=a.target.resolve()
    source.relative_to(ROOT/'outputs');folder.relative_to(ROOT/'outputs')
    if a.worker:build(source,folder);return
    assert not folder.exists();folder.mkdir(parents=True)
    files=[source/n for n in ['mesh_process.json','mesh_report.json','hybrid_volume.mat','mesh_input.mat']]+[
        Path(__file__).resolve(),ROOT/'scripts/buildSU2NoseEnvelope.py',ROOT/'scripts/meshSU2HybridVolume.py',
        ROOT/'docs/validation/v134/MATURE_SURFACE_PLAN.md']+[RUNTIME/n for n in ['runtime.json','gmsh.py','gmsh-4.15.dll']]
    locks={str(p):sha(p) for p in files};start=time.monotonic();receipt=folder/'prepare_process.json'
    record=dict(status='running',maximumSeconds=180,sourceHashes=locks,sourceMeshDirectory=str(source),nativeFlowRuns=0,supervisorPID=os.getpid())
    write(receipt,record)
    with (folder/'prepare.log').open('wb') as log:
        proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker','--source',str(source),'--target',str(folder)],
            stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,
            env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',PYTHONIOENCODING='utf-8'))
        record['processID']=proc.pid;write(receipt,record);timed=False
        try:code=proc.wait(timeout=180)
        except subprocess.TimeoutExpired:timed=True;proc.kill();code=proc.wait()
    stable={p:sha(Path(p))==value for p,value in locks.items()}
    record.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start,
        sourceStableAfter=stable,sourceStable=all(stable.values()),inputStable=all(stable.values()),
        outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p!=receipt})
    write(receipt,record);print(json.dumps({k:record[k] for k in ['status','exitCode','elapsedSeconds']}))
    raise SystemExit(124 if timed else code)

if __name__=='__main__':main()
