"""保留近壁棱柱，Gmsh重填完整外域并自动连接四边侧面；不运行流场。"""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.spatial import ConvexHull
from meshSU2HybridVolume import prepare_layers,wedge_jacobian_min,RUNTIME
from auditSU2HybridTopology import audit,polyvolume
from optimizeSU2FanExterior import sicn
from prepareSU2CheckpointContinuation import ROOT,read,sha

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')

def oriented_hull(V):
    F=ConvexHull(V).simplices.copy();p=V[F]
    flip=np.sum(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])*(p.mean(axis=1)-V.mean(axis=0)),axis=1)<0
    F[flip]=F[flip][:,[0,2,1]];return F

def boundary(W,F):
    """从保留棱柱的实际共享面导出空腔，侧四边面保持完整，禁止三角化偷换。"""
    tri=np.concatenate([W[:,[0,2,1]],W[:,[3,4,5]]]);quad=np.concatenate([W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]]])
    result=[]
    for faces in [tri,quad]:
        _,index,count=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_counts=True)
        assert np.all(count<=2);result.append(faces[index[count==1]])
    tri,quad=result;wall={tuple(v) for v in np.sort(F,axis=1)}
    keep=np.array([tuple(v) not in wall for v in np.sort(tri,axis=1)])
    assert sum(~keep)==len(F)
    return tri[keep][:,[0,2,1]],quad[:,[0,3,2,1]]

def prepare(source,fixture):
    if fixture:
        V=np.array([[1.,0,0],[-1.,0,0],[0,1.,0],[0,-1.,0],[0,0,1.],[0,0,-1.]])
        F=oriented_hull(V);h=np.array([0.,.1,.3,.6]);P,W,_,_,_,_=prepare_layers(V,F,dict(layerHeights=h))
        p=V[F];cap=np.mean(p[:,:,0],axis=1)>0;retained=1
        farV=np.array([[x,y,z] for x in [-2.,2.] for y in [-2.,2.] for z in [-2.,2.]]);farF=oriented_hull(farV)
        far=farF+len(P);P=np.concatenate([P,farV]);report=dict(normalHeightsM=h.tolist(),normalLayers=3)
        sourcePrisms=len(W);keep=~np.concatenate([cap & (j>=retained) for j in range(3)])
        opts=dict(outerMinSize=.15,farSize=.65,layerHeights=h.tolist())
    else:
        r=read(source/'mesh_process.json');assert r['exitCode']==0 and r['sourceStable'] and r['inputStable']
        for name in ['hybrid_volume.mat','mesh_report.json']:assert sha(source/name)==r['outputSHA256'][name]
        d=loadmat(str(source/'hybrid_volume.mat'));report=read(source/'mesh_report.json')
        V=d['sourceVertices'];F=d['sourceFaces'].astype(int)-1;P=d['vertices'];W=d['wedges'].astype(int)-1;far=d['farfield'].astype(int)-1
        h=np.array(report['normalHeightsM']);retained=int(np.flatnonzero(h<=.001)[-1]);assert retained==19
        cut=float(min(V[V[:,0]>=.03,0]));cap=np.all(V[F,0]>=cut,axis=1)
        assert np.any(cap) and len(V)==9878 and len(F)==19752
        sourcePrisms=len(W);keep=np.r_[~np.concatenate([cap & (j>=retained) for j in range(len(h)-1)]),np.ones(len(W)-(len(h)-1)*len(F),bool)]
        opts=json.loads(str(loadmat(str(source/'mesh_input.mat'))['optionsJSON'][0]))
        report['noseCutSourceXM']=cut;report['sourcePrismsBeforeTermination']=len(W)
    W=W[keep];tri,quad=boundary(W,F);assert len(quad)>0
    report.update(noseRetainedLayers=retained,noseRetainedHeightM=float(h[retained]),noseCapFaces=int(sum(cap)),
        terminatedPrisms=int(sourcePrisms-len(W)),transitionQuadrilateralFaces=len(quad),protectedWallAndNearWallUnchanged=True)
    return P,W,V,F,far,tri,quad,report,opts,keep

def build(source,folder,fixture):
    P,W,V,F,far,tri,quad,report,opts,keep=prepare(source,fixture)
    assert min(wedge_jacobian_min(P[W]))>0
    ids=np.unique(np.r_[tri.ravel(),quad.ravel(),far.ravel()]);innerIDs=np.unique(np.r_[tri.ravel(),quad.ravel()]);farIDs=np.unique(far)
    assert not np.intersect1d(innerIDs,farIDs).size
    sys.path.insert(0,str(RUNTIME));import gmsh
    assert gmsh.__version__=='4.15.2';gmsh.initialize();gmsh.option.setNumber('General.NumThreads',2)
    try:
        gmsh.model.add('RocketAeroHeat-prism-termination')
        for tag,nodes in [(1,innerIDs),(2,farIDs)]:
            gmsh.model.addDiscreteEntity(2,tag);gmsh.model.mesh.addNodes(2,tag,nodes+1,P[nodes].ravel())
        gmsh.model.mesh.addElementsByType(1,2,[],(tri+1).ravel());gmsh.model.mesh.addElementsByType(1,3,[],(quad+1).ravel())
        gmsh.model.mesh.addElementsByType(2,2,[],(far+1).ravel())
        inner=gmsh.model.geo.addSurfaceLoop([1]);outer=gmsh.model.geo.addSurfaceLoop([2])
        vol=gmsh.model.geo.addVolume([outer,inner]);gmsh.model.geo.synchronize()
        for name,value in [('Mesh.MeshOnlyEmpty',1),('Mesh.Algorithm3D',1),('Mesh.Renumber',0),('Mesh.Optimize',1),
            ('Mesh.MeshSizeMin',opts['outerMinSize']),('Mesh.MeshSizeMax',opts['farSize']),('Mesh.MeshSizeExtendFromBoundary',1)]:gmsh.option.setNumber(name,value)
        gmsh.model.mesh.generate(3)
        tags,xyz,_=gmsh.model.mesh.getNodes();xyz=xyz.reshape(-1,3)
        types,et,conn=gmsh.model.mesh.getElements(3,vol)
        raw=dict(nodeTags=tags,vertices=xyz,elementTypes=np.array(types))
        for typ,nodes in zip(types,conn):raw['elements'+str(typ)]=nodes
        savemat(str(folder/'gmsh_raw.mat'),raw,do_compression=True)
        assert set(types)<=set([4,7]) and 7 in types,'需要真实四面体/金字塔连接'
        position={int(tag):i for i,tag in enumerate(tags)}
        assert np.array_equal(xyz[[position[int(i+1)] for i in ids]],P[ids]),'边界节点变化'
        for entity,typ,original in [(1,2,tri),(1,3,quad),(2,2,far)]:
            types2,_,conn2=gmsh.model.mesh.getElements(2,entity);actual=conn2[list(types2).index(typ)].reshape(-1,original.shape[1])-1
            assert np.array_equal(np.unique(np.sort(actual,axis=1),axis=0),np.unique(np.sort(original,axis=1),axis=0))
        newMask=~np.isin(tags,ids+1);newP=np.concatenate([P,xyz[newMask]])
        mapping=np.full(int(max(tags))+1,-1,int);mapping[ids+1]=ids;mapping[tags[newMask]]=np.arange(len(P),len(newP))
        T=np.empty((0,4),int);Q=np.empty((0,5),int)
        for typ,nodes in zip(types,conn):
            if typ==4:T=mapping[nodes.reshape(-1,4)]
            else:Q=mapping[nodes.reshape(-1,5)]
        assert min(T.ravel())>=0 and min(Q.ravel())>=0
    finally:gmsh.finalize()
    used=np.unique(np.r_[W.ravel(),T.ravel(),Q.ravel(),far.ravel()]);remap=np.full(len(newP),-1,int);remap[used]=np.arange(len(used))
    assert np.array_equal(used[:len(V)],np.arange(len(V)))
    result=dict(vertices=newP[used],wedges=remap[W]+1,tetrahedra=remap[T]+1,pyramids=remap[Q]+1,
        farfield=remap[far]+1,sourceVertices=V,sourceFaces=F+1)
    savemat(str(folder/'hybrid_volume.mat'),result,do_compression=True)
    auditResult=audit(result);write(folder/'geometry_audit.json',auditResult)
    quality=sicn(result['vertices'],result['tetrahedra']-1)
    areas=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    report.update(nodes=len(used),prisms=len(W),pyramids=len(Q),tetrahedra=len(T),gmshVersion='4.15.2',
        minimumTetSICN=float(min(quality)),tetsBelowSICN01=int(sum(quality<.1)),minimumTetVolume=auditResult['minimumTetVolume'],
        minimumPrismJacobian=float(min(wedge_jacobian_min(result['vertices'][result['wedges']-1]))),
        bodyArea=float(sum(areas)),bodyVolume=polyvolume(V,F),farfieldVolume=polyvolume(result['vertices'],result['farfield']-1),
        originalFaces=len(F),originalVertices=len(V),outerBoundaryFaces=len(far),physicsQualified=False,newQualifiedArea=0,
        constructionOperation='closed-convex-rim-prism-termination-gmsh-pyramid-transition' if not fixture else 'octahedron-prism-termination-gmsh-pyramid-transition',
        originalFacesUnchanged=True,interfaceExactlyPreserved=True,exteriorTetrahedraExactlyPreserved=False,
        originalWallNodesExactlyPreserved=True,retainedPrismCoordinatesExactlyPreserved=True,
        nearWallResolutionValidated=False,farfieldAdequacyValidated=False,qualityGatePassed=bool(min(quality)>=.1))
    if not fixture:
        originalRegular=report['originalPrisms'];report['originalPrisms']=int(sum(keep[:originalRegular]))
        for name in ['sideNodeIDs','baseNodeIDs']:report[name]=remap[np.array(report[name])].tolist()
        assert min(np.array(report['sideNodeIDs']).ravel())>=0 and min(np.array(report['baseNodeIDs']).ravel())>=0
    write(folder/'mesh_report.json',report)
    savemat(str(folder/'mesh_input.mat'),dict(V=V,F=F+1,optionsJSON=json.dumps(opts)),do_compression=True)
    assert len(used)<=(20000 if fixture else 400000) and len(W)+len(T)+len(Q)<=(100000 if fixture else 900000)
    if not fixture:assert report['qualityGatePassed'],'实际外域SICN未过，不采纳'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path);parser.add_argument('--target',type=Path,required=True)
    parser.add_argument('--fixture',action='store_true');parser.add_argument('--worker',action='store_true');a=parser.parse_args()
    source=a.source.resolve() if a.source else None;folder=a.target.resolve();folder.relative_to(ROOT/'outputs')
    assert a.fixture or source is not None
    if a.worker:build(source,folder,a.fixture);return
    assert not folder.exists();folder.mkdir(parents=True);maximum=60 if a.fixture else 180
    files=[Path(__file__).resolve(),ROOT/'scripts/su2PyramidGeometry.py',ROOT/'scripts/auditSU2HybridTopology.py',ROOT/'scripts/meshSU2HybridVolume.py',
        ROOT/'scripts/optimizeSU2FanExterior.py',ROOT/'docs/validation/v135/PLAN.md']+[RUNTIME/n for n in ['runtime.json','gmsh.py','gmsh-4.15.dll']]
    if source:files += [source/n for n in ['hybrid_volume.mat','mesh_input.mat','mesh_report.json','mesh_process.json']]
    locks={str(p):sha(p) for p in files};record=dict(status='running',maximumSeconds=maximum,sourceHashes=locks,fixture=a.fixture,supervisorPID=os.getpid())
    receipt=folder/'mesh_process.json';write(receipt,record);start=time.monotonic()
    command=[sys.executable,str(Path(__file__).resolve()),'--worker','--target',str(folder)]
    command+=['--fixture'] if a.fixture else ['--source',str(source)]
    with (folder/'gmsh.log').open('wb') as log:
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,
            env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',PYTHONIOENCODING='utf-8'))
        record['processID']=p.pid;write(receipt,record);timed=False
        try:code=p.wait(timeout=maximum)
        except subprocess.TimeoutExpired:timed=True;p.kill();code=p.wait()
    stable={path:sha(Path(path))==h for path,h in locks.items()}
    record.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start,
        sourceStableAfter=stable,sourceStable=all(stable.values()),inputStable=all(stable.values()),
        input=str(folder/'mesh_input.mat'),inputSHA256=sha(folder/'mesh_input.mat') if (folder/'mesh_input.mat').exists() else None,
        outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p!=receipt})
    write(receipt,record);print(json.dumps({k:record[k] for k in ['status','exitCode','elapsedSeconds']}))
    raise SystemExit(124 if timed else code)

if __name__=='__main__':main()
