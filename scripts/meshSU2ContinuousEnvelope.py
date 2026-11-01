"""v141单次连续鼻部外包络；完整壳审核后，限时调用一次固定Gmsh填外域。"""
from pathlib import Path
import json,time,traceback,sys,argparse
import numpy as np
from scipy.io import loadmat,savemat
from redistributePrismEnvelope import redistribute
from terminatePrismLayers import graph_layers
from meshSU2HybridVolume import wedge_jacobian_min,RUNTIME
from auditSU2HybridTopology import audit
from reviewSU2ExplicitQuality import quality
from mmgLocalInterface import sorted_cells
from runCGALSurfaceAudit import audit_surface
from buildCGALRemeshingRuntime import run
from prepareCGALRemeshingRuntime import ROOT,sha

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def exterior(folder):
    """独立限时子进程：一次Delaunay，先保存原始结果，再做完整网格门槛。"""
    d=loadmat(str(folder/'prepared_envelope.mat'));P=d['vertices'];W=d['wedges'].astype(np.int64)-1
    top=d['farfield'].astype(np.int64)-1;far=d['exteriorFarfield'].astype(np.int64)-1
    opts=json.loads(str(d['optionsJSON'][0]));ids=np.unique(np.r_[top.ravel(),far.ravel()]);sys.path.insert(0,str(RUNTIME));import gmsh
    assert gmsh.__version__=='4.15.2';gmsh.initialize();gmsh.option.setNumber('General.NumThreads',2)
    try:
        gmsh.model.add('RocketAeroHeat-continuous-envelope')
        for tag,faces in [(1,top[:,[0,2,1]]),(2,far)]:
            ns=np.unique(faces);gmsh.model.addDiscreteEntity(2,tag);gmsh.model.mesh.addNodes(2,tag,ns+1,P[ns].ravel())
            gmsh.model.mesh.addElementsByType(tag,2,[],(faces+1).ravel())
        a=gmsh.model.geo.addSurfaceLoop([1]);b=gmsh.model.geo.addSurfaceLoop([2]);vol=gmsh.model.geo.addVolume([b,a]);gmsh.model.geo.synchronize()
        for k,v in [('Mesh.MeshOnlyEmpty',1),('Mesh.Algorithm3D',1),('Mesh.Renumber',0),('Mesh.Optimize',1),
            ('Mesh.MeshSizeMin',opts['outerMinSize']),('Mesh.MeshSizeMax',opts['farSize']),('Mesh.MeshSizeExtendFromBoundary',1)]:gmsh.option.setNumber(k,v)
        gmsh.model.mesh.generate(3);tags,xyz,_=gmsh.model.mesh.getNodes();xyz=xyz.reshape(-1,3);types,_,conn=gmsh.model.mesh.getElements(3,vol)
        raw=dict(nodeTags=tags,vertices=xyz,elementTypes=types)
        for typ,nodes in zip(types,conn):raw['elements'+str(typ)]=nodes
        savemat(str(folder/'gmsh_raw.mat'),raw,do_compression=True);assert list(types)==[4]
        position={int(v):i for i,v in enumerate(tags)};assert np.array_equal(xyz[[position[int(i+1)] for i in ids]],P[ids])
        for tag,expected in [(1,top),(2,far)]:
            typ,_,cells=gmsh.model.mesh.getElements(2,tag);assert list(typ)==[2]
            assert np.array_equal(sorted_cells(cells[0].reshape(-1,3)-1),sorted_cells(expected))
        add=~np.isin(tags,ids+1);allP=np.r_[P,xyz[add]];mapping=np.full(int(max(tags))+1,-1,np.int64)
        mapping[ids+1]=ids;mapping[tags[add]]=np.arange(len(P),len(allP));T=mapping[conn[0].reshape(-1,4)];assert T.min()>=0
    finally:gmsh.finalize()
    used=np.unique(np.r_[W.ravel(),T.ravel(),far.ravel()]);remap=np.full(len(allP),-1,np.int64);remap[used]=np.arange(len(used))
    data=dict(vertices=allP[used],wedges=remap[W]+1,tetrahedra=remap[T]+1,farfield=remap[far]+1,
        sourceVertices=d['sourceVertices'],sourceFaces=d['sourceFaces'],normalHeightsByColumnM=d['normalHeightsByColumnM'],
        sourcePreparedNodeIDsOneBased=used+1)
    assert np.array_equal(data['vertices'][data['wedges']-1],P[W]) and np.array_equal(data['vertices'][data['farfield']-1],P[far])
    savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True)
    geom=audit(data);write(folder/'geometry_audit.json',geom)
    q=quality(data['vertices'][data['tetrahedra']-1]);bad=np.flatnonzero(q<.1)
    report=dict(nodes=len(used),prisms=len(W),tetrahedra=len(T),totalElements=len(W)+len(T),minimumTetSICN=float(min(q)),
        tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=q[bad].tolist(),
        badTetCentersM=data['vertices'][data['tetrahedra'][bad]-1].mean(axis=1).tolist(),
        originalFaces=len(d['sourceFaces']),protectedNearWallExactlyPreserved=True,retainedTailFanExactlyPreserved=True,
        newPrismEnvelopeExactlyPreserved=True,originalFarfieldExactlyPreserved=True,normalHeightMode='per-column-continuous-envelope',
        geometryPassed=True,qualityGatePassed=bool(min(q)>=.1),sizeGatePassed=bool(len(used)<=400000 and len(W)+len(T)<=900000),
        exteriorAlgorithm='Gmsh4.15.2 Delaunay, fixed original farfield and continuous prism envelope',physicsQualified=False,newQualifiedArea=0)
    write(folder/'candidate_report.json',report)
    assert report['qualityGatePassed'] and report['sizeGatePassed'],'full exterior quality or size gate failed'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--exterior-worker',type=Path);args=parser.parse_args()
    if args.exterior_worker:exterior(args.exterior_worker.resolve());return
    start=time.monotonic();folder=ROOT/'outputs/v141_validation/actual';folder.mkdir(parents=True,exist_ok=False)
    source=ROOT/'outputs/v133_validation/fan_workbench';small=ROOT/'docs/validation/v141/fixture/report.json'
    paths=[small,ROOT/'docs/validation/v141/plan.md',Path(__file__).resolve()]+[source/n for n in ['mesh_process.json','mesh_report.json','mesh_input.mat','hybrid_volume.mat']]
    paths += [ROOT/'scripts'/n for n in ['redistributePrismEnvelope.py','terminatePrismLayers.py','meshSU2HybridVolume.py',
        'auditSU2HybridTopology.py','su2PyramidGeometry.py','reviewSU2ExplicitQuality.py','mmgLocalInterface.py',
        'runCGALSurfaceAudit.py','cgalSurfaceIntersections.py','buildCGALRemeshingRuntime.py']]
    paths += [RUNTIME/n for n in ['gmsh.py','gmsh-4.15.dll','runtime.json']]
    locked={str(p):sha(p) for p in paths};error=None;code=1;gmshStarted=False
    try:
        assert json.loads(small.read_text())['passed'];sr=json.loads((source/'mesh_process.json').read_text())
        assert sr['exitCode']==0 and sr['sourceStable'] and sr['inputStable']
        for n in ['mesh_report.json','hybrid_volume.mat']:assert sha(source/n)==sr['outputSHA256'][n]
        old=loadmat(str(source/'hybrid_volume.mat'));prior=json.loads((source/'mesh_report.json').read_text())
        P=old['vertices'];V=old['sourceVertices'];F=old['sourceFaces'].astype(np.int64)-1;W=old['wedges'].astype(np.int64)-1
        far=old['farfield'].astype(np.int64)-1;h=np.array(prior['normalHeightsM']);assert len(h)==33 and len(F)==19752 and len(V)==9878
        cut=.031232520723574634;_,distance=graph_layers(F,len(V),V[:,0]>=cut,19,32);join=float(max(V[distance>=13,0]))
        Q,H,r=redistribute(P,V,h,19,join,cut);top=F+32*len(V);tail=W[32*len(F):]
        assert np.array_equal(Q[W[:19*len(F)]],P[W[:19*len(F)]]) and np.array_equal(Q[tail],P[tail]) and np.array_equal(Q[far],P[far])
        opts=str(loadmat(str(source/'mesh_input.mat'))['optionsJSON'][0])
        data=dict(vertices=Q,wedges=W+1,tetrahedra=np.empty((0,4),np.int64),farfield=top+1,exteriorFarfield=far+1,
            sourceVertices=V,sourceFaces=F+1,normalHeightsByColumnM=H,originalNormalHeightsM=h,optionsJSON=opts)
        savemat(str(folder/'prepared_envelope.mat'),data,do_compression=True)
        geom=audit(data);wj=float(min(wedge_jacobian_min(Q[W])));assert wj>0
        write(folder/'shell_geometry_audit.json',dict(geom,minimumPrismJacobian=wj))
        r.update(prisms=len(W),unchangedNearWallPrisms=19*len(F),unchangedTailFanPrisms=len(tail),
            originalWallAndFarfieldExactlyPreserved=True,allPrismConnectivityExactlyPreserved=True,
            preparationSeconds=time.monotonic()-start,minimumPrismJacobian=wj)
        write(folder/'preparation_report.json',r);assert r['preparationSeconds']<=60
        intersect=audit_surface(Q,top,folder/'boundary_audit',dict(source='prepared_envelope.mat',sha256=sha(folder/'prepared_envelope.mat')))
        assert intersect['selfIntersectionFree'],'new envelope self intersects'
        gmshStarted=True;inputs=dict(locked);inputs[str(folder/'prepared_envelope.mat')]=sha(folder/'prepared_envelope.mat')
        native=run([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--exterior-worker',str(folder)],folder/'exterior_process',180,inputs)
        assert native['exitCode']==0 and all(native['sourceStableAfter'].values()),'exterior worker failed; retain original report'
        code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        stable={p:sha(Path(p))==h for p,h in locked.items()};receipt=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,
            error=error,sourceSHA256=locked,sourceStableAfter=stable,sourceStable=all(stable.values()),inputStable=all(stable.values()),
            gmshStarted=gmshStarted,outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},physicsQualified=False,newQualifiedArea=0)
        write(folder/'mesh_process.json',receipt);print(json.dumps(receipt,indent=2));assert all(stable.values())
    raise SystemExit(code)

if __name__=='__main__':main()
