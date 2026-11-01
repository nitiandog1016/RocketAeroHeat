"""v140唯一渐变终止候选：原壁和近壁列保持，预检合格才由Gmsh填完整外域。"""
from pathlib import Path
import json,time,traceback,sys
import numpy as np
from scipy.io import loadmat,savemat
from terminatePrismLayers import graph_layers,terminate
from meshSU2HybridVolume import wedge_jacobian_min,RUNTIME
from auditSU2HybridTopology import audit,polyvolume
from reviewSU2ExplicitQuality import quality
from runCGALSurfaceAudit import audit_surface
from prepareCGALRemeshingRuntime import ROOT,sha

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def main():
    start=time.monotonic();folder=ROOT/'outputs/v140_validation/actual';folder.mkdir(parents=True,exist_ok=False)
    source=ROOT/'outputs/v133_validation/fan_workbench';small=ROOT/'docs/validation/v140/fixture/report.json'
    paths=[small,ROOT/'docs/validation/v140/plan.md',Path(__file__).resolve()]+[source/n for n in ['mesh_process.json','mesh_report.json','mesh_input.mat','hybrid_volume.mat']]
    paths += [ROOT/'scripts'/n for n in ['terminatePrismLayers.py','meshSU2HybridVolume.py','auditSU2HybridTopology.py',
        'su2PyramidGeometry.py','reviewSU2ExplicitQuality.py','runCGALSurfaceAudit.py','cgalSurfaceIntersections.py']]
    paths += [RUNTIME/n for n in ['gmsh.py','gmsh-4.15.dll','runtime.json']]
    locked={str(p):sha(p) for p in paths};error=None;code=1;gmshStarted=False
    try:
        assert json.loads(small.read_text())['passed']
        sr=json.loads((source/'mesh_process.json').read_text());assert sr['exitCode']==0 and sr['sourceStable'] and sr['inputStable']
        for n in ['mesh_report.json','hybrid_volume.mat']:assert sha(source/n)==sr['outputSHA256'][n]
        old=loadmat(str(source/'hybrid_volume.mat'));prior=json.loads((source/'mesh_report.json').read_text())
        P=old['vertices'];V=old['sourceVertices'];F=old['sourceFaces'].astype(np.int64)-1;W0=old['wedges'].astype(np.int64)-1
        far=old['farfield'].astype(np.int64)-1;nl=prior['normalLayers'];assert nl==32 and len(F)==19752 and len(V)==9878
        cut=.031232520723574634;K,d=graph_layers(F,len(V),V[:,0]>=cut,19,32);r=terminate(F,W0,K,nl)
        W,Q,T=[r[n] for n in ['wedges','pyramids','tetrahedra']];top=r['outerFaces']
        # 最低19层逐连接/坐标保持；原尾缘扇全部保留，不能改几何取得质量通过。
        assert np.array_equal(W[:19*len(F)],W0[:19*len(F)]) and np.all(K[np.array(prior['rimNodeIDs'])]==32)
        tail=W0[nl*len(F):];assert np.array_equal(W[-len(tail):],tail)
        data=dict(vertices=P,wedges=W+1,pyramids=Q+1,tetrahedra=T+1,farfield=top+1,sourceVertices=V,sourceFaces=F+1)
        save=dict(data,exteriorFarfield=far+1,columnLayers=K,surfaceGraphDistance=d,
            wedgeSourceFaceLayer=r['wedgeSourceFaceLayer']+1,pyramidSourceFaceLayer=r['pyramidSourceFaceLayer']+1,
            tetSourceFaceLayer=r['tetSourceFaceLayer']+1)
        savemat(str(folder/'prepared_transition.mat'),save,do_compression=True)
        geom=audit(data);wj=float(min(wedge_jacobian_min(P[W])));assert wj>0
        write(folder/'shell_geometry_audit.json',dict(geom,minimumPrismJacobian=wj))
        tq=quality(P[T]);bad=np.flatnonzero(tq<.1)
        prep=dict(nodesInOriginalArray=len(P),usedShellNodes=len(np.unique(np.r_[W.ravel(),Q.ravel(),T.ravel()])),
            prisms=len(W),pyramids=len(Q),terminationTets=len(T),columnLayerHistogram={str(k):int(sum(K==k)) for k in np.unique(K)},
            minimumTetSICN=float(min(tq)),tetsBelowSICN01=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=tq[bad].tolist(),
            badTetCentersM=P[T[bad]].mean(axis=1).tolist(),badSourceFaceLayerOneBased=(r['tetSourceFaceLayer'][bad]+1).tolist(),
            noseCutSourceXM=cut,transitionWallXRangeM=[float(min(V[(K>19)&(K<32),0])),float(max(V[(K>19)&(K<32),0]))],
            fixedMinimumLayers=19,unchangedTailFanPrisms=len(tail),originalWallExactlyPreserved=True,allOriginalCoordinatesUnchanged=True,
            protectedNearWallExactlyPreserved=True,maximumNeighborLayerDifference=1,newVertices=0,geometryPassed=True,
            qualityGatePassed=bool(min(tq)>=.1),preparationSeconds=time.monotonic()-start,physicsQualified=False,newQualifiedArea=0)
        write(folder/'preparation_report.json',prep)
        assert prep['preparationSeconds']<=60,'preparation budget exceeded'
        assert prep['qualityGatePassed'],'new termination tetrahedra fail SICN 0.1; do not start exterior meshing'
        # 仅在本轮新壳正性/体积/质量均通过后启动一次完整自交检查与一次Delaunay。
        intersect=audit_surface(P,top,folder/'boundary_audit',dict(source='prepared_transition.mat',sha256=sha(folder/'prepared_transition.mat')))
        assert intersect['selfIntersectionFree']
        opts=json.loads(str(loadmat(str(source/'mesh_input.mat'))['optionsJSON'][0]));ids=np.unique(np.r_[top.ravel(),far.ravel()])
        sys.path.insert(0,str(RUNTIME));import gmsh
        assert gmsh.__version__=='4.15.2';gmshStarted=True;gmsh.initialize();gmsh.option.setNumber('General.NumThreads',2)
        try:
            gmsh.model.add('RocketAeroHeat-graded-termination')
            for tag,faces in [(1,top[:,[0,2,1]]),(2,far)]:
                ns=np.unique(faces);gmsh.model.addDiscreteEntity(2,tag);gmsh.model.mesh.addNodes(2,tag,ns+1,P[ns].ravel())
                gmsh.model.mesh.addElementsByType(tag,2,[],(faces+1).ravel())
            a=gmsh.model.geo.addSurfaceLoop([1]);b=gmsh.model.geo.addSurfaceLoop([2]);vol=gmsh.model.geo.addVolume([b,a]);gmsh.model.geo.synchronize()
            for k,v in [('Mesh.MeshOnlyEmpty',1),('Mesh.Algorithm3D',1),('Mesh.Renumber',0),('Mesh.Optimize',1),
                ('Mesh.MeshSizeMin',opts['outerMinSize']),('Mesh.MeshSizeMax',opts['farSize']),('Mesh.MeshSizeExtendFromBoundary',1)]:gmsh.option.setNumber(k,v)
            gmsh.model.mesh.generate(3);tags,xyz,_=gmsh.model.mesh.getNodes();xyz=xyz.reshape(-1,3);types,_,conn=gmsh.model.mesh.getElements(3,vol)
            savemat(str(folder/'gmsh_raw.mat'),dict(nodeTags=tags,vertices=xyz,elementTypes=types,elements4=conn[0]),do_compression=True)
            assert list(types)==[4]
            position={int(v):i for i,v in enumerate(tags)};assert np.array_equal(xyz[[position[int(i+1)] for i in ids]],P[ids])
            add=~np.isin(tags,ids+1);allP=np.r_[P,xyz[add]];mapping=np.full(int(max(tags))+1,-1,np.int64)
            mapping[ids+1]=ids;mapping[tags[add]]=np.arange(len(P),len(allP));TT=np.r_[T,mapping[conn[0].reshape(-1,4)]]
        finally:gmsh.finalize()
        used=np.unique(np.r_[W.ravel(),Q.ravel(),TT.ravel(),far.ravel()]);remap=np.full(len(allP),-1,np.int64);remap[used]=np.arange(len(used))
        data.update(vertices=allP[used],wedges=remap[W]+1,pyramids=remap[Q]+1,tetrahedra=remap[TT]+1,farfield=remap[far]+1)
        savemat(str(folder/'hybrid_volume.mat'),data,do_compression=True);full=audit(data);write(folder/'geometry_audit.json',full)
        q=quality(data['vertices'][data['tetrahedra']-1]);report=dict(prep,nodes=len(used),tetrahedra=len(TT),totalElements=len(W)+len(Q)+len(TT),
            minimumTetSICN=float(min(q)),tetsBelowSICN01=int(sum(q<.1)),qualityGatePassed=bool(min(q)>=.1))
        write(folder/'candidate_report.json',report)
        assert len(used)<=400000 and report['totalElements']<=900000 and report['qualityGatePassed']
        code=0
    except Exception as e:
        error=repr(e);(folder/'failure.log').write_text(traceback.format_exc());print(error,flush=True)
    finally:
        stable={p:sha(Path(p))==h for p,h in locked.items()}
        receipt=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,error=error,sourceSHA256=locked,
            sourceStableAfter=stable,sourceStable=all(stable.values()),inputStable=all(stable.values()),gmshStarted=gmshStarted,
            outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()},physicsQualified=False,newQualifiedArea=0)
        write(folder/'mesh_process.json',receipt);print(json.dumps(receipt,indent=2));assert all(stable.values())
    raise SystemExit(code)

if __name__=='__main__':main()
