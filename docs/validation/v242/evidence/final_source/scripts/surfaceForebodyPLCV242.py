"""模块职责：成熟Gmsh只重建人工外表面；原4434壁面和46点环严格保持。"""
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from reviewSU2ExplicitQuality import bound

BASE=ROOT/'outputs/v242_validation';RUNTIME=ROOT/'tools/runtime/gmsh/4.15.2'


def shell_audit(p,f):
    """封闭有向单壳与散度体积；不是物理守恒或网格无交叠的充分证明。"""
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);keys,inv,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(count==2),'PLC_EDGE_PAIRING'
    assert np.all(np.bincount(inv,weights=np.sign(edges[:,1]-edges[:,0]))==0),'PLC_ORIENTATION'
    assert np.array_equal(np.unique(f),np.arange(len(p)))
    graph=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(p),len(p))).tocsr()
    assert connected_components(graph,directed=False,return_labels=False)==1,'PLC_ONE_SHELL'
    assert len(np.unique(np.sort(f,axis=1),axis=0))==len(f),'PLC_DUPLICATE_FACE'
    volume=float(np.sum(np.einsum('ij,ij->i',p[f[:,0]],np.cross(p[f[:,1]],p[f[:,2]])))/6)
    assert volume>0
    return dict(closedOrientedSingleShell=True,edges=len(keys),signedVolumeM3=volume)


def surface():
    prep=read(BASE/'prepare01/report.json');assert prep['adoptSingleShell']
    rt=read(RUNTIME/'runtime.json')
    for row in rt['files']:assert sha(RUNTIME/row['path'])==row['sha256']
    z=np.load(str(BASE/'prepare01/wall.npz'));p=z['vertices'];wall=z['wall'];ring=z['ring'];lo,hi=z['bounds'];cut=float(lo[0])
    h0=float(np.median(np.linalg.norm(p[ring]-np.roll(p[ring],-1,axis=0),axis=1)));hmax=float(np.max(hi-lo)/16)
    radius=float(np.mean(np.linalg.norm(p[ring,1:],axis=1)))
    out=BASE/'surface01';out.mkdir(exist_ok=False);sys.path.insert(0,str(RUNTIME));import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.NumThreads',1);gmsh.option.setNumber('Mesh.MaxNumThreads2D',1)
        gmsh.option.setNumber('Mesh.Algorithm',6);gmsh.option.setNumber('Mesh.ElementOrder',1)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0);gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0);gmsh.model.add('forebody-original-wall-open-shell')
        for node in ring:gmsh.model.geo.addPoint(*p[node],meshSize=h0,tag=int(node+1))
        corners=np.array([[a,b,c] for a in (lo[0],hi[0]) for b in (lo[1],hi[1]) for c in (lo[2],hi[2])])
        for i,point in enumerate(corners):gmsh.model.geo.addPoint(*point,meshSize=hmax,tag=len(p)+i+1)
        edge_map={};ringlines=[]
        def line(a,b):
            pair=tuple(sorted((a,b)))
            if pair not in edge_map:edge_map[pair]=gmsh.model.geo.addLine(*pair)
            return edge_map[pair]*(1 if a<b else -1)
        for a,b in zip(ring,np.roll(ring,-1)):ringlines.append(line(int(a+1),int(b+1)))
        inner=gmsh.model.geo.addCurveLoop(ringlines)
        definitions=[([0,1,3,2],0,-1),([4,6,7,5],0,1),([0,4,5,1],1,-1),([2,3,7,6],1,1),([0,2,6,4],2,-1),([1,5,7,3],2,1)]
        surfaces=[]
        for idx,(nodes,axis,sign) in enumerate(definitions):
            tags=[len(p)+n+1 for n in nodes];loop=gmsh.model.geo.addCurveLoop([line(a,b) for a,b in zip(tags,tags[1:]+tags[:1])])
            surf=gmsh.model.geo.addPlaneSurface([loop,inner] if idx==0 else [loop]);surfaces.append(surf)
        gmsh.model.geo.synchronize()
        for tag in ringlines:gmsh.model.mesh.setTransfiniteCurve(abs(tag),2)
        gmsh.model.mesh.setSizeCallback(lambda dim,tag,x,y,z,lc:min(hmax,h0+.25*np.hypot(x-cut,np.hypot(y,z)-radius)))
        gmsh.model.mesh.generate(2)
        tags,coords,_=gmsh.model.mesh.getNodes();coords=coords.reshape(-1,3);faces=[];surface_ids=[]
        for i,surf in enumerate(surfaces):
            kinds,_,conn=gmsh.model.mesh.getElements(2,surf);assert list(kinds)==[2]
            f=conn[0].reshape(-1,3);faces.append(f);surface_ids.extend([i]*len(f))
        faces=np.concatenate(faces);surface_ids=np.asarray(surface_ids)
        np.savez_compressed(str(out/'native_surface.npz'),nodeTags=tags,vertices=coords,faceNodeTags=faces,surfaceID=surface_ids)
        lookup={int(n):i for i,n in enumerate(tags)};mapping={};ring_error=0.
        for node in ring:
            ids,xyz,_=gmsh.model.mesh.getNodes(0,int(node+1));assert len(ids)==1
            error=float(np.max(abs(xyz-p[node])));ring_error=max(ring_error,error);assert error==0
            mapping[int(ids[0])]=int(node)
        for tag in ringlines:
            _,_,conn=gmsh.model.mesh.getElements(1,abs(tag));assert len(conn)==1 and len(conn[0])==2
        free=[int(n) for n in tags if int(n) not in mapping]
        q=np.r_[p,np.array([coords[lookup[n]] for n in free])]
        mapping.update({n:len(p)+i for i,n in enumerate(free)})
        outer=np.array([[mapping[int(n)] for n in f] for f in faces],dtype=np.int64)
        for i,(_,axis,sign) in enumerate(definitions):
            mask=np.flatnonzero(surface_ids==i);f=outer[mask];normal=np.cross(q[f[:,1]]-q[f[:,0]],q[f[:,2]]-q[f[:,0]])
            dot=normal[:,axis]*sign;assert np.all(dot!=0);flip=mask[dot<0]
            outer[flip,1],outer[flip,2]=outer[flip,2].copy(),outer[flip,1].copy()
            expected=(lo if sign<0 else hi)[axis];assert np.max(abs(q[f,axis]-expected))<5e-14
        allfaces=np.r_[wall,outer];markers=np.r_[np.ones(len(wall),dtype=int),np.where(surface_ids==0,2,3)]
        shell=shell_audit(q,allfaces);error=abs(shell['signedVolumeM3']/float(z['expectedVolume'])-1);assert error<2e-10
        ub=bound(q[allfaces])[0];minimum=float(ub.min());assert minimum>=.1,('FIXED_FACE_BOUND',minimum)
        assert np.array_equal(q[wall],p[wall]) and len(wall)==4434
        np.savez_compressed(str(out/'plc.npz'),vertices=q,faces=allfaces,markers=markers,wall=wall,
          originalFaceIDs=z['originalFaceIDs'],originalArea=z['originalArea'],expectedVolume=z['expectedVolume'])
        with (out/'plc.dat').open('w',encoding='ascii',newline='\n') as f:
            f.write('%d %d 0\n'%(len(q),len(allfaces)));np.savetxt(f,np.c_[q,np.arange(1,len(q)+1)],fmt=['%.17g']*3+['%d'])
            np.savetxt(f,np.c_[allfaces,markers],fmt='%d')
        save(out/'report.json',dict(passed=True,vertices=len(q),facets=len(allfaces),originalWallFaces=len(wall),
          outletFaces=int(sum(markers==2)),farfieldFaces=int(sum(markers==3)),shell=shell,volumeRelativeDifference=error,
          originalWallExactlyPreserved=True,originalRingExactlyPreserved=True,maximumRingDifferenceM=ring_error,
          minimumSingleFaceMaxSICN=minimum,h0M=h0,hMaximumM=hmax,growth=.25,algorithm='Gmsh4.15.2 Frontal-Delaunay2D/6',
          gmshCPUThreads=1,physicalQualification=False,sourceSHA256={str(BASE/'prepare01/wall.npz'):sha(BASE/'prepare01/wall.npz')},
          runtimeManifestSHA256=sha(RUNTIME/'runtime.json'),plcSHA256=sha(out/'plc.dat')))
    finally:gmsh.finalize()


if __name__=='__main__':raise SystemExit('v242 fixed task closed; use NEXT_DECISION')
