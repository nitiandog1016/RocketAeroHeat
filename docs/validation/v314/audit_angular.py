"""读回 96 档体域及壁面，用交点枚举独立检验原面上的线性积分权重。"""
import argparse
import importlib.util
from itertools import combinations
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha,write_json
from thermistorOuterDomain import load_mesh

spec=importlib.util.spec_from_file_location("tip_topology",HERE.parent/"v313/audit_mesh.py")
tip=importlib.util.module_from_spec(spec);spec.loader.exec_module(tip)


def cross(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def independent_overlap(old,fine,parents,children):
    """在细三角形坐标系内枚举内含顶点和边交点，对有序凸多边形求面积与一阶矩。"""
    p=fine[children];t=old[parents]
    basis=p[:,1:]-p[:,:1]
    gram=basis@basis.transpose(0,2,1)
    xy=np.linalg.solve(gram,basis@(t-p[:,:1]).transpose(0,2,1)).transpose(0,2,1)
    unit=np.broadcast_to(np.array([[0.,0.],[1.,0.],[0.,1.]]),xy.shape)
    edge=np.roll(xy,-1,axis=1)-xy
    native_edge=np.roll(unit,-1,axis=1)-unit
    tol=2e-12
    old_in=(xy.min(axis=2)>=-tol)&(xy.sum(axis=2)<=1+tol)
    side=cross(edge[:,:,None,:],unit[:,None,:,:]-xy[:,:,None,:])
    fine_in=np.all(side>=-tol,axis=1)
    distance=unit[:,None,:,:]-xy[:,:,None,:]
    divisor=cross(edge[:,:,None,:],native_edge[:,None,:,:])
    regular=abs(divisor)>1e-13
    den=np.where(regular,divisor,1)
    along_old=cross(distance,native_edge[:,None,:,:])/den
    along_fine=cross(distance,edge[:,:,None,:])/den
    both=regular&(along_old>=-tol)&(along_old<=1+tol)&(along_fine>=-tol)&(along_fine<=1+tol)
    crossing=xy[:,:,None,:]+along_old[...,None]*edge[:,:,None,:]
    pts=np.concatenate((xy,unit,crossing.reshape(-1,9,2)),axis=1)
    valid=np.column_stack((old_in,fine_in,both.reshape(-1,9)))
    count=valid.sum(axis=1)
    center=np.sum(np.where(valid[...,None],pts,0),axis=1)/np.maximum(count,1)[:,None]
    first=pts[np.arange(len(pts)),np.argmax(valid,axis=1)]
    pts=np.where(valid[...,None],pts,first[:,None,:])
    d=pts-center[:,None,:];angle=np.arctan2(d[...,1],d[...,0])
    pts=np.take_along_axis(pts,np.argsort(angle,axis=1)[:,:,None],axis=1)
    end=np.roll(pts,-1,axis=1)
    frac=cross(pts-center[:,None,:],end-center[:,None,:])
    frac[count<3]=0.
    mass=frac.sum(axis=1)
    moment=np.einsum("nv,nvi->ni",frac,(center[:,None,:]+pts+end)/3)
    return np.column_stack((mass-moment.sum(axis=1),moment))


def audit(folder):
    locks=read(folder/"inputs.json");base=Path(locks["base"]);seed=Path(locks["seed"])
    assert all(sha(folder/k)==v for k,v in locks["files"].items())
    old,origin,marks=load_mesh(base/"mesh.su2",("far","wall"))
    cells,xyz,fine_marks=load_mesh(folder/"mesh.su2",("far","wall"))
    vol0,_=tip.topology(old,origin,marks)
    vol,min_cone=tip.topology(cells,xyz,fine_marks)
    assert len(cells)==2*len(old) and all(len(fine_marks[k])==2*len(v) for k,v in marks.items())
    vol_error=float(np.max(abs(vol.reshape(-1,2)/vol0[:,None]-.5)))
    assert vol_error<2e-12
    with np.load(folder/"nodes.npz") as z,np.load(base/"nodes.npz") as ref:
        edges,rays,wall,near,orig,first,area0,area,pi,ci,weights=(z[k] for k in (
            "angularEdges","ids","wallMap","firstMap","originalWallMap","originalFirstMap",
            "originalArea","area","overlapParent","overlapChild","shapeAreaM2"))
        np.testing.assert_array_equal(z["xyz"],xyz);np.testing.assert_array_equal(rays[::2],ref["ids"])
        for x,k in ((orig,"wallMap"),(first,"firstMap"),(area0,"area")):np.testing.assert_array_equal(x,ref[k])
    np.testing.assert_array_equal(xyz[:len(origin)],origin)
    np.testing.assert_array_equal(xyz[len(origin):],origin[edges].mean(axis=1))
    axis=rays[0]==rays[1]
    assert len(rays)==97 and np.array_equal(rays[0],rays[-1]) and np.all(rays[:,axis]==rays[0,axis])
    assert len(np.unique(rays[:-1,~axis]))==96*(~axis).sum() and len(np.unique(rays))==len(xyz)
    # 物理边界上只加原边中点；子面有向面积向量的和须逐父面复原。
    for tag in marks:
        for i,f in enumerate(marks[tag]):
            u,v=fine_marks[tag][2*i:2*i+2]
            verts=set(u)|set(v)
            assert {k for k in verts if k<len(origin)}==set(f)
            assert all(set(edges[k-len(origin)])<=set(f) for k in verts if k>=len(origin))
            p=origin[list(f)];n=tip.normals(p[None])[0]
            nsub=tip.normals(xyz[np.asarray([u,v])])
            assert np.max(np.abs(nsub-n/2))<1e-12*np.linalg.norm(n)+1e-15
    body=read(folder/"body.json");c=read(folder/"case.json")
    assert c["angularCells"]==96 and c["reportAngularCells"]==48 and c["outerRadiusOverD"]==20
    coarse_tri=xyz[orig]*c["diameterM"];fine_tri=xyz[wall]*c["diameterM"]
    np.testing.assert_allclose(coarse_tri,np.asarray(body["V"])[np.asarray(body["F"])-1],rtol=0,atol=1e-16)
    # 所有报告三角形和细三角形必须恰好铺满其各自的原生壁片。
    owners=[]
    for world,tri_id,tri,areas in ((marks,orig,coarse_tri,area0),(fine_marks,wall,fine_tri,area)):
        faces=world["wall"]
        assign={tuple(sorted(t)):i for i,f in enumerate(faces) for t in combinations(f,3)}
        ids=np.asarray([assign[tuple(sorted(t))] for t in tri_id]);owners.append(ids)
        normal=tip.normals(tri);native=np.asarray([tip.normals(xyz[np.array([f])])[0] for f in faces])*c["diameterM"]**2
        assert np.all(np.einsum("ij,ij->i",normal,native[ids])<0)
        np.testing.assert_allclose(np.linalg.norm(normal,axis=1),areas,rtol=2e-12,atol=0)
        np.testing.assert_allclose(np.bincount(ids,weights=areas),np.linalg.norm(native,axis=1),rtol=2e-12,atol=0)
        assert len(np.unique(np.sort(tri_id,axis=1),axis=0))==len(tri_id)
    groups=[[] for _ in marks["wall"]]
    for i,k in enumerate(owners[0]):groups[k].append(i)
    pairs=np.array([(p,child) for child,k in enumerate(owners[1]) for p in groups[k//2]])
    pa,ch=pairs.T
    fractions=independent_overlap(coarse_tri,fine_tri,pa,ch)
    stored={(int(p),int(j)):w/area[j] for p,j,w in zip(pi,ci,weights)}
    assert len(stored)==len(weights) and weights.min()>=-1e-25
    comparisons=np.array([stored.pop((p,j),[0.,0.,0.]) for p,j in pairs]);assert not stored
    shape_error=float(abs(comparisons-fractions).max());assert shape_error<2e-11
    target_area=np.bincount(pi,weights=weights.sum(axis=1))
    area_error=float(np.max(abs(target_area/area0-1)));assert area_error<2e-12
    for k in range(3):
        np.testing.assert_allclose(np.bincount(ci,weights=weights[:,k]),area/3,rtol=2e-12,atol=1e-25)
    reconstructed=np.zeros((len(orig),3))
    np.add.at(reconstructed,pi,np.einsum("ik,ikj->ij",weights,fine_tri[ci]))
    centroid_error=float(np.linalg.norm(reconstructed/area0[:,None]-coarse_tri.mean(axis=1),axis=1).max())
    assert centroid_error<1e-15
    # 两档实际首层顶点保持相同位置，新增采样顶点沿相应层的周向弦取中点。
    to_first=np.full(len(xyz),-1,dtype=int);to_first[orig]=first
    ids=np.unique(wall);added=ids[ids>=len(origin)]
    edge_index={tuple(v):i+len(origin) for i,v in enumerate(edges)}
    for v in added:
        u,b=to_first[edges[v-len(origin)]]
        to_first[v]=u if u==b else edge_index[(u,b)]
    np.testing.assert_array_equal(to_first[wall],near)
    state=np.loadtxt(folder/"initial_flow.csv",skiprows=1,delimiter=",")
    prev=np.loadtxt(seed/"restart.csv",skiprows=1,delimiter=",")
    np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_array_equal(state[:,1:4],xyz)
    np.testing.assert_array_equal(state[:len(prev),4:],prev[:,4:])
    np.testing.assert_array_equal(state[len(prev):,4:],prev[edges,4:].mean(axis=1))
    assert np.isfinite(state).all() and min(state[:,9])>=-1e-12 and max(state[:,9])<=1+1e-12
    assert np.all(state[ids,9]==1) and np.max(abs(state[ids,5:8]))<1e-10
    assert not (folder/"mesh_audit.json").exists()
    write_json(folder/"mesh_audit.json",dict(nodes=len(xyz),cells=len(cells),volumeD3=vol.sum(),minFaceConeD3=min_cone,
               maxHalfVolumeFractionError=vol_error,oldNodes=len(origin),oldCells=len(old),nativeWall=len(fine_marks["wall"]),
               nativeFar=len(fine_marks["far"]),nativeWallTriangles=len(area),reportFaces=len(orig),wallNodes=len(ids),
               originalAreaM2=sum(area0),overlaps=len(ci),maxNativeShapeIntegralFractionError=shape_error,
               maxOriginalRelativeAreaError=area_error,maxCentroidErrorM=centroid_error,commonPositionsExact=True,
               axisAndSeamChecked=True,baseMeshSHA256=locks["baseMeshSHA256"],seedSHA256=locks["seedSHA256"],files=locks["files"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("folder",type=Path);audit(p.parse_args().folder.resolve())
