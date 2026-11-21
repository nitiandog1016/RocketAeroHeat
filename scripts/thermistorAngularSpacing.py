"""在带线 B 型热敏珠的 20D 多面体域内将 48 个周向间隔二分，建立原三角面的积分权重。"""
import argparse
from itertools import combinations
from pathlib import Path
import shutil

import numpy as np

from thermistorBead import ROOT, read, sha, write_json
from thermistorOuterDomain import load_mesh, write_mesh
from thermistorReaction import terminal


def split_faces(faces, side, forward, backward):
    """顶点顺序沿用父单元；周向轴上的节点同时属于左右子单元。"""
    result = []
    for face in faces:
        sector_ids = {side[i] for i in face} - {-1}
        assert len(sector_ids) == 2
        left, right = sorted(sector_ids)
        if left == 0 and right == 47:
            left, right = 47, 0
        assert right == (left + 1) % 48
        result.append(tuple(backward[i] if side[i] == right else i for i in face))
        result.append(tuple(forward[i] if side[i] == left else i for i in face))
    return result


def overlap_weights(target, native, area):
    """在原面三条边内裁剪细三角形，对交叠多边形精确积分三个细面线性形函数。"""
    basis = (target[1:] - target[0]).T
    uv = np.linalg.lstsq(basis, (native-target[0]).T, rcond=None)[0].T
    bary = np.column_stack((1-uv.sum(axis=1),uv))
    bary[np.abs(bary)<1e-13] = 0
    poly = list(np.column_stack((bary,np.eye(3))))
    for k in range(3):
        new = []
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if a[k]>=0:
                new.append(a)
            if (a[k]<0<b[k]) or (b[k]<0<a[k]):
                new.append(a+(b-a)*(a[k]/(a[k]-b[k])))
        poly = new
    weights = np.zeros(3)
    if len(poly)>=3:
        p = np.asarray(poly)[:,3:]
        for j in range(1,len(p)-1):
            t = p[[0,j,j+1]]
            fraction = np.linalg.det(t)
            assert fraction >= -1e-13
            if fraction>0:
                weights += area*fraction*t.mean(axis=0)
    return weights


def surface_mapping(xyz, scale, patches, original_wall, fine_wall):
    """同一物理壁片内的三角形相交，避免跨接头、跨周向接缝及重复计面积。"""
    owner = {tuple(sorted(v)):i for i,f in enumerate(patches) for v in combinations(f,3)}
    groups = [[] for _ in patches]
    for i,t in enumerate(original_wall):
        groups[owner[tuple(sorted(t))]].append(i)
    old = xyz[original_wall]*scale
    fine = xyz[fine_wall]*scale
    area = np.linalg.norm(np.cross(fine[:,1]-fine[:,0],fine[:,2]-fine[:,0]),axis=1)/2
    parents, children, weights = [], [], []
    for group in groups:
        assert len(group) in (1,2)
        for parent in group:
            for sibling in group:
                for child in (2*sibling,2*sibling+1):
                    w = overlap_weights(old[parent], fine[child], area[child])
                    if w.sum() > 1e-14*area[child]:
                        parents.append(parent);children.append(child);weights.append(w)
    return dict(area=area, overlapParent=np.array(parents), overlapChild=np.array(children),
                shapeAreaM2=np.array(weights))


def project(node_values, wall, original_area, parents, children, weights):
    """从原生壁节点的面密度或表压，以实际交叠面积返回各原三角面的平均值。"""
    integral = np.einsum("ij,ij->i", node_values[wall[children]], weights)
    return np.bincount(parents,weights=integral,minlength=len(original_area))/original_area


def build_geometry(base, dest):
    """只生成热敏珠体网、原面和节点映射，供工作台在确认热边界前使用。"""
    cells, coords, marks=load_mesh(base/"mesh.su2",("wall","far"))
    with np.load(base/"nodes.npz") as z:
        rays,orig,first,old_area=(z[k] for k in ("ids","wallMap","firstMap","area"))
        np.testing.assert_array_equal(z["xyz"],coords)
    assert len(rays)==49 and np.array_equal(rays[-1],rays[0])
    off = rays[0]!=rays[1]
    lo, hi = rays[:-1,off], rays[1:,off]
    edges=np.column_stack((lo.ravel(),hi.ravel()))
    new_ids=np.arange(len(coords),len(coords)+len(edges)).reshape(lo.shape)
    xyz=np.vstack((coords,coords[edges].mean(axis=1)))
    mid=rays[:-1].copy();mid[:,off]=new_ids
    fine_rays=np.empty((97,rays.shape[1]),dtype=int)
    fine_rays[::2]=rays;fine_rays[1::2]=mid
    side=np.full(len(coords),-1,dtype=int);side[lo]=np.arange(48)[:,None]
    forward=np.arange(len(coords));backward=forward.copy()
    forward[lo]=new_ids;backward[hi]=new_ids
    side,forward,backward=(a.tolist() for a in (side,forward,backward))
    def halves(items):
        return split_faces(items,side,forward,backward)
    dest.mkdir(parents=True,exist_ok=False)
    fine_marks={k:halves(v) for k,v in marks.items()}
    write_mesh(dest/"mesh.su2",halves(cells),xyz,fine_marks)
    fine_wall=np.asarray(halves(orig))
    spec=read(base/"case.json")
    assert spec["outerRadiusOverD"]==20 and spec["angularCells"]==48 and spec["velocityMps"]==[0,-3,0]
    m=surface_mapping(xyz,spec["diameterM"],marks["wall"],orig,fine_wall)
    np.savez_compressed(dest/"nodes.npz",xyz=xyz,ids=fine_rays,wallMap=fine_wall,
                        firstMap=np.asarray(halves(first)),originalWallMap=orig,originalFirstMap=first,
                        originalArea=old_area,angularEdges=edges,**m)
    for name in ("case.cfg","body.json"):shutil.copyfile(base/name,dest/name)
    spec.update(angularCells=96,reportAngularCells=48,angularConstruction="chord-midpoint")
    write_json(dest/"case.json",spec)
    return coords, xyz, edges


def prepare(base, seed, dest):
    key = read(ROOT/"cases/thermistor_b/sidewindTipSpacing.json")
    base_audit = read(base/"mesh_audit.json")
    for name,h in base_audit["files"].items():
        assert sha(base/name)==h
    assert base_audit["files"]["mesh.su2"]==key["meshSHA256"] and base_audit["files"]["case.json"]==key["caseSHA256"]
    assert base_audit["files"]["body.json"]==key["bodySHA256"] and base_audit["commonPositionNodesExact"]
    receipt=terminal(seed);s=read(seed/"surface/summary.json")
    frame=sha(seed/"restart.csv")
    assert frame==receipt["outputSHA256"]["restart.csv"]==s["sourceFrameSHA256"]==s["energy"]["inputFrameSHA256"]
    assert s["meshSHA256"]==sha(seed/"mesh.su2")==key["meshSHA256"] and s["caseSHA256"]==key["caseSHA256"]
    assert sha(seed/"nodes.npz")==base_audit["files"]["nodes.npz"]
    assert max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    assert abs(s["energy"]["relativeDifference"])<.005 and s["relativeChange20"]<.002
    coords,xyz,edges=build_geometry(base,dest)
    state=np.loadtxt(seed/"restart.csv",skiprows=1,delimiter=",")
    np.testing.assert_array_equal(state[:,0],np.arange(len(coords)))
    np.testing.assert_allclose(state[:,1:4],coords,rtol=0,atol=1e-13)
    initial=np.column_stack((np.arange(len(xyz)),xyz,np.vstack((state[:,4:],state[edges,4:].mean(axis=1)))))
    with (base/"initial_flow.csv").open(encoding="ascii") as f:header=f.readline().strip()
    np.savetxt(dest/"initial_flow.csv",initial,fmt=["%d"]+["%.17g"]*9,delimiter=",",comments="",header=header)
    write_json(dest/"inputs.json",dict(base=str(base),seed=str(seed),seedSHA256=frame,
               baseMeshSHA256=key["meshSHA256"],files={k:sha(dest/k) for k in base_audit["files"]}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("base","seed","dest"):p.add_argument(name,type=Path)
    a=p.parse_args();dest=a.dest.resolve();dest.relative_to(ROOT/"outputs")
    prepare(a.base.resolve(),a.seed.resolve(),dest)
