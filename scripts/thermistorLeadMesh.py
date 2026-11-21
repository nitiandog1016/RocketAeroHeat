"""生成两端带引线的 B 型热敏珠原面和轴向稳态扇区体网，长度在求解器中按横径归一。"""
from collections import Counter
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

from thermistorBead import A, B, D, ROOT, read, sha, write_json

LABELS = ["后引线端面", "后引线", "后接头", "珠面", "前接头", "前引线", "前引线端面"]
SECTORS = 48
WEDGE = 4


SMOOTH_LABELS = ["后极冠", "共同珠面", "前极冠"]


def meridian(smooth=False):
    """椭圆与直线之间构造双切圆弧。各段只保留暴露外壁，切点由两侧共享。"""
    a, b = A / D, B / D
    wire = fillet = 0.027e-3 / D
    length = 0.300e-3 / D

    def normal(t):
        n = np.array([b*np.cos(t), a*np.sin(t)])
        return n / np.linalg.norm(n)

    t = brentq(lambda s: b*np.sin(s) + fillet*normal(s)[1] - wire - fillet, 0, np.pi/2)
    tangent = np.array([a*np.cos(t), b*np.sin(t)])
    center = tangent + fillet*normal(t)
    angle = np.linspace(np.arctan2(-normal(t)[1], -normal(t)[0]), -np.pi/2, 7)
    joint = center + fillet*np.column_stack((np.cos(angle), np.sin(angle)))
    shell = np.column_stack((np.linspace(center[0], center[0]+length, 13), np.full(13, wire)))
    cap = np.column_stack((np.full(5, shell[-1, 0]), np.linspace(wire, 0, 5)))
    polar = np.linspace(np.pi-t, t, 65)
    bead = np.column_stack((a*np.cos(polar), b*np.sin(polar)))
    mirror = np.array([-1, 1])
    parts = [cap[::-1]*mirror, shell[::-1]*mirror, joint[::-1]*mirror, bead, joint, shell, cap]
    for left, right in zip(parts, parts[1:]):
        np.testing.assert_allclose(left[-1], right[0], rtol=0, atol=1e-13)
    wall = np.vstack([part[:-1] for part in parts] + [parts[-1][-1:]])
    region = np.repeat(np.arange(7), [len(part)-1 for part in parts])
    if smooth:
        # 原封保留切点间的每个顶点，包括从接头端取得的最后一排坐标。
        indices = np.flatnonzero(region == 3)
        common = wall[indices[0]:indices[-1]+2]
        count = int(np.ceil(t / (polar[0]-polar[1])))
        angles = np.linspace(0, t, count+1)[:-1]
        end = np.column_stack((a*np.cos(angles), b*np.sin(angles)))
        wall = np.vstack((end*mirror, common, end[::-1]))
        region = np.repeat(np.arange(3), [count, len(indices), count])
    return wall, region


def rotate(xy, rays, closed=False):
    """真实轴线只有一个节点。全周原面首尾共用同一组编号。"""
    axis = xy[:, 1] == 0
    xyz, ids = [], []
    count = 0
    for i, angle in enumerate(rays):
        if closed and i == len(rays)-1:
            ids.append(ids[0])
            continue
        mask = ~axis if i else np.ones(len(xy), dtype=bool)
        row = np.empty(len(xy), dtype=int)
        if i:
            row[axis] = ids[0][axis]
        row[mask] = np.arange(count, count+mask.sum())
        count += mask.sum()
        x, r = xy[mask].T
        xyz.append(np.column_stack((x, r*np.cos(angle), r*np.sin(angle))))
        ids.append(row)
    return np.vstack(xyz), np.array(ids)


def original(smooth=False):
    wall, region = meridian(smooth)
    xyz, rings = rotate(wall, np.linspace(0, 2*np.pi, SECTORS+1), True)
    faces, groups, native = [], [], []
    for k in range(SECTORS):
        folded = lambda i: min(i % (2*WEDGE), 2*WEDGE - i % (2*WEDGE))
        for j, group in enumerate(region):
            corners = [(k,j), (k+1,j), (k+1,j+1), (k,j+1)]
            split = [(0,1,2), (0,2,3)] if folded(k)<folded(k+1) else [(0,1,3), (1,2,3)]
            for tri in split:
                corners_tri = [corners[i] for i in tri]
                f = [rings[i,jj] for i,jj in corners_tri]
                if len(set(f)) < 3:
                    continue
                faces.append(f)
                groups.append(group)
                native.append([(folded(i),jj) for i,jj in corners_tri])
    f = np.asarray(faces)
    tri = xyz[f]*D
    cross = np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area = np.linalg.norm(cross,axis=1)/2
    edges = np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    assert np.all(np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)[1]==2)
    assert len(np.unique(edges,axis=0))==len(edges)
    assert np.all(area>0) and np.all(np.einsum("ij,ij->i",cross,tri.mean(axis=1))>0)
    return dict(V=(xyz*D).tolist(),F=(f+1).tolist(),regionID=(np.array(groups)+1).tolist(),
                regionNames=SMOOTH_LABELS if smooth else LABELS), np.asarray(native), area


def outer_mesh(front, extent):
    """Gmsh 只剖分近壁层以外的子午域，保留层间接口每一条确定的边。"""
    runtime = ROOT/"tools/runtime/gmsh/4.15.2"
    for record in read(ROOT/"third_party/gmsh/source_lock.json")["files"]:
        assert sha(runtime/record["path"])==record["sha256"]
    sys.path.insert(0,str(runtime))
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setNumber("Mesh.Algorithm",6)
        gmsh.model.add("B-thermistor-with-leads")
        geo=gmsh.model.geo
        h=np.linalg.norm(np.diff(front,axis=0),axis=1)
        points=[geo.addPoint(x,r,0,max(.012,min(.08,step))) for (x,r),step in zip(front,np.r_[h,h[-1]])]
        lines=[geo.addLine(a,b) for a,b in zip(points,points[1:])]
        center=geo.addPoint(0,0,0)
        poles=[geo.addPoint(x,r,0,2) for x,r in ((extent,0),(0,extent),(-extent,0))]
        arcs=[geo.addCircleArc(poles[i],center,poles[i+1]) for i in (0,1)]
        loop=lines+[geo.addLine(points[-1],poles[0]),*arcs,geo.addLine(poles[-1],points[0])]
        surface=geo.addPlaneSurface([geo.addCurveLoop(loop)])
        geo.synchronize()
        for line in lines:
            gmsh.model.mesh.setTransfiniteCurve(line,2)
        gmsh.model.mesh.generate(2)
        tags,coords,_=gmsh.model.mesh.getNodes()
        xy=np.array(coords).reshape(-1,3)[:,:2]
        numbering={int(t):i for i,t in enumerate(tags)}
        def elements(dim,tag,kind,n):
            types,_,nodes=gmsh.model.mesh.getElements(dim,tag)
            assert list(types)==[kind]
            return np.array([numbering[int(t)] for t in nodes[0]],int).reshape(-1,n)
        triangles=elements(2,surface,2,3)
        far=np.vstack([elements(1,arc,1,2) for arc in arcs])
        interface=np.array([numbering[int(gmsh.model.mesh.getNodes(0,p)[0][0])] for p in points])
        active=np.unique(triangles)
        return xy,triangles,far,interface,active
    finally:
        gmsh.finalize()


def wall_layers(wall, first):
    d=np.diff(wall,axis=0)
    n=np.column_stack((-d[:,1],d[:,0]))/np.linalg.norm(d,axis=1)[:,None]
    bisector=n[:-1]+n[1:]
    normals=np.vstack([n[0],bisector/np.sum(bisector*n[:-1],axis=1)[:,None],n[-1]])
    levels=np.array([0,.002,.0045,.0075,.011,.015])
    if first==.001:
        levels=np.sort(np.r_[levels,.001,.00325])
    assert first==levels[1]
    layers=np.array([wall+offset*normals for offset in levels])
    layers[:,[0,-1],1]=0
    return layers


def section(first,extent=20,smooth=False,midbody=False,third_grid=False):
    wall,_=meridian(smooth)
    layers=wall_layers(wall,first)
    width=len(wall)
    if smooth:
        from thermistorCapMesh import replace_caps
        wired,_=meridian()
        front=wall_layers(wired,first)[-1]
        background=replace_caps(front,layers[-1],outer_mesh(front,extent))
    else:
        background=outer_mesh(layers[-1],extent)
    if midbody or third_grid:
        from thermistorMidbodyMesh import refine_midbody
        background=refine_midbody(background,wall,third_grid)
    xy,tri,far,interface,active=background
    extra=np.setdiff1d(active,interface)
    translate=np.full(len(xy),-1,int)
    translate[interface]=np.arange(width)+(len(layers)-1)*width
    translate[extra]=np.arange(len(extra))+layers.shape[0]*width
    np.testing.assert_allclose(xy[interface],layers[-1],rtol=0,atol=1e-13)
    points=np.vstack((layers.reshape(-1,2),xy[extra]))
    quads=[]
    for k in range(len(layers)-1):
        for j in range(width-1):
            a=k*width+j
            quads.append([a,a+1,a+width+1,a+width])
    cells=[*quads,*translate[tri].tolist()]
    area_moment=0.
    for c in cells:
        p=points[c];q=np.roll(p,-1,axis=0)
        twice=p[:,0]*q[:,1]-q[:,0]*p[:,1]
        assert twice.sum()>0
        area_moment+=np.sum((p[:,1]+q[:,1])*twice)/6
    assert np.min(points[:,1])>=0 and np.min(translate[far])>=0
    return points,cells,translate[far],area_moment


# 面序取 VTK/SU2 的外向法线，同时用于检查每个内面的相反绕向。
SIDES={4:((0,2,1),(0,1,3),(1,2,3),(2,0,3)),
       5:((0,3,2,1),(0,1,4),(1,2,4),(2,3,4),(3,0,4)),
       6:((0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)),
       8:((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))}
VTK={4:10,5:14,6:13,8:12}


def sweep(c, low, high):
    bottom,top=low[c],high[c]
    on=bottom==top
    if not on.any():
        return [*bottom,*top]
    if len(c)==3:
        if on.sum()==1:
            i=int(np.flatnonzero(on)[0]); a,b=(i+1)%3,(i+2)%3
            return [bottom[a],top[a],top[b],bottom[b],bottom[i]]
        i=int(np.flatnonzero(~on)[0])
        return [*bottom,top[i]]
    a=next(i for i in range(4) if on[i] and on[(i+1)%4])
    b,c,d=[(a+i)%4 for i in (1,2,3)]
    return [bottom[a],bottom[d],top[d],bottom[b],bottom[c],top[c]]


def prepare_mesh(out,first,smooth=False,midbody=False,third_grid=False):
    from collections import defaultdict
    wall,_=meridian(smooth)
    xy,section_cells,far_edges,moment=section(first,smooth=smooth,midbody=midbody,third_grid=third_grid)
    angles=np.arange(WEDGE+1)*2*np.pi/SECTORS
    xyz,ids=rotate(xy,angles)
    cells=[sweep(c,lo,hi) for lo,hi in zip(ids,ids[1:]) for c in section_cells]
    marker={}
    for name,edges in (("wall",np.column_stack((np.arange(len(wall)-1),np.arange(1,len(wall))))),("far",far_edges)):
        for lo,hi in zip(ids,ids[1:]):
            for a,b in edges:
                marker[tuple(sorted(set([lo[a],lo[b],hi[a],hi[b]])))]=name
    for name,row in (("sym0",ids[0]),("sym1",ids[-1])):
        for c in section_cells:
            marker[tuple(sorted(row[c]))]=name
    boundaries=defaultdict(list)
    unmatched={}
    volume=0.;minimum=float("inf")
    # 分类型向量化检查各面相对单元质心的正投影和旋转体积。
    for width,sides in SIDES.items():
        group=np.array([c for c in cells if len(c)==width],int)
        if not len(group):continue
        pts=xyz[group];center=pts.mean(axis=1)
        v=np.zeros(len(group))
        for side in sides:
            faces=group[:,side];p=xyz[faces]
            normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])/2
            if len(side)==4:normal+=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0])/2
            dv=np.einsum("ij,ij->i",normal,p.mean(axis=1)-center)/3
            if np.min(dv)<=0:raise ValueError(f"非正单元锥体积：{width} 节点，面 {side}，最小 {dv.min()}")
            v+=dv
            for face in faces.tolist():
                key=tuple(sorted(face))
                if key in marker:
                    boundaries[marker.pop(key)].append(face)
                elif key in unmatched:
                    old=unmatched.pop(key)
                    index=face.index(old[0])
                    assert all(old[j]==face[(index-j)%len(face)] for j in range(len(face)))
                else:unmatched[key]=face
        volume+=float(v.sum());minimum=min(minimum,float(v.min()))
    assert not unmatched and not marker and np.isclose(volume,moment*WEDGE*np.sin(angles[1]),rtol=2e-12)
    with (out/"mesh.su2").open("x",encoding="ascii") as f:
        f.write(f"NDIME= 3\nNELEM= {len(cells)}\n")
        for cell in cells:f.write(" ".join(map(str,[VTK[len(cell)],*cell]))+"\n")
        f.write(f"NPOIN= {len(xyz)}\n");np.savetxt(f,xyz,fmt="%.17g")
        f.write(f"NMARK= {len(boundaries)}\n")
        for name,faces in sorted(boundaries.items()):
            f.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(faces)}\n")
            for face in faces:f.write(" ".join(map(str,[5 if len(face)==3 else 9,*face]))+"\n")
    body,mapping,areas=original(smooth)
    write_json(out/"body.json",body)
    np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=ids,wallMap=ids[mapping[:,:,0],mapping[:,:,1]],
                        firstMap=ids[mapping[:,:,0],mapping[:,:,1]+len(wall)],area=areas)
    stats=dict(nodes=len(xyz),cells=len(cells),meridianIntervals=len(wall)-1,cellTypes=dict(Counter(map(len,cells))),
               boundaryFaces={k:len(v) for k,v in boundaries.items()},volumeD3=volume,minCellVolumeD3=minimum,
               firstLayerM=first*D,originalAreaM2=float(areas.sum()),originalFaces=len(areas),
               bodySHA256=sha(out/"body.json"),meshSHA256=sha(out/"mesh.su2"),nodeMapSHA256=sha(out/"nodes.npz"))
    write_json(out/"geometry.json",stats)
    return stats
