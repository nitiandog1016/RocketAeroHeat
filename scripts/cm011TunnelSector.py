"""生成 CM-011 低速等温试件的固定三维对照，沿用同工况的原面和物性。"""
import argparse
from collections import Counter
import json
from pathlib import Path
import re

import numpy as np

from cm011Tunnel import D, ROOT, meridian, volume, write_json
from readCM011Tunnel import evaluate
from runEmbeddedSU2 import sha

ANGLE = np.pi / 30
RAYS = np.array([0., ANGLE/2, ANGLE])
PARENT = ROOT/"outputs/v271_validation/fine02"
REFERENCE = ROOT/"docs/validation/v272/evidence/fine02.json"
HEX_FACES = ((3,2,1,0),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
PRISM_FACES = ((2,1,0),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5))


def mesh(rays=RAYS, level=2, layers=42, wall=None, extent=12., diameter=D):
    """合并轴线节点，沿母线法向和周向生成体单元，并独立核对面邻接。"""
    if wall is None:
        wall, _ = meridian(level)
    xy, quads, _ = volume(wall, layers, extent, level)
    steps = len(rays)-1
    angle = rays[-1]-rays[0]
    delta = angle/steps
    assert 0 < delta < np.pi and np.allclose(np.diff(rays),delta,atol=1e-14,rtol=0)
    width = len(wall)
    closed_ring = abs(angle-2*np.pi) < 1e-12
    ids = np.empty((len(rays), len(xy)), dtype=int)
    xyz, source_node, angles = [], [], []
    for k, phi in enumerate(rays):
        # 全周计算只保留一组接缝节点，首末两排体单元共享真实内面。
        if closed_ring and k == steps:
            ids[k] = ids[0]
            continue
        for i, (x, r) in enumerate(xy):
            if k and i % width in (0, width-1):
                ids[k,i] = ids[0,i]
            else:
                ids[k,i] = len(xyz)
                xyz.append((x,r*np.cos(phi),r*np.sin(phi)))
                source_node.append(i)
                angles.append(phi)
    xyz = np.asarray(xyz)
    source_node = np.asarray(source_node)
    cells, owners, boundary = [], {}, {}
    total_volume, min_volume = 0., float("inf")
    for k in range(steps):
        for q in quads:
            a,b,c,d = ids[k,q]
            e,f,g,h = ids[k+1,q]
            if a == e:
                nodes = [a,f,b,d,g,c]
            elif b == f:
                nodes = [b,a,e,c,d,h]
            else:
                nodes = [a,b,c,d,e,f,g,h]
            surface = PRISM_FACES if len(nodes)==6 else HEX_FACES
            p = xyz[nodes]
            center = p.mean(axis=0)
            cell_volume = 0.
            for index in surface:
                face = tuple(nodes[j] for j in index)
                pts = xyz[list(face)]
                area_vector = .5*np.cross(pts[1]-pts[0],pts[-1]-pts[0])
                if len(face)==4:
                    area_vector += .5*np.cross(pts[2]-pts[1],pts[3]-pts[1])
                    # 两个三角形取一致对角线。首三角形为 0,1,3。
                projected = np.dot(area_vector, pts.mean(axis=0)-center)
                if projected <= 0:
                    raise ValueError(f"翻转的三维单元 {len(cells)}，原节点 {q.tolist()}")
                cell_volume += projected/3
                key = tuple(sorted(face))
                if key in owners:
                    old = owners[key]
                    if old is None or not np.allclose(old,-area_vector,rtol=1e-10,atol=1e-14):
                        raise ValueError("三维内面邻接或方向错误")
                    owners[key] = None
                    del boundary[key]
                else:
                    owners[key] = area_vector
                    boundary[key] = face
            total_volume += cell_volume
            min_volume = min(min_volume,cell_volume)
            cells.append(nodes)
    # 闭合母线四边形的一阶矩乘以周向弦面因子，独立得到几何体积。
    p = xy[quads]
    end = np.roll(p,-1,axis=1)
    moment = np.sum((p[:,:,1]+end[:,:,1])*(p[:,:,0]*end[:,:,1]-end[:,:,0]*p[:,:,1]))/6
    reference_volume = moment*steps*np.sin(delta)
    assert np.isclose(total_volume,reference_volume,rtol=1e-12,atol=0)
    markers = {key:[] for key in (("wall","far") if closed_ring else ("wall","far","sym0","sym1"))}
    plane0, plane1 = set(ids[0]), set(ids[-1])
    for face in boundary.values():
        layer = source_node[list(face)]//width
        if np.all(layer==0): name="wall"
        elif np.all(layer==layers): name="far"
        elif not closed_ring and plane0.issuperset(face): name="sym0"
        elif not closed_ring and plane1.issuperset(face): name="sym1"
        else: raise ValueError("出现未归属的三维边界面")
        markers[name].append(face)
    if closed_ring:
        assert all(owners[tuple(sorted(ids[0,q]))] is None for q in quads)
    else:
        assert len(markers["sym0"])==len(markers["sym1"])==len(quads)
    assert len(markers["wall"])==len(markers["far"])==steps*(width-1)
    counts = Counter(len(c) for c in cells)
    counts_by_marker = {k:len(v) for k,v in markers.items()}
    wall_area = sum(np.linalg.norm(owners[tuple(sorted(f))]) for f in markers["wall"])
    mid_r = .5*(wall[1:,1]+wall[:-1,1])
    dx,dr = np.diff(wall,axis=0).T
    chord_area = np.sum(2*steps*np.sin(delta/2)*mid_r*np.hypot(dx,dr*np.cos(delta/2)))
    assert np.isclose(wall_area,chord_area,rtol=1e-12,atol=0)
    copies = round(2*np.pi/angle)
    assert abs(angle*copies-2*np.pi)<1e-12
    stats = dict(nodes=len(xyz), hexahedra=counts[8], prisms=counts[6], sharedAxisNodes=2*(layers+1),
                 boundaryFaces=counts_by_marker, minCellVolumeD3=min_volume, volumeD3=total_volume,
                 meridianMomentVolumeD3=reference_volume, chordWallAreaM2=wall_area*diameter**2,
                 fullChordWallAreaM2=copies*wall_area*diameter**2,
                 fullRotatedWallAreaM2=2*np.pi*np.sum(mid_r*np.hypot(dx,dr))*diameter**2)
    return xy,xyz,source_node,np.array(angles),ids,cells,markers,stats


def write_volume(folder,xyz,cells,markers):
    with (folder/"mesh.su2").open("x",encoding="ascii") as f:
        f.write(f"NDIME= 3\nNELEM= {len(cells)}\n")
        for cell in cells:
            f.write(" ".join(map(str,[13 if len(cell)==6 else 12,*cell]))+"\n")
        f.write(f"NPOIN= {len(xyz)}\n")
        np.savetxt(f,xyz,fmt="%.17g")
        f.write(f"NMARK= {len(markers)}\n")
        for marker,faces in markers.items():
            f.write(f"MARKER_TAG= {marker}\nMARKER_ELEMS= {len(faces)}\n")
            for face in faces:
                f.write(" ".join(map(str,[5 if len(face)==3 else 9,*face]))+"\n")


def prepare(folder):
    refs = json.loads(REFERENCE.read_bytes())
    for n,digest in refs["sourceSHA256"].items():
        assert sha(PARENT/n)==digest, n
    parent_stats,_ = evaluate(PARENT)
    assert parent_stats["acceptedNumericalWall"] and parent_stats["originalFaces"]==14140
    folder.mkdir(parents=True,exist_ok=False)
    xy,xyz,source,phi,ids,cells,markers,stats = mesh()
    state = np.loadtxt(PARENT/"restart.csv",delimiter=",",skiprows=1)
    assert np.array_equal(state[:,0],np.arange(len(xy)))
    np.testing.assert_allclose(state[:,1:3],xy,atol=1e-12,rtol=0)
    pressure,u,v,h = state[source,3:].T
    write_volume(folder,xyz,cells,markers)
    fields = np.column_stack((np.arange(len(xyz)),xyz,pressure,u,v*np.cos(phi),v*np.sin(phi),h))
    np.savetxt(folder/"initial_flow.csv",fields,fmt=["%d"]+8*["%.17g"],delimiter=",",comments="",
               header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy"')
    config = (PARENT/"case.cfg").read_text("ascii")
    replacements = {"AXISYMMETRIC":"NO","MARKER_SYM":"(sym0, sym1)","ITER":"16000",
                    "CONV_RESIDUAL_MINVAL":"-8.5","OUTPUT_WRT_FREQ":"(1000, 1000, 16000)"}
    for key,value in replacements.items():
        config,num = re.subn("^"+key+r"= .*?$",key+"= "+value,config,flags=re.M)
        assert num==1
    config = config.replace("RMS_VELOCITY-Y, RMS_ENTHALPY","RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY")
    (folder/"case.cfg").write_text(config,encoding="ascii")
    body = json.loads((PARENT/"body.json").read_bytes())
    orig = np.array(body["V"])/D
    wall,_ = meridian(2)
    np.testing.assert_allclose(np.column_stack((orig[:,0],np.hypot(orig[:,1],orig[:,2]))),
                               wall[body["nodeRow"]],rtol=0,atol=1e-12)
    stats.update(parentBodySHA256=sha(PARENT/"body.json"), parentRestartSHA256=sha(PARENT/"restart.csv"),
                 originalFaces=len(body["F"]), originalVertices=len(orig), sectorDegrees=6, rotationCopies=60)
    write_json(folder/"geometry.json",stats)
    np.savez_compressed(folder/"nodes.npz",xyz=xyz,ids=ids)
    locked = [Path(__file__),REFERENCE,PARENT/"body.json",folder/"geometry.json",folder/"nodes.npz"]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               maximumSeconds=720,threads=1,sourceFiles=[p.resolve().relative_to(ROOT).as_posix() for p in locked]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_output",type=Path)
    path = parser.parse_args().new_output.resolve()
    path.relative_to(ROOT/"outputs")
    prepare(path)
