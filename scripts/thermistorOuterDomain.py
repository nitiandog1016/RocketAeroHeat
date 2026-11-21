"""从已生成的 20D 带线热敏珠体网外接 30D 远场，保留原生单元及节点编号。"""
import argparse
from collections import Counter
import os
from pathlib import Path
import uuid

import numpy as np

from runThermistorLeadCase import configuration
from thermistorLeadMesh import SIDES,VTK
from thermistorBead import ENGINE,ENGINE_SHA,ROOT,read,sha,write_json


def load_mesh(path):
    """读取本项目写出的 SU2 网格，保留各面和各体单元的有序连接。"""
    with path.open(encoding="ascii") as f:
        def count(key):
            name,value=next(f).split("=");assert name.strip()==key
            return int(value)
        assert count("NDIME")==3
        cells=[]
        for _ in range(count("NELEM")):
            kind,*cell=map(int,next(f).split());assert VTK[len(cell)]==kind
            cells.append(tuple(cell))
        xyz=np.asarray([list(map(float,next(f).split())) for _ in range(count("NPOIN"))])
        boundary={}
        for _ in range(count("NMARK")):
            tag,name=next(f).split("=");name=name.strip()
            assert tag=="MARKER_TAG" and name not in boundary
            faces=[]
            for _ in range(count("MARKER_ELEMS")):
                kind,*ids=map(int,next(f).split());assert kind=={3:5,4:9}[len(ids)]
                faces.append(tuple(ids))
            boundary[name]=faces
        assert not f.read().strip()
    assert set(boundary)=={"far","sym0","sym1","wall"}
    return cells,xyz,boundary


def write_mesh(path,cells,xyz,boundary):
    with path.open("x",encoding="ascii") as f:
        f.write(f"NDIME= 3\nNELEM= {len(cells)}\n")
        for cell in cells:f.write(" ".join(map(str,(VTK[len(cell)],*cell)))+"\n")
        f.write(f"NPOIN= {len(xyz)}\n");np.savetxt(f,xyz,fmt="%.17g")
        f.write(f"NMARK= {len(boundary)}\n")
        for name,faces in sorted(boundary.items()):
            f.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(faces)}\n")
            for face in faces:f.write(" ".join(map(str,({3:5,4:9}[len(face)],*face)))+"\n")


def append_shell(cells,xyz,boundary,ids):
    far=boundary["far"]
    source=np.unique([v for face in far for v in face]);positions={v:i for i,v in enumerate(source)}
    np.testing.assert_allclose(np.linalg.norm(xyz[source],axis=1),20,rtol=0,atol=1e-11)
    columns=np.flatnonzero(np.isin(ids[0],source))
    assert np.array_equal(np.unique(ids[:,columns]),source)
    interface_columns=np.searchsorted(source,ids[:,columns])
    symmetry={name:set(v for face in boundary[name] for v in face) for name in ("sym0","sym1")}
    new_cells=list(cells);new_boundary={k:list(v) for k,v in boundary.items()}
    nodes=[xyz];meridians=[ids];lower=source
    for layer in range(1,6):
        start=len(xyz)+(layer-1)*len(source);upper=np.arange(start,start+len(source))
        nodes.append(xyz[source]*((20+2*layer)/20))
        meridians.append(upper[interface_columns]);outer=[]
        for face in far:
            local=[positions[v] for v in face];lo=tuple(lower[local]);hi=tuple(upper[local])
            new_cells.append(lo+hi);outer.append(hi)
            for i,a in enumerate(face):
                j=(i+1)%len(face);b=face[j]
                for name,symmetric in symmetry.items():
                    if a in symmetric and b in symmetric:
                        new_boundary[name].append((lo[i],lo[j],hi[j],hi[i]))
        lower=upper
    new_boundary["far"]=outer
    return new_cells,np.vstack(nodes),new_boundary,np.column_stack(meridians)


def ordered(face):
    j=face.index(min(face))
    return face[j:]+face[:j]


def audit(base,expanded):
    """从写出的文件重新读取，逐体检查有序节点、内面双侧相反和标记外向。"""
    old,coords,marks=load_mesh(base/"mesh.su2");cells,xyz,boundary=load_mesh(expanded/"mesh.su2")
    assert cells[:len(old)]==old and boundary["wall"]==marks["wall"]
    np.testing.assert_array_equal(coords,xyz[:len(coords)])
    owners={}
    for cell in cells:
        for side in SIDES[len(cell)]:
            face=tuple(cell[i] for i in side)
            owners.setdefault(tuple(sorted(face)),[]).append(ordered(face))
    for faces in boundary.values():
        for face in faces:assert owners.pop(tuple(sorted(face)))==[ordered(face)]
    for face in marks["far"]:assert len(owners[tuple(sorted(face))])==2
    assert all(len(pair)==2 and pair[0]==ordered(pair[1][::-1]) for pair in owners.values())
    volumes=[]
    for cell in cells[len(old):]:
        center=xyz[list(cell)].mean(axis=0);volume=0.
        for side in SIDES[len(cell)]:
            p=xyz[[cell[i] for i in side]]-center
            for i in range(1,len(p)-1):
                dv=p[0]@np.cross(p[i],p[i+1])/6
                assert dv>0
                volume+=dv
        volumes.append(volume)
    envelope=0.
    for face in marks["far"]:
        p=coords[list(face)]
        envelope+=sum(p[0]@np.cross(p[i],p[i+1])/6 for i in range(1,len(p)-1))
    expected=envelope*((30/20)**3-1)
    assert np.isclose(sum(volumes),expected,rtol=1e-12,atol=0)
    np.testing.assert_allclose(np.linalg.norm(xyz[np.unique([v for f in boundary["far"] for v in f])],axis=1),30,rtol=0,atol=2e-11)
    with np.load(base/"nodes.npz") as a,np.load(expanded/"nodes.npz") as b:
        np.testing.assert_array_equal(b["xyz"],xyz)
        np.testing.assert_array_equal(a["ids"],b["ids"][:,:a["ids"].shape[1]])
        for k in ("wallMap","firstMap","area"):np.testing.assert_array_equal(a[k],b[k])
        np.testing.assert_array_equal(np.unique(b["ids"]),np.arange(len(xyz)))
    return dict(unchangedCells=len(old),unchangedNodeIDs=len(coords),shellCells=len(volumes),shellNodes=len(xyz)-len(coords),
                shellVolumeD3=float(sum(volumes)),expectedShellVolumeD3=float(expected),minShellCellVolumeD3=float(min(volumes)),
                formerFarFacesPaired=len(marks["far"]),newFarFaces=len(boundary["far"]),interiorFaces=len(owners),radiiD=[20,22,24,26,28,30])


def prepare(base,out):
    contract=read(ROOT/"cases/thermistor_b/leadsMidbody3.json")
    fixed={"mesh.su2":"meshSHA256","body.json":"bodySHA256","nodes.npz":"nodeMapSHA256","case.json":"caseSHA256",
           "case.cfg":"configSHA256","initial_flow.csv":"initialSHA256"}
    assert all(sha(base/k)==contract[v] for k,v in fixed.items())
    assert sha(ENGINE)==ENGINE_SHA and (base/"case.cfg").read_text("ascii")==configuration()
    geom=read(base/"geometry.json");spec=read(base/"case.json")
    assert spec["outerRadiusOverD"]==20 and geom["cells"]==contract["cells"] and geom["nodes"]==contract["nodes"]
    out.mkdir(parents=True,exist_ok=False)
    old,xyz,boundary=load_mesh(base/"mesh.su2")
    with np.load(base/"nodes.npz") as z:
        np.testing.assert_array_equal(z["xyz"],xyz)
        cells,new_xyz,marks,ids=append_shell(old,xyz,boundary,z["ids"])
        np.savez_compressed(out/"nodes.npz",xyz=new_xyz,ids=ids,wallMap=z["wallMap"],firstMap=z["firstMap"],area=z["area"])
    write_mesh(out/"mesh.su2",cells,new_xyz,marks)
    checked=audit(base,out);write_json(out/"shell.json",checked)
    # 新建全部远区的均匀冷态。原域从原始冷态读入，原终场不作初值。
    initial=np.loadtxt(base/"initial_flow.csv",skiprows=1,delimiter=",")
    extra=np.zeros((len(new_xyz)-len(xyz),10));extra[:,0]=np.arange(len(xyz),len(new_xyz))
    extra[:,1:4]=new_xyz[len(xyz):];extra[:,5]=-1;extra[:,8]=1-298.15
    with (base/"initial_flow.csv").open() as f:header=f.readline().strip()
    np.savetxt(out/"initial_flow.csv",np.vstack((initial,extra)),fmt=["%d"]+["%.17g"]*9,delimiter=",",header=header,comments="")
    cold=np.loadtxt(out/"initial_flow.csv",skiprows=1,delimiter=",")
    np.testing.assert_array_equal(cold[:len(xyz)],initial);np.testing.assert_array_equal(cold[len(xyz):],extra)
    for name in ("body.json","case.cfg"):os.link(base/name,out/name)
    spec["outerRadiusOverD"]=30;write_json(out/"case.json",spec)
    geom.update(nodes=len(new_xyz),cells=len(cells),cellTypes=dict(Counter(map(len,cells))),boundaryFaces={k:len(v) for k,v in marks.items()},
                minCellVolumeD3=min(geom["minCellVolumeD3"],checked["minShellCellVolumeD3"]),volumeD3=geom["volumeD3"]+checked["shellVolumeD3"],
                meshSHA256=sha(out/"mesh.su2"),nodeMapSHA256=sha(out/"nodes.npz"))
    write_json(out/"geometry.json",geom)
    sources=[Path(__file__).resolve(),ROOT/"scripts/thermistorLeadMesh.py",ROOT/"scripts/runThermistorLeadCase.py",ROOT/"docs/validation/v303/PLAN.md",
             ROOT/"cases/thermistor_b/leadsMidbody3.json",*(out/n for n in ("geometry.json","shell.json","body.json","case.json","initial_flow.csv","nodes.npz"))]
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",threads=1,maximumSeconds=720,
               jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in sources]))
    write_json(out/"inputs.json",dict(base=base.relative_to(ROOT).as_posix(),sourceSHA256={k:sha(base/k) for k in fixed},
               inputs={k:sha(out/k) for k in (*fixed,"geometry.json","shell.json","engine_request.json")},
               executableSHA256=ENGINE_SHA,registeredSources={str(p.relative_to(ROOT)):sha(p) for p in sources}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("base",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();base=a.base.resolve();out=a.output.resolve()
    for folder in (base,out):folder.relative_to(ROOT/"outputs")
    prepare(base,out)
