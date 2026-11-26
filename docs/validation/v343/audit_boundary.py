"""按原 STEP 的共享边识别封盖内侧壁面，用同工况参考坐标复核人工封盖。"""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "docs/validation/v342")]
from check_reference import COORDINATES, distances
from prepareRFZGeometry import LOCAL_STEP, digest
import installGmshRuntime

REFERENCE = ROOT / "outputs/rfz_reference/author_data_20261124/files/Aerodynamics/Stage_1/EUCASS_2023"
FIELDS = {"p": "c9d34f49e3344cfeef89710da27cd1bc678b954316b572593670ca01729cfab5",
          "q": "0fc2a2542406426271ca23071efbeded202231b6145c761cbde9151e8d0fbabd",
          "T": "e95f0553012466cad15b86e72b2d315c0b625d0a75499975a8ef95b4483b8da1"}


def cad_pockets(step, rows):
    """取作者末端的九个平面环。切开内圈共享边后，用 CAD 拓扑遍历每个腔体。"""
    installGmshRuntime.main()
    sys.path.insert(0, str(ROOT / "tools/runtime/gmsh" / installGmshRuntime.VERSION))
    import gmsh

    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setString("Geometry.OCCTargetUnit","M")
        gmsh.model.occ.importShapes(str(step.resolve()))
        gmsh.model.occ.synchronize()
        entities = [tag for _, tag in gmsh.model.getEntities(2)]
        assert entities == [int(row["cadFace"]) for row in rows]
        adjacency = {i:set() for i in entities}
        cuts = set()
        rings = {}
        for row in rows:
            face = int(row["cadFace"])
            occ = gmsh.model.occ
            assert abs(occ.getMass(2,face)-float(row["cadAreaM2"])) < 1e-8
            if row["surfaceType"] != "Plane":
                continue
            center = np.array(occ.getCenterOfMass(2,face))
            _, loops = occ.getCurveLoops(face)
            # 两个出口平面已由冻结 CAD 的原点与作者图 3、4 的末端方位识别。
            if min(abs(center[0]),abs(center[0]-.1)) > 1e-6 or len(loops)!=2:
                continue
            inner = min(loops,key=lambda loop:sum(occ.getMass(1,int(e)) for e in loop))
            adjacent = set()
            for edge in inner:
                cuts.add(int(edge))
                adjacent.update(gmsh.model.getAdjacencies(1,int(edge))[0].tolist())
            adjacent.remove(face)
            assert len(adjacent)==1
            rings[face] = (adjacent.pop(),center)
        assert len(rings)==9
        for _,edge in gmsh.model.getEntities(1):
            pair = np.unique(gmsh.model.getAdjacencies(1,edge)[0])
            if edge not in cuts and len(pair)==2:
                a,b=map(int,pair)
                adjacency[a].add(b)
                adjacency[b].add(a)
        pockets={}
        for rim,(seed,center) in rings.items():
            todo=[seed]
            seen=set()
            while todo:
                tag=todo.pop()
                if tag not in seen:
                    seen.add(tag)
                    todo.extend(adjacency[tag]-seen)
            assert not seen.intersection(rings)
            pockets[rim]=(seen,center)
        all_inner=set.union(*(v[0] for v in pockets.values()))
        assert len(all_inner)==sum(len(v[0]) for v in pockets.values())
        return pockets
    finally:
        gmsh.finalize()


def verify(stl, step, reference):
    spec=json.loads((ROOT/"cases/rfz_rec/geometry.json").read_bytes())
    case_file=ROOT/"cases/rfz_rec/deg190_reference.json"
    case=json.loads(case_file.read_bytes())
    assert digest(reference/"readme")==case["readmeSHA256"]
    assert digest(reference/"visualise_data.py")==case["fieldKeySHA256"]
    ledger=ROOT/"cases/rfz_rec/cad_faces.csv"
    assert digest(ledger)==spec["regionsSHA256"] and digest(stl)==spec["stlSHA256"] and digest(step)==spec["stepSHA256"]
    with ledger.open(encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    pockets=cad_pockets(step,rows)
    dtype=[("normal","<f4",3),("corners","<f4",(3,3)),("attr","<u2")]
    tri=np.fromfile(stl,dtype=dtype,offset=84)["corners"].astype(float)
    verts,indices=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True)
    faces=indices.reshape(-1,3)
    region=np.repeat([int(r["cadFace"]) for r in rows],[int(r["lastFace"])-int(r["firstFace"])+1 for r in rows])
    edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    owners=np.tile(region,3)
    edges.sort(axis=1)
    order=np.lexsort((edges[:,1],edges[:,0]))
    edges=edges[order].reshape(-1,2,2)
    assert np.array_equal(edges[:,0],edges[:,1])
    owners=owners[order].reshape(-1,2)
    areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert abs(areas.sum()-spec["stlAreaM2"])<1e-8

    source=reference/"deg190"
    def array(key,expected):
        path=source/(key+"_data.npy")
        assert digest(path)==expected
        a=np.load(path,allow_pickle=False).reshape(-1)
        assert len(a)==72706 and np.isfinite(a).all()
        return a
    points=np.column_stack([array(k,h) for k,h in zip("XYZ",COORDINATES)])
    values={k:array(k,h) for k,h in FIELDS.items()}
    far=distances(points,tri)>.05
    assert int(far.sum())==1407 and points[:,1].min()>-1e-10
    hits=np.zeros(len(points),dtype=int)
    caps=[]
    for rim,(inner,center) in sorted(pockets.items()):
        # 在独立 STL 拓扑上重找腔体与其余面的全部交界；每一条必须接到同一实体环。
        inside=np.isin(owners,list(inner))
        seam=inside[:,0]!=inside[:,1]
        assert np.all(owners[seam][~inside[seam]]==rim)
        pairs=edges[seam,0]
        nodes,degree=np.unique(pairs,return_counts=True)
        assert len(nodes)>3 and (degree==2).all()
        xyz=verts[nodes]
        assert np.ptp(xyz[:,0])<1e-7 and abs(xyz[0,0]-center[0])<1e-6
        hull=ConvexHull(xyz[:,1:])
        assert len(hull.vertices)==len(nodes)
        cycle=nodes[hull.vertices]
        expected=np.sort(np.column_stack((cycle,np.roll(cycle,-1))),axis=1)
        assert set(map(tuple,expected))==set(map(tuple,pairs))
        near=np.flatnonzero(abs(points[:,0]-center[0])<=.00025)
        signed=points[near,1:]@hull.equations[:,:2].T+hull.equations[:,2]
        # 留出 0.25 mm 的边缘带，避免把唇边上的参考点硬分到人工面。
        within=near[signed.max(axis=1)<-.00025]
        hits[within]+=1
        far_ids=within[far[within]]
        caps.append(dict(rimCADFace=rim,innerCADFaces=sorted(inner),polygonAreaM2=float(hull.volume),
                         referenceInteriorPoints=len(within),farPoints=len(far_ids),
                         maxAxialOffsetFarM=float(np.abs(points[far_ids,0]-center[0]).max()) if len(far_ids) else None,
                         minInPlaneMarginFarM=float(-signed[np.isin(near,far_ids)].max()) if len(far_ids) else None))
    assert hits.max()==1
    groups={"cavity":sorted(set.union(*(v[0] for v in pockets.values()))),"rim":sorted(pockets)}
    groups["other"]=sorted(set(region)-set(groups["rim"])-set(groups["cavity"]))
    totals={}
    for name,ids in groups.items():
        selected=np.isin(region,ids)
        totals[name]=dict(cadFaces=len(ids),faces=int(selected.sum()),areaM2=float(areas[selected].sum()))
    assert sum(v["faces"] for v in totals.values())==spec["faces"]
    field_proof=dict(temperatureMaxDeviationFrom300K=float(abs(values["T"]-300).max()),
                     pressureFinite=int(np.isfinite(values["p"]).sum()),heatFluxFinite=int(np.isfinite(values["q"]).sum()),
                     heatFluxSignsPresent=np.unique(np.sign(values["q"])).astype(int).tolist())
    return dict(schema="rfz-rec-reference-boundary-v1",dataset=case["dataset"],caseSHA256=digest(case_file),stepSHA256=spec["stepSHA256"],
                stlSHA256=spec["stlSHA256"],cadTableSHA256=spec["regionsSHA256"],cadFaces=spec["cadFaces"],
                faces=spec["faces"],areaM2=spec["stlAreaM2"],caps=caps,groups=totals,
                referencePoints=len(points),farPoints=int(far.sum()),farOnCaps=int(np.sum(far & (hits==1))),
                unassignedFarPoints=int(np.sum(far & (hits==0))),planeBandM=.00025,rimBandM=.00025,
                fieldSHA256=FIELDS,fields=field_proof,projectComputedFaces=0,jointExperimentalAreaM2=0.)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--step",type=Path,default=LOCAL_STEP)
    parser.add_argument("--stl",type=Path,default=ROOT/"outputs/rfz_rec_geometry/RFZ-REC_m.stl")
    parser.add_argument("--reference",type=Path,default=REFERENCE)
    args=parser.parse_args()
    result=verify(args.stl,args.step,args.reference)
    with args.output.open("x",encoding="utf8",newline="\n") as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print("CAD/STL closure and reference identity checked:",result["farOnCaps"],"/",result["farPoints"])
