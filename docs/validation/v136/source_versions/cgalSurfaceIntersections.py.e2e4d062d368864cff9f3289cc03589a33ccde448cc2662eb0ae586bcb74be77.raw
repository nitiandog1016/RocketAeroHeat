"""固定CGAL进程只读检查三角表面，保持输入节点/面身份，不修补或改变坐标。"""
from pathlib import Path
import argparse,json,hashlib,time
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL import CGAL_Polygon_mesh_processing as PMP

def audit(path):
    start=time.monotonic();lines=path.read_text().splitlines();assert lines[0]=='OFF'
    nv,nf,_=map(int,lines[1].split());points=[tuple(map(float,row.split())) for row in lines[2:2+nv]]
    faces=[tuple(map(int,row.split()[1:])) for row in lines[2+nv:]]
    assert len(faces)==nf and all(len(f)==3 for f in faces)
    pointIDs={p:i for i,p in enumerate(points)};faceIDs={tuple(sorted(f)):i for i,f in enumerate(faces)}
    assert len(pointIDs)==nv and len(faceIDs)==nf
    mesh=Polyhedron_3(str(path));assert mesh.is_valid() and mesh.size_of_vertices()==nv and mesh.size_of_facets()==nf
    for vertex in mesh.vertices():
        p=vertex.point();key=(p.x(),p.y(),p.z());assert key in pointIDs;vertex.set_id(pointIDs[key])
    seen=set()
    for face in mesh.facets():
        half=face.halfedge();ids=[]
        for j in range(3):ids.append(half.vertex().id());half=half.next()
        assert half==face.halfedge();index=faceIDs[tuple(sorted(ids))];original=faces[index]
        assert any(tuple(ids[k:]+ids[:k])==original for k in range(3));assert index not in seen;seen.add(index);face.set_id(index)
    assert len(seen)==nf
    pairs=[];PMP.self_intersections(mesh,pairs)
    bad=sorted({tuple(sorted((p.first.id(),p.second.id()))) for p in pairs})
    return dict(method='CGAL 6.0.1 Polygon_mesh_processing::self_intersections',vertices=nv,faces=nf,
        originalCoordinatesAndOrientedFacesExactlyMatched=True,closed=mesh.is_closed(),pairsZeroBased=bad,pairCount=len(bad),
        selfIntersectionFree=not bad,seconds=time.monotonic()-start,inputSHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
        physicsQualified=False,newQualifiedArea=0)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert not args.output.exists();result=audit(args.input);args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps(result))
if __name__=='__main__':main()
