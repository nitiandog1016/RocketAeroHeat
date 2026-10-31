"""独立检查v108原生SU2单元/边界方向、真实三维首层距离和原面面积。"""
from pathlib import Path
from scipy.io import loadmat
import numpy as np
import json
import argparse
import hashlib


def determinant(p):
    return np.dot(np.cross(p[1]-p[0],p[2]-p[0]),p[3]-p[0])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');args=parser.parse_args();folder=Path(args.folder).resolve()
    file=folder/'normal_volume.mat';data=loadmat(str(file),squeeze_me=True,struct_as_record=False)
    v=data['volume'];s=data['surface'];mesh=folder/'native_input/mesh.su2'
    lines=mesh.read_text().splitlines();assert lines[0]=='NDIME= 3'
    n=int(lines[1].split('=')[1]);cells=[np.fromstring(line,sep=' ',dtype=int) for line in lines[2:2+n]]
    row=2+n;nv=int(lines[row].split('=')[1]);row+=1
    points=np.array([[float(x) for x in line.split()[:3]] for line in lines[row:row+nv]]);row+=nv
    assert np.array_equal(points,v.vertices)
    assert lines[row]=='NMARK= 2';row+=1;markers={}
    for _ in range(2):
        name=lines[row].split('=')[1].strip();nb=int(lines[row+1].split('=')[1]);row+=2
        markers[name]=[np.fromstring(line,sep=' ',dtype=int)[1:]+1 for line in lines[row:row+nb]];row+=nb
    # SU2 v8.6.0 Check_IntElem_Orientation的参考四面体，独立于项目体积分解。
    minima=[];face_owners={}
    hex_faces=((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    prism_faces=((0,1,2),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5))
    for i,cell in enumerate(cells):
        kind=cell[0];ids=cell[1:]+1;p=points[ids-1]
        tests=((0,2,1,3),(3,4,5,2)) if kind==13 else ((0,1,2,5),(0,2,3,7),(4,6,5,1),(4,7,6,3))
        values=[determinant(p[list(t)]) for t in tests];assert min(values)>0;minima.extend(values)
        for local in prism_faces if kind==13 else hex_faces:
            key=tuple(sorted(ids[list(local)]));face_owners.setdefault(key,[]).append(i)
    boundary={key:owners[0] for key,owners in face_owners.items() if len(owners)==1}
    assert all(len(owners) in (1,2) for owners in face_owners.values())
    declared={tuple(sorted(f)) for faces in markers.values() for f in faces}
    assert declared==set(boundary)
    for faces in markers.values():
        for ids in faces:
            p=points[ids-1];cell=cells[boundary[tuple(sorted(ids))]][1:]
            center=points[cell].mean(axis=0)
            for j in range(1,len(ids)-1):
                assert np.dot(np.cross(p[j]-p[0],p[j+1]-p[0]),center-p[0])>0
    wall_ids=np.unique(np.concatenate(markers['wall']));layer_nodes=len(wall_ids)
    assert np.array_equal(wall_ids,np.arange(1,layer_nodes+1))
    distances=[];alignment=[];wall_area=0
    for ids in markers['wall']:
        p=points[ids-1];step=points[ids-1+layer_nodes]-p
        normal=np.cross(p[1]-p[0],p[2]-p[0]);normal/=np.linalg.norm(normal)
        distances.extend(step @ normal)
        move=step.mean(axis=0);alignment.append(np.degrees(np.arcsin(np.clip(move @ normal/np.linalg.norm(move),-1,1))))
        for j in range(1,len(ids)-1):wall_area+=np.linalg.norm(np.cross(p[j]-p[0],p[j+1]-p[0]))/2
    distances=np.array(distances);assert np.min(distances)>1.98e-6 and np.max(distances)<2.001e-6
    assert abs(wall_area/np.sum(s.area)-1)<1e-12
    assert np.array_equal(s.V,v.sourceVertices) and np.array_equal(s.F,v.sourceFaces)
    result=dict(nativeCells=len(cells),nativeNodes=nv,allSU2OrientationTestsPositive=True,boundaryIntoFluid=True,
                minimumSU2ReferenceDeterminant=float(min(minima)),wallArea=float(wall_area),originalFaces=len(s.F),
                actual3DFirstPlaneDistanceMin=float(min(distances)),actual3DFirstPlaneDistanceMax=float(max(distances)),
                actual3DMinimumWallAlignmentDeg=float(min(alignment)),
                meridionalFirstNormalDistance=2e-6,circumferentialChordEffectRetained=True,
                physicsQualified=False,newQualifiedArea=0,inputSHA256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (file,mesh,Path(__file__).resolve())})
    (folder/'independent_geometry_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='inputSHA256'}))


if __name__=='__main__':main()
