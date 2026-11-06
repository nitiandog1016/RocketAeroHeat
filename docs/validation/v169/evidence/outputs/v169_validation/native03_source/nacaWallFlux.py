"""模块用途：多等温壁面逐marker复核，保留交界点两侧的热边界贡献。"""
import numpy as np
from analyzeSU2DualVolumes import FACES


def geometry(points,cells,markers):
    """模块一：按中位对偶面构造法向面积，按原棱边选择最大法向余弦邻点。"""
    neighbors=[set() for _ in points]
    for cell in cells:
        for local in FACES[cell[0]]:
            ids=[cell[1+i] for i in local]
            for a,b in zip(ids,ids[1:]+ids[:1]):
                neighbors[a].add(b);neighbors[b].add(a)
    rows={}
    for tag,faces in markers.items():
        if not tag.startswith('wall_'):continue
        normal={}
        for face in faces:
            ids=face[1:];q=points[ids];center=q.mean(axis=0)
            for k,i in enumerate(ids):
                a=q[k];b=(a+q[(k+1)%len(ids)])/2;c=(a+q[(k-1)%len(ids)])/2
                n=-(np.cross(b-a,center-a)+np.cross(center-a,c-a))/2
                normal[i]=normal.get(i,np.zeros(3))+n
        ids=np.array(sorted(normal));normals=np.array([normal[int(i)] for i in ids])
        area=np.linalg.norm(normals,axis=1);unit=normals/area[:,None]
        selected=[]
        for i,n in zip(ids,unit):
            adjacent=np.array(sorted(neighbors[int(i)]))
            displacement=points[adjacent]-points[i]
            selected.append(int(adjacent[np.argmax(displacement@n/np.linalg.norm(displacement,axis=1))]))
        selected=np.array(selected)
        distance=np.abs(np.sum((points[selected]-points[ids])*unit,axis=1))
        assert min(area)>0 and min(distance)>0
        rows[tag]=dict(ids=ids,normal=unit,area=area,normalNeighbor=selected,distance=distance)
    return rows


def evaluate(geometry,temperature,conductivity,wall_values):
    """模块二：与原生弱等温边界一致，正号入壁；不把共享点折成最后一个marker。"""
    assert set(geometry)==set(wall_values)
    result={}
    for tag,g in geometry.items():
        q=conductivity[g['ids']]*(temperature[g['normalNeighbor']]-wall_values[tag])/g['distance']
        result[tag]=dict(qIntoWall=q,dualArea=g['area'],ids=g['ids'],
            powerIntoWall=float(np.dot(q,g['area'])),givenTw=wall_values[tag])
    return result
