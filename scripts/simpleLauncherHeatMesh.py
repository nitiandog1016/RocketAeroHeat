"""风洞显热网格按距壁距离布点；原速度场保持其自身离散。"""
from pathlib import Path
import sys

import matplotlib.tri as mtri
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]


def heat_grid(xy, ids, spacing):
    sys.path.insert(0,str(ROOT/"tools/runtime/gmsh/4.15.2"))
    import gmsh
    boundary=np.r_[ids[0],ids[1:,-1],ids[-1,-2::-1],ids[-2:0:-1,0]]
    nw=ids.shape[1]-1; nl=ids.shape[0]-1
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal",0)
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.model.add("wind_tunnel_temperature")
        geo=gmsh.model.geo
        pts=[geo.addPoint(*xy[i],0) for i in boundary]
        edges=[geo.addLine(int(a),int(b)) for a,b in zip(pts,np.roll(pts,-1))]
        domain=geo.addPlaneSurface([geo.addCurveLoop(edges)])
        geo.synchronize()
        field=gmsh.model.mesh.field
        dist=field.add("Distance");field.setNumbers(dist,"CurvesList",edges[:nw]);field.setNumber(dist,"Sampling",12)
        size=field.add("MathEval");field.setString(size,"F",f"{spacing}+0.085*F{dist}+0.002*F{dist}^2")
        field.setAsBackgroundMesh(size)
        for key,value in {"Mesh.MeshSizeFromPoints":0,"Mesh.MeshSizeFromCurvature":0,
                          "Mesh.MeshSizeExtendFromBoundary":0,"Mesh.Algorithm":5,
                          "Mesh.MeshSizeMax":.65,"Mesh.Smoothing":3}.items():
            gmsh.option.setNumber(key,value)
        gmsh.model.mesh.generate(2)
        tags,xyz,_=gmsh.model.mesh.getNodes()
        points=xyz.reshape(-1,3)[:,:2]
        lookup=np.full(int(max(tags))+1,-1,int);lookup[tags.astype(int)]=np.arange(len(points))
        triangles=lookup[gmsh.model.mesh.getElementsByType(2)[1].reshape(-1,3)]
        lines=[]
        for k,line in enumerate(edges):
            t,_,_=gmsh.model.mesh.getNodes(1,line,True)
            ns=lookup[t.astype(int)];start=xy[boundary[k]];direction=xy[boundary[(k+1)%len(boundary)]]-start
            ns=ns[np.argsort((points[ns]-start)@direction)]
            lines.append(ns)
        def chain(ls):
            assert all(a[-1]==b[0] for a,b in zip(ls[:-1],ls[1:]))
            return np.r_[np.concatenate([a[:-1] for a in ls]),ls[-1][-1]]
        wet=chain(lines[:nw]);far=chain(lines[nw+nl:nw+nl+nw])[::-1]
        axis=np.unique(np.concatenate(lines[nw:nw+nl]+lines[nw+nl+nw:]))
        boundary_edges=np.sort(np.vstack([np.column_stack((a[:-1],a[1:])) for a in lines]),axis=1)
        te=np.sort(np.vstack([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]]),axis=1)
        unique,counts=np.unique(te,axis=0,return_counts=True)
        assert np.all(counts<=2) and set(map(tuple,unique[counts==1]))==set(map(tuple,boundary_edges))
        wall_length=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xy[ids[0]],axis=0),axis=1))]
        wet_length=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(points[wet],axis=0),axis=1))]
        assert abs(wet_length[-1]-wall_length[-1])<1e-8
        station=[]
        for i,ns in enumerate(lines[:nw]):
            a,b=xy[ids[0,i:i+2]];delta=b-a
            t=(points[ns]-a)@delta/np.dot(delta,delta)
            assert abs(t[0])<1e-9 and abs(t[-1]-1)<1e-9 and np.all(np.diff(t)>0)
            t[0]=0.;station.extend(i+t[:-1])
        wet_station=np.r_[station,float(nw)]
        return points,triangles,wet,far,axis,wet_station
    finally:
        gmsh.finalize()


def sampler(points,triangles,targets):
    """在包含查询点的真实流体三角形中线性插值，界面舍入只接受同位置节点。"""
    mesh=mtri.Triangulation(points[:,0],points[:,1],triangles)
    owner=mesh.get_trifinder()(targets[:,0],targets[:,1])
    missing=owner<0
    corners=triangles[np.maximum(owner,0)]
    a,b,c=np.moveaxis(points[corners],1,0)
    e1=b-a;e2=c-a;delta=targets-a;det=np.cross(e1,e2)
    s=np.cross(delta,e2)/det;t=np.cross(e1,delta)/det
    weights=np.column_stack((1-s-t,s,t))
    if missing.any():
        distance,node=cKDTree(points).query(targets[missing]);assert distance.max()<1e-8
        corners[missing]=node[:,None];weights[missing]=[1.,0.,0.]
    assert weights.min()>-1e-7
    def sample(field):
        return np.einsum('ij,ij...->i...',weights,field[corners])
    return sample
