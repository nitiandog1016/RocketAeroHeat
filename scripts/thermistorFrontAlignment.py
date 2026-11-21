"""在带线 B 型热敏珠第三档体网的前中段，保持连接关系并沿壁切向重排节点。"""
import argparse
import os
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[variable]="1"

import numpy as np
from scipy.optimize import lsq_linear
from scipy.sparse import coo_matrix, eye, vstack

from thermistorBead import ROOT, read, sha, write_json
from thermistorLeadMesh import meridian, SECTORS
from thermistorOuterDomain import load_mesh, write_mesh

START, STOP = 54, 70
LOW, HIGH = .015, .100
KEYS = {"mesh.su2":"meshSHA256","nodes.npz":"nodeMapSHA256","body.json":"bodySHA256",
        "case.cfg":"configSHA256","case.json":"caseSHA256","initial_flow.csv":"initialSHA256"}


def wall_frame(points,wall):
    """求距实际分段母线最近的位置、无量纲壁距及局部单位切向。"""
    delta=np.diff(wall,axis=0)
    length=np.linalg.norm(delta,axis=1)
    tangent=delta/length[:,None]
    dist2=np.full(len(points),np.inf)
    station=np.zeros(len(points));segment=np.zeros(len(points),int)
    for j,side in enumerate(tangent):
        offset=points-wall[j]
        along=np.clip(offset@side,0,length[j])
        diff=offset-along[:,None]*side
        sq=np.einsum("ij,ij->i",diff,diff)
        hit=sq<dist2
        dist2[hit]=sq[hit];segment[hit]=j
        station[hit]=j+along[hit]/length[j]
    return station,np.sqrt(dist2),tangent[segment]


def meridian_faces(marks,ids,count):
    index=np.full(count,-1,int)
    index[ids[0]]=np.arange(ids.shape[1])
    tri=index[np.array([f for f in marks["sym0"] if len(f)==3])]
    quad=index[np.array([f for f in marks["sym0"] if len(f)==4])]
    assert tri.min()>=0 and quad.min()>=0
    return tri,quad


def doubled_areas(xy,tri):
    verts=xy[tri]
    return np.cross(verts[:,1]-verts[:,0],verts[:,2]-verts[:,0])


def relocate(xy,tri,wall):
    """最小化目标单元中一条近法向边的切向分量，限制位移并保护全部有向面积。"""
    stations,dist,tangent=wall_frame(xy,wall)
    taper=np.clip(np.minimum(stations-START,STOP-stations),0,1)
    taper*=np.clip(np.minimum(dist-LOW,HIGH-dist)/.008,0,1)
    taper[:8*len(wall)]=0
    free=np.flatnonzero(taper>1e-7)
    edge_pairs=tri[:,[[0,1],[1,2],[2,0]]]
    edge_len=np.linalg.norm(xy[edge_pairs[:,:,1]]-xy[edge_pairs[:,:,0]],axis=2)
    shortest=np.full(len(xy),np.inf)
    for end in (0,1):np.minimum.at(shortest,edge_pairs[:,:,end].ravel(),edge_len.ravel())
    scale=np.zeros(len(xy));scale[free]=.30*shortest[free]*taper[free]
    assert len(free)>0 and np.isfinite(scale).all()
    station_c,dist_c,tan_c=wall_frame(xy[tri].mean(axis=1),wall)
    keep=(station_c>=START)&(station_c<STOP)&(dist_c>=LOW)&(dist_c<=HIGH)
    pairs=edge_pairs[keep];unit=(xy[pairs[:,:,1]]-xy[pairs[:,:,0]])/edge_len[keep,:,None]
    best=np.argmin(abs(np.einsum("ted,td->te",unit,tan_c[keep])),axis=1)
    pairs=pairs[np.arange(len(best)),best]
    lengths=edge_len[keep][np.arange(len(best)),best];tan=tan_c[keep]
    target=-np.einsum("ij,ij->i",xy[pairs[:,1]]-xy[pairs[:,0]],tan)/lengths
    # 未移动节点的系数为零。轻度正则项令无约束方向保持原来的位置。
    size=len(pairs);columns=np.full(len(xy),-1,int);columns[free]=np.arange(len(free))
    rows=[];cols=[];coeff=[]
    for k,sign in ((0,-1),(1,1)):
        node=pairs[:,k];mask=columns[node]>=0
        rows.extend(np.flatnonzero(mask));cols.extend(columns[node[mask]])
        a=sign*scale[node]*np.einsum("ij,ij->i",tangent[node],tan)/lengths
        coeff.extend(a[mask])
    matrix=coo_matrix((coeff,(rows,cols)),shape=(size,len(free))).tocsr()
    objective=vstack((matrix,.025*eye(len(free))),format="csr")
    fit=lsq_linear(objective,np.r_[target,np.zeros(len(free))],bounds=(-1.,1.),tol=1e-8,lsmr_tol=1e-10,max_iter=150)
    if not fit.success:raise ValueError(f"几何最小二乘未结束：{fit.message}")
    move=np.zeros_like(xy);move[free]=(scale[free]*fit.x)[:,None]*tangent[free]
    old=doubled_areas(xy,tri);step=1.
    while np.min(doubled_areas(xy+step*move,tri)/old)<.5:
        step*=.5
    return xy+step*move,dict(freeMeridianNodes=len(free),targetTriangles=size,areaProtectionStep=step,
                             fitStatus=int(fit.status),fitIterations=fit.nit,optimality=float(fit.optimality))


def prepare(base,out):
    key=read(ROOT/"cases/thermistor_b/leadsThird3.json")
    assert all(sha(base/k)==key[v] for k,v in KEYS.items())
    spec=read(base/"case.json")
    assert spec["caseKind"]=="thermistor-B-axial-two-leads" and spec["outerRadiusOverD"]==20
    cells,xyz,marks=load_mesh(base/"mesh.su2")
    with np.load(base/"nodes.npz") as data:
        np.testing.assert_array_equal(data["xyz"],xyz)
        ids=data["ids"];tri,_=meridian_faces(marks,ids,len(xyz))
        xy=xyz[ids[0],:2];new_xy,stats=relocate(xy,tri,meridian()[0])
        changed=np.any(new_xy!=xy,axis=1)
        updated=xyz.copy()
        for k,ring in enumerate(ids):
            angle=2*np.pi*k/SECTORS
            x,r=new_xy[changed].T
            updated[ring[changed]]=np.column_stack((x,r*np.cos(angle),r*np.sin(angle)))
        out.mkdir(parents=True,exist_ok=False)
        write_mesh(out/"mesh.su2",cells,updated,marks)
        np.savez_compressed(out/"nodes.npz",xyz=updated,ids=ids,area=data["area"],wallMap=data["wallMap"],firstMap=data["firstMap"])
    stats.update(baseMeshSHA256=key["meshSHA256"],meshSHA256=sha(out/"mesh.su2"),nodeMapSHA256=sha(out/"nodes.npz"),
                 generatorSHA256=sha(Path(__file__).resolve()),movedMeridianNodes=int(changed.sum()))
    write_json(out/"placement.json",stats)
    return stats


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("base",type=Path);cli.add_argument("output",type=Path)
    a=cli.parse_args();base=a.base.resolve();out=a.output.resolve()
    for folder in (base,out):folder.relative_to(ROOT/"outputs")
    print(prepare(base,out))
