"""计算民用热敏珠第三档固定节点单轮换边的可行收益和可逐项验证的上界。"""
import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import sys
import time

for variable in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[variable]="1"

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorFrontAlignment import KEYS,START,STOP,LOW,HIGH,meridian_faces,wall_frame
from thermistorLeadMesh import meridian
from thermistorOuterDomain import load_mesh
from thermistorBead import read,sha

OUTPUTS=ROOT/"outputs/v307_validation"
DEFAULT_BASE=ROOT/"outputs/v304_validation/wired01"
TICKS=10**12


def save(path,value):
    with path.open("x",encoding="utf8") as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False)


def costs(points,tri,wall):
    """返回三角形面积与最近壁法向边的面积加权夹角，延续 v306 的完整目标分母。"""
    p=points[tri]
    station,depth,t=wall_frame(p.mean(axis=1),wall)
    vec=p[:,[1,2,0]]-p
    projection=np.abs(np.einsum("ted,td->te",vec,t))
    normal=np.column_stack((-t[:,1],t[:,0]))
    angle=np.rad2deg(np.arctan2(projection,abs(np.einsum("ted,td->te",vec,normal)))).min(axis=1)
    area=abs(np.cross(vec[:,0],-vec[:,2]))/2
    return area,area*angle,station,depth


def possible_ray_lock(points,quads,wall):
    """仅剔除确定改变交点的换边。四边形边界已有的交点和临界情况从宽保留。"""
    p=points[quads]
    ends=np.array([[0,1],[2,3],[0,2],[2,1],[1,3],[3,0]])
    segments=p[:,ends]
    keep=np.ones(len(quads),bool)
    for j in range(START,STOP):
        t=wall[j+1]-wall[j];t/=np.linalg.norm(t)
        n=np.array([-t[1],t[0]]);origin=(wall[j+1]+wall[j])/2
        s=(segments-origin)@t;h=(segments-origin)@n
        den=s[:,:,1]-s[:,:,0]
        u=np.divide(-s[:,:,0],den,out=np.zeros_like(den),where=abs(den)>1e-14)
        cut=h[:,:,0]+u*(h[:,:,1]-h[:,:,0])
        hit=(abs(den)>1e-14)&(u>=-1e-10)&(u<=1+1e-10)&(cut>LOW+1e-9)&(cut<HIGH-1e-9)
        # 若射线与边重合或只擦过节点，保守上界容纳此四边形。
        tangent_touch=np.any(abs((p-origin)@t)<1e-10,axis=1)
        for source,other in ((0,1),(1,0)):
            extra=hit[:,source] & ~tangent_touch & (cut[:,source]>LOW+8e-9) & (cut[:,source]<HIGH-8e-9)
            for k in (other,2,3,4,5):
                extra &= ~(hit[:,k] & (abs(cut[:,source]-cut[:,k])<=8e-9))
            keep &= ~extra
    return keep


def enumerate_flips(xy,tri,wall):
    area,cost,sc,dc=costs(xy,tri,wall)
    target=(sc>=START)&(sc<STOP)&(dc>=LOW)&(dc<=HIGH)
    denom=math.fsum(cost[target])
    owners=defaultdict(list)
    for k,(a,b,c) in enumerate(tri):
        for u,v,w in ((a,b,c),(b,c,a),(c,a,b)):
            owners[tuple(sorted((u,v)))].append((k,w))
    pairs=[];quads=[]
    for (u,v),rows in sorted(owners.items()):
        if len(rows)!=2:continue
        (l,c),(r,d)=rows
        if target[l] and target[r] and tuple(sorted((c,d))) not in owners:
            pairs.append((l,r));quads.append((u,v,c,d))
    pairs=np.array(pairs);quads=np.array(quads)
    sv,dv,_=wall_frame(xy,wall)
    # 四顶点位于目标区是完全包含的必要条件，上界可包含边界未严格证明的组合。
    inside=(sv>=START-1e-10)&(sv<=STOP+1e-10)&(dv>=LOW-1e-10)&(dv<=HIGH+1e-10)
    vertices=inside[quads].all(axis=1)
    p=xy[quads];a,b,c,d=p.transpose(1,0,2)
    orient1=np.cross(b-a,c-a)*np.cross(b-a,d-a)
    orient2=np.cross(d-c,a-c)*np.cross(d-c,b-c)
    convex=(orient1<0)&(orient2<0)
    before=cost[pairs].sum(axis=1)
    a1,c1,_,_=costs(xy,quads[:,[2,3,0]],wall)
    a2,c2,_,_=costs(xy,quads[:,[3,2,1]],wall)
    np.testing.assert_allclose((a1+a2)[convex],area[pairs].sum(axis=1)[convex],rtol=1e-11,atol=1e-16)
    fraction=(before-c1-c2)/denom
    locked=possible_ray_lock(xy,quads,wall)
    candidate=vertices&convex&locked
    summary=dict(targetTriangles=int(target.sum()),sampleAreaD2=float(area[target].sum()),baseCost=denom,
                 baseMeanDeg=denom/math.fsum(area[target]),adjacentPairs=len(pairs),
                 fourVerticesInBand=int(vertices.sum()),convexInBand=int((vertices&convex).sum()),
                 rayCompatibleUpperSet=int(candidate.sum()),beneficialWithoutRayLock=int((vertices&convex&(fraction>0)).sum()),
                 beneficialWithRayLock=int((candidate&(fraction>0)).sum()))
    return pairs,quads,fraction,vertices,convex,locked,summary


def certified_cover(pairs,gain):
    """对顶点分配非负势，检查每条边的整数收益上界均被两端覆盖。"""
    # 一个刻度为完整目标区原始代价的 10^-12；向有利于换边方向取整并另留一个刻度。
    ticks=np.maximum(0,np.ceil(gain*TICKS).astype(np.int64)+1)
    active=ticks>0
    vertices,compressed=np.unique(pairs[active],return_inverse=True)
    compressed=compressed.reshape(-1,2);size=len(compressed)
    if size==0:return vertices,np.array([],dtype=np.int64),ticks
    matrix=coo_matrix((np.ones(2*size),(np.repeat(np.arange(size),2),compressed.ravel())),shape=(size,len(vertices))).tocsr()
    factor=float(ticks.max())
    fit=linprog(np.ones(len(vertices)),A_ub=-matrix,b_ub=-ticks[active]/factor,bounds=(0,None),method="highs")
    if not fit.success:raise ValueError("三角形势上界优化失败："+fit.message)
    potential=np.ceil(np.maximum(fit.x,0)*factor).astype(np.int64)
    for (l,r),w in zip(compressed,ticks[active]):
        short=int(w-potential[l]-potential[r])
        if short>0:potential[l]+=short
    assert np.all(potential>=0) and np.all(potential[compressed].sum(axis=1)>=ticks[active])
    return vertices,potential,ticks


def run(base,out):
    start=time.perf_counter()
    contract=read(ROOT/"cases/thermistor_b/leadsThird3.json")
    for name,key in KEYS.items():assert sha(base/name)==contract[key],name
    cells,xyz,marks=load_mesh(base/"mesh.su2")
    with np.load(base/"nodes.npz") as data:
        ids=data["ids"];np.testing.assert_array_equal(data["xyz"],xyz)
    tri,_=meridian_faces(marks,ids,len(xyz));xy=xyz[ids[0],:2];wall=meridian()[0]
    np.testing.assert_array_equal(xy[:len(wall)],wall)
    pairs,quads,gain,inside,convex,locked,summary=enumerate_flips(xy,tri,wall)
    eligible=inside&convex&locked
    v,p,w=certified_cover(pairs[eligible],gain[eligible])
    # 用原三角形编号保存全部势，包括孤立三角形。与候选的配对可独立验证。
    prices=np.zeros(len(tri),np.int64);prices[v]=p
    upper=int(prices.sum())
    possible=np.flatnonzero(eligible & (gain>0))
    used=set();selected=[]
    for k in sorted(possible,key=lambda k:(-gain[k],int(k))):
        l,r=map(int,pairs[k])
        if l not in used and r not in used:
            selected.append(k);used.update((l,r))
    chosen=np.array(selected,dtype=int)
    out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/"certificate.npz",triNodeIDs=ids[0,tri],pairs=pairs,quadNodeIDs=ids[0,quads],
                        gainFraction=gain,inside=inside,convex=convex,rayLockUpperSet=locked,
                        prices=prices,eligibleUpperTicks=w,greedyPairs=chosen)
    summary.update(vertices=len(xyz),volumeCells=len(cells),ticksPerUnit=TICKS,coverTicks=upper,upperFraction=upper/TICKS,
                   upperMeanImprovementDeg=summary["baseMeanDeg"]*upper/TICKS,
                   positiveGainSum=float(np.maximum(gain[eligible],0).sum()),
                   greedyPairs=len(chosen),greedyOptimisticFraction=float(gain[chosen].sum()),meetsUpperGate=upper>=TICKS//5,
                   baseMeshSHA256=contract["meshSHA256"],baseMapSHA256=contract["nodeMapSHA256"],
                   certificateSHA256=sha(out/"certificate.npz"),planSHA256=sha(HERE/"PLAN.md"),scriptSHA256=sha(Path(__file__).resolve()),
                   scanSeconds=time.perf_counter()-start)
    save(out/"bound.json",summary)
    print(json.dumps(summary,ensure_ascii=True,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path);p.add_argument("--base",type=Path,default=DEFAULT_BASE)
    a=p.parse_args();out=a.output.resolve();base=a.base.resolve();out.relative_to(OUTPUTS);base.relative_to(ROOT/"outputs")
    run(base,out)
