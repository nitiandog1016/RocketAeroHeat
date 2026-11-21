"""用精确二进制坐标定向、独立交点公式及整数不等式复核热敏珠换边决定。"""
import argparse
import io
import json
from pathlib import Path
import sys
import zipfile

from flip_bound import DEFAULT_BASE,HERE,OUTPUTS,START,STOP,LOW,HIGH,TICKS,meridian,wall_frame,load_mesh,np,read,sha,save


def exact_orient(xy,quads):
    """输入 double 的分母均为二的幂，转为共同尺度的任意精度整数后计算行列式符号。"""
    vertices=np.unique(quads)
    denom=max(x.as_integer_ratio()[1] for x in xy[vertices].ravel())
    coords={i:tuple(num*(denom//div) for num,div in (x.as_integer_ratio() for x in xy[i])) for i in vertices}
    def turn(i,j,k):
        a,b,c=(coords[v] for v in (i,j,k))
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    valid=[];small=[]
    for a,b,c,d in quads:
        first,second=turn(c,d,a),turn(c,d,b)
        valid.append(turn(a,b,c)*turn(a,b,d)<0 and first*second<0)
        small.append(min(abs(first),abs(second))/denom**2)
    return np.array(valid),np.array(small)


def measure(xy,tri,wall):
    p=xy[tri];center=p.mean(axis=1);_,_,t=wall_frame(center,wall)
    vec=p[:,[1,2,0]]-p
    # atan2 单独求三条边与切向线的角，转成相对法向的锐角。
    dots=abs(np.einsum("ijk,ik->ij",vec,t))
    crosses=abs(np.cross(vec,t[:,None]))
    angle=90-np.rad2deg(np.arctan2(crosses,dots)).max(axis=1)
    area=abs(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]))/2
    return angle*area


def crossing(xy,pair,origin,normal):
    start=xy[pair[:,0]];edge=xy[pair[:,1]]-start
    den=np.cross(normal,edge);valid=abs(den)>1e-14
    fraction=np.divide(np.cross(start-origin,normal),den,where=valid,out=np.zeros(len(edge)))
    depth=np.divide(np.cross(start-origin,edge),den,where=valid,out=np.zeros(len(edge)))
    hit=valid & (fraction>=-1e-10) & (fraction<=1+1e-10) & (depth>LOW+1e-9) & (depth<HIGH-1e-9)
    return depth,hit


def verify_local_rays(xy,quads,wall):
    """每个四边形的原、候选交点分别与其四条周界的交点合并。"""
    locked=np.ones(len(quads),bool)
    edges=quads[:,[[0,1],[2,3],[0,2],[2,1],[1,3],[3,0]]].reshape(-1,2)
    for strip in range(START,STOP):
        t=wall[strip+1]-wall[strip];n=np.array([-t[1],t[0]])/np.linalg.norm(t)
        depths,hits=crossing(xy,edges,(wall[strip]+wall[strip+1])/2,n)
        d=depths.reshape(-1,6);h=hits.reshape(-1,6)
        for k in np.flatnonzero(locked & h[:,:2].any(axis=1)):
            def merged(j):
                cols=[j,2,3,4,5];v=np.sort(d[k,cols][h[k,cols]])
                return v if len(v)<2 else v[np.r_[True,np.diff(v)>1e-9]]
            a,b=merged(0),merged(1)
            locked[k]=len(a)==len(b) and np.allclose(a,b,rtol=0,atol=1e-9)
    return locked


def verify(base,bound,candidate,selection,audit,out):
    summary=read(bound/"bound.json");chosen=read(selection);geometry=read(audit/"geometry.json")
    assert sha(base/"mesh.su2")==summary["baseMeshSHA256"] and sha(candidate/"mesh.su2")==geometry["candidateMeshSHA256"]
    assert sha(bound/"certificate.npz")==summary["certificateSHA256"]==chosen["certificateSHA256"]
    assert sha(selection)==geometry["selectionSHA256"] and sha(audit/"rays.json")==geometry["raysSHA256"]
    old,xyz,markers=load_mesh(base/"mesh.su2");new,nxyz,new_markers=load_mesh(candidate/"mesh.su2")
    np.testing.assert_array_equal(nxyz,xyz)
    wall=meridian()[0];xy=xyz[:,:2]
    with zipfile.ZipFile(HERE.parent/"v306/evidence/measurements.zip") as z:
        with np.load(io.BytesIO(z.read("baseline/measurements.npz"))) as b:
            target={tuple(sorted(t)) for t in b["triangleNodeIDs"]}
            base_cost=float(np.dot(b["closestEdgeAngleDeg"],b["triangleAreaD2"]))
        original_faces=z.read("original_front_mid_faces.csv")
    with np.load(bound/"certificate.npz") as f:
        tri,pairs,q,prices,gains=(f[k] for k in ("triNodeIDs","pairs","quadNodeIDs","prices","gainFraction"))
        np.testing.assert_array_equal(tri,np.array([p for p in markers["sym0"] if len(p)==3]))
        # 独立枚举全目标邻接，按原节点编号而非脚本内部压缩编号配对。
        owners={}
        for i,t in enumerate(tri):
            for a,b in zip(t,np.roll(t,-1)):
                owners.setdefault(tuple(sorted((a,b))),[]).append(i)
        subset={i for i,t in enumerate(tri) if tuple(sorted(t)) in target}
        expected=set()
        for (a,b),adj in owners.items():
            if len(adj)!=2 or not set(adj).issubset(subset):continue
            l,r=adj;c,=set(tri[l])-{a,b};d,=set(tri[r])-{a,b}
            if tuple(sorted((c,d))) not in owners:expected.add((l,r,a,b,c,d))
        assert expected==set(map(tuple,np.column_stack((pairs,q)))) and len(expected)==len(pairs)
        convex,area2=exact_orient(xy,q);np.testing.assert_array_equal(convex,f["convex"])
        costs=measure(xy,tri,wall)
        gain=(costs[pairs].sum(axis=1)-measure(xy,q[:,[2,3,0]],wall)-measure(xy,q[:,[3,2,1]],wall))/base_cost
        assert np.max(abs(gain-gains))<1/TICKS
        node_ids,back=np.unique(q,return_inverse=True)
        s,h,_=wall_frame(xy[node_ids],wall)
        inside=((s>=START-1e-10)&(s<=STOP+1e-10)&(h>=LOW-1e-10)&(h<=HIGH+1e-10))[back].reshape(q.shape).all(axis=1)
        np.testing.assert_array_equal(inside,f["inside"])
        lock=verify_local_rays(xy,q,wall);allowed=inside&convex&lock
        upper_set=f["inside"]&f["convex"]&f["rayLockUpperSet"]
        assert not (allowed&~upper_set).any() and np.all(prices>=0)
        ticks=f["eligibleUpperTicks"]
        assert np.all(ticks>=np.maximum(0,np.ceil(gain[upper_set]*TICKS))) and np.all(prices[pairs[upper_set]].sum(axis=1)>=ticks)
        assert sum(map(int,prices))==summary["coverTicks"] and summary["upperFraction"]==summary["coverTicks"]/TICKS
        ids=np.asarray(chosen["pairs"],int)
        assert np.all(allowed[ids]) and len(np.unique(pairs[ids]))==2*len(ids)
        assert abs(gain[ids].sum()-geometry["reduction"])<1e-11
        # 通过棱柱上下投影和三角形原编号，定位全部面锥失败是否来自同一个换边。
        with np.load(base/"nodes.npz") as nodes:
            rings=nodes["ids"]
        fold=np.full(len(xyz),-1,int)
        for ring in rings:fold[ring]=rings[0]
        thin=set()
        for failure in geometry["quality"]["nonpositiveCones"]:
            cell=new[failure["cellZeroBased"]];assert list(cell)==failure["nodeIDs"]
            bottom=set(fold[list(cell)])
            matches=[k for k in ids if bottom.issubset(set(q[k]))]
            assert len(matches)==1;thin.add(int(matches[0]))
        minimum=float(area2[sorted(thin)].min()) if thin else None
    ray_rows=read(audit/"rays.json")["rows"]
    assert len(ray_rows)==48
    for marks,key in ((markers,"baseCutsD"),(new_markers,"newCutsD")):
        edges=np.array(sorted({tuple(sorted((l,r))) for f in marks["sym0"] for l,r in zip(f,f[1:]+f[:1])}))
        for row in ray_rows:
            strip=row["strip"];side=wall[strip+1]-wall[strip]
            d,hit=crossing(xy,edges,wall[strip]+row["fraction"]*side,np.array([-side[1],side[0]])/np.linalg.norm(side))
            v=np.r_[LOW,np.sort(d[hit]),HIGH];v=v[np.r_[True,np.diff(v)>1e-9]]
            np.testing.assert_allclose(v,row[key],rtol=0,atol=1e-11)
    hold=[r for r in ray_rows if r["fraction"]!=.5]
    worst=max(hold,key=lambda r:abs(r["widthRatio"]-1))
    out.mkdir(parents=True,exist_ok=False)
    save(out/"verification.json",dict(upperFraction=summary["upperFraction"],candidateFraction=geometry["reduction"],
         completeAdjacencies=len(expected),exactConvexPairs=int((convex&inside).sum()),strictRaysAndVertices=int(allowed.sum()),
         certifiedTicks=summary["coverTicks"],selectedPairs=len(ids),disjointSourceTriangles=len(np.unique(pairs[ids])),
         thinPairIDs=sorted(thin),smallestDoubledMeridianAreaD2=minimum,failedCells=len(geometry["quality"]["nonpositiveCones"]),
         worstWidth={k:worst[k] for k in ("strip","fraction","baseCount","newCount","baseMaxD","newMaxD","widthRatio")},
         raysChecked=len(ray_rows),sourceArchiveSHA256=sha(HERE.parent/"v306/evidence/measurements.zip"),
         geometrySHA256=sha(audit/"geometry.json"),verifierSHA256=sha(Path(__file__).resolve())))
    with zipfile.ZipFile(out/"measurements.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for p in (bound/"certificate.npz",bound/"bound.json",selection,audit/"geometry.json",audit/"rays.json",
                  audit/"measured/measurements.npz",audit/"measured/metrics.json",HERE/"PLAN.md",HERE/"CANDIDATE.md"):
            z.write(p,p.name)
        z.writestr("original_front_mid_faces.csv",original_faces)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("bound","candidate","selection","audit","output"):p.add_argument(k,type=Path)
    p.add_argument("--base",type=Path,default=DEFAULT_BASE);a=p.parse_args()
    out=a.output.resolve();out.relative_to(OUTPUTS)
    verify(a.base.resolve(),a.bound.resolve(),a.candidate.resolve(),a.selection.resolve(),a.audit.resolve(),out)
