"""从第三档与候选成品 SU2 逐边量取角度，在十六条固定壁法线上比较穿越间距。"""
import argparse
import json
import os
from pathlib import Path
import sys

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[variable]="1"
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorFrontAlignment import START,STOP,LOW,HIGH,KEYS,meridian_faces,wall_frame
from thermistorLeadMesh import SIDES,meridian
from thermistorOuterDomain import load_mesh
from thermistorBead import read,sha
from thermistorReaction import terminal

OLD=ROOT/"outputs/v304_validation/wired01"
BASE=ROOT/"outputs/v306_validation"


def save(path,value):
    with path.open("x",encoding="utf8") as f:json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False)


def ray_intervals(xy,edges,wall):
    a=xy[edges[:,0]];direction=xy[edges[:,1]]-a
    rays=[]
    for j in range(START,STOP):
        side=wall[j+1]-wall[j]
        origin=(wall[j+1]+wall[j])/2
        normal=np.array([-side[1],side[0]])/np.linalg.norm(side)
        den=np.cross(normal,direction);use=np.abs(den)>1e-14
        away=a[use]-origin
        depth=np.cross(away,direction[use])/den[use]
        fraction=np.cross(away,normal)/den[use]
        hits=np.sort(depth[(fraction>=-1e-10)&(fraction<=1+1e-10)&(depth>LOW+1e-9)&(depth<HIGH-1e-9)])
        cuts=np.r_[LOW,hits,HIGH];cuts=cuts[np.r_[True,np.diff(cuts)>1e-9]]
        intervals=np.diff(cuts)
        assert len(intervals)>0 and abs(intervals.sum()-(HIGH-LOW))<1e-12
        rays.append(dict(strip=j,originD=origin.tolist(),normal=normal.tolist(),count=len(intervals),
                         maximumD=float(intervals.max()),intervalsD=intervals.tolist()))
    return rays


def metrics(folder,out):
    cells,xyz,marks=load_mesh(folder/"mesh.su2")
    with np.load(folder/"nodes.npz") as z:
        np.testing.assert_array_equal(xyz,z["xyz"]);ids=z["ids"]
    xy=xyz[ids[0],:2];wall=meridian()[0]
    np.testing.assert_array_equal(wall,xy[:len(wall)])
    tri,quad=meridian_faces(marks,ids,len(xyz))
    assert len(tri)+len(quad)==len(marks["sym0"])
    sides=tri[:,[[0,1],[1,2],[2,0]]]
    edges=np.unique(np.sort(np.vstack((sides.reshape(-1,2),quad[:,[[0,1],[1,2],[2,3],[3,0]]].reshape(-1,2))),axis=1),axis=0)
    pos,depth,tan=wall_frame(xy[edges].mean(axis=1),wall)
    inside=(pos>=START)&(pos<STOP)&(depth>=LOW)&(depth<=HIGH)
    e=edges[inside];vec=xy[e[:,1]]-xy[e[:,0]]
    sine=abs(np.einsum("ij,ij->i",vec,tan[inside]))/np.linalg.norm(vec,axis=1)
    edge_angle=np.rad2deg(np.arcsin(np.clip(sine,0,1)))
    pos,depth,tan=wall_frame(xy[tri].mean(axis=1),wall)
    samples=(pos>=START)&(pos<STOP)&(depth>=LOW)&(depth<=HIGH)
    t=tri[samples];delta=xy[t[:,[1,2,0]]]-xy[t]
    sine=abs(np.einsum("ted,td->te",delta,tan[samples]))/np.linalg.norm(delta,axis=2)
    best=np.rad2deg(np.arcsin(np.clip(sine.min(axis=1),0,1)))
    area=abs(np.cross(delta[:,0],-delta[:,2]))/2
    rays=ray_intervals(xy,edges,wall)
    widths=np.concatenate([r["intervalsD"] for r in rays])
    out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/"measurements.npz",edgeNodeIDs=ids[0,e],edgeAngleDeg=edge_angle,triangleNodeIDs=ids[0,t],
                        closestEdgeAngleDeg=best,triangleAreaD2=area)
    stats=dict(meshSHA256=sha(folder/"mesh.su2"),mapSHA256=sha(folder/"nodes.npz"),nodes=len(xyz),cells=len(cells),
               targetEdges=len(e),edgeQ10Q50Q90Deg=np.percentile(edge_angle,[10,50,90]).tolist(),
               edgeFractionWithin15Deg=float(np.mean(edge_angle<=15)),targetTriangles=len(t),sampleAreaD2=float(area.sum()),
               closestEdgeMeanDeg=float(np.average(best,weights=area)),closestEdgeMedianDeg=float(np.median(best)),
               intervalMedianD=float(np.median(widths)),intervalP90D=float(np.percentile(widths,90)),rays=rays,
               measurementSHA256=sha(out/"measurements.npz"))
    save(out/"metrics.json",stats)
    return stats


def volumes(xyz,cells):
    total=0.;minimum=float("inf")
    for width,sides in SIDES.items():
        selected=np.array([c for c in cells if len(c)==width])
        if not len(selected):continue
        p=xyz[selected];center=p.mean(axis=1);vol=np.zeros(len(selected))
        for side in sides:
            f=p[:,side]
            n=np.cross(f[:,1]-f[:,0],f[:,2]-f[:,0])/2
            if len(side)==4:n+=np.cross(f[:,2]-f[:,0],f[:,3]-f[:,0])/2
            dv=np.einsum("ij,ij->i",n,f.mean(axis=1)-center)/3
            assert np.all(dv>0),"非正面锥体积"
            vol+=dv
        total+=vol.sum();minimum=min(minimum,vol.min())
    return float(total),float(minimum)


def compare(base,new,reference,out):
    old=read(reference/"metrics.json")
    assert old["meshSHA256"]==sha(base/"mesh.su2") and old["measurementSHA256"]==sha(reference/"measurements.npz")
    curr=metrics(new,out)
    oc,ox,om=load_mesh(base/"mesh.su2");nc,nx,nm=load_mesh(new/"mesh.su2")
    assert oc==nc and om==nm and ox.shape==nx.shape
    with np.load(base/"nodes.npz") as a,np.load(new/"nodes.npz") as b:
        for k in ("ids","wallMap","firstMap","area"):np.testing.assert_array_equal(a[k],b[k])
        rings=a["ids"];wall=meridian()[0];protected=rings[:,:8*len(wall)]
        np.testing.assert_array_equal(ox[protected],nx[protected])
    for label in ("wall","far"):
        ids=np.unique([v for face in om[label] for v in face]);np.testing.assert_array_equal(ox[ids],nx[ids])
    moved=np.any(ox!=nx,axis=1);changed=np.array([moved[list(c)].any() for c in oc])
    vertices=np.unique([v for c,hit in zip(oc,changed) if hit for v in c])
    def footprint(xyz,ids):
        xy=np.column_stack((xyz[ids,0],np.linalg.norm(xyz[ids,1:],axis=1)))
        s,d,_=wall_frame(xy,wall)
        return dict(station=[float(s.min()),float(s.max())],wallDistanceD=[float(d.min()),float(d.max())])
    moving=footprint(ox,np.flatnonzero(moved))
    assert moving["station"][0]>=START and moving["station"][1]<=STOP
    assert moving["wallDistanceD"][0]>=LOW and moving["wallDistanceD"][1]<=HIGH
    # 检查每个周向副本的角度和轴向/径向位移；邻接单元可以跨过目标带边界。
    for k,ring in enumerate(rings):
        hits=moved[ring];np.testing.assert_array_equal(hits,moved[rings[0]])
        p=nx[ring[hits]];q=nx[rings[0,hits]]
        np.testing.assert_allclose(p[:,0],q[:,0],rtol=0,atol=1e-15)
        np.testing.assert_allclose(p[:,1:],q[:,1,None]*[[np.cos(k*np.pi/24),np.sin(k*np.pi/24)]],rtol=0,atol=1e-15)
    volume,minimum=volumes(nx,nc);old_volume,old_min=volumes(ox,oc)
    assert np.isclose(volume,old_volume,rtol=2e-12,atol=0)
    assert np.isclose(old_volume,read(base/"geometry.json")["volumeD3"],rtol=2e-12,atol=0)
    count_ratio=[c["count"]/o["count"] for o,c in zip(old["rays"],curr["rays"])]
    width_ratio=[c["maximumD"]/o["maximumD"] for o,c in zip(old["rays"],curr["rays"])]
    checks=dict(direction=curr["closestEdgeMeanDeg"]<=.8*old["closestEdgeMeanDeg"],
                normalEdgeFraction=curr["edgeFractionWithin15Deg"]>old["edgeFractionWithin15Deg"],
                pooledMedian=.85<=curr["intervalMedianD"]/old["intervalMedianD"]<=1.15,
                pooledP90=.85<=curr["intervalP90D"]/old["intervalP90D"]<=1.15,
                rayCounts=all(.80<=r<=1.20 for r in count_ratio),rayMaxima=all(.75<=r<=1.25 for r in width_ratio))
    report=dict(passed=all(checks.values()),checks=checks,countRatios=count_ratio,maxIntervalRatios=width_ratio,
                movedNodes=int(moved.sum()),movedMeridianNodes=int(moved[rings[0]].sum()),changedCells=int(changed.sum()),
                unchangedCells=int((~changed).sum()),movedNodeRange=moving,changedCellRangeBefore=footprint(ox,vertices),
                changedCellRangeAfter=footprint(nx,vertices),volumeD3=volume,oldVolumeD3=old_volume,minCellVolumeD3=minimum,
                oldMinCellVolumeD3=old_min,baseMesh=old["meshSHA256"],candidateMesh=curr["meshSHA256"],
                baseMetricsSHA256=sha(reference/"metrics.json"),candidateMetricsSHA256=sha(out/"metrics.json"),
                planSHA256=sha(HERE/"PLAN.md"),auditSHA256=sha(Path(__file__).resolve()))
    save(out/"geometry.json",report)
    print("Geometry decision:",checks,"passed",report["passed"])
    print("Angles, normal fractions:",[(m["closestEdgeMeanDeg"],m["edgeFractionWithin15Deg"]) for m in (old,curr)])
    print("Ray count ratios",min(count_ratio),max(count_ratio),"maximum width ratios",min(width_ratio),max(width_ratio))
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("action",choices=("baseline","compare"))
    p.add_argument("output",type=Path);p.add_argument("--candidate",type=Path);p.add_argument("--reference",type=Path)
    p.add_argument("--base",type=Path,default=OLD,help="已登记的第三档网格或同一输入的重建目录")
    a=p.parse_args();base=a.base.resolve();base.relative_to(ROOT/"outputs");out=a.output.resolve();out.relative_to(BASE)
    key=read(ROOT/"cases/thermistor_b/leadsThird3.json")
    assert all(sha(base/k)==key[v] for k,v in KEYS.items())
    if (base/"engine_process.json").exists():terminal(base)
    if a.action=="baseline":
        s=metrics(base,out);print(s["closestEdgeMeanDeg"],s["edgeFractionWithin15Deg"],s["intervalMedianD"],s["intervalP90D"])
    else:
        if a.candidate is None or a.reference is None:p.error("compare 需要 --candidate 和 --reference")
        report=compare(base,a.candidate.resolve(),a.reference.resolve(),out)
        if not report["passed"]:raise SystemExit(2)
