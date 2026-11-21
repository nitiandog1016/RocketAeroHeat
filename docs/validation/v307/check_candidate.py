"""从唯一一张换边成品读取全部体连接和 48 条实际测线，判断能否启动冷态计算。"""
import argparse
from pathlib import Path
import sys
import zipfile
import json

from flip_bound import DEFAULT_BASE,HERE,OUTPUTS,ROOT,START,STOP,LOW,HIGH,meridian,load_mesh,read,save,sha,np
sys.path.insert(0,str(HERE.parent/"v306"))
from audit_alignment import metrics
from thermistorLeadMesh import SIDES
from thermistorOuterDomain import ordered


def topology(cells,marks):
    """按实际有序顶点核对每个边界面的唯一归属、所有内部面的双侧相反绕序。"""
    boundary={}
    for name,faces in marks.items():
        for face in faces:
            key=tuple(sorted(face));assert key not in boundary
            boundary[key]=ordered(face)
    inner={};closed=set()
    for cell in cells:
        for side in SIDES[len(cell)]:
            f=tuple(cell[j] for j in side);key=tuple(sorted(f))
            assert key not in closed
            if key in boundary:
                assert boundary.pop(key)==ordered(f);closed.add(key)
            elif key in inner:
                assert inner.pop(key)==ordered(f[::-1]);closed.add(key)
            else:inner[key]=ordered(f)
    assert not inner and not boundary
    return len(closed)-sum(map(len,marks.values()))


def volume_report(xyz,cells):
    """保留每个失败单元的原编号及面锥体积，同一张网其余条件仍可独立测量。"""
    total=0.;min_cell=float("inf");min_cone=float("inf");failures=[]
    for width,sides in SIDES.items():
        idx=np.array([j for j,c in enumerate(cells) if len(c)==width])
        if not len(idx):continue
        points=xyz[np.array([cells[i] for i in idx])];center=points.mean(axis=1)
        parts=[]
        for side in sides:
            f=points[:,side];n=np.cross(f[:,1]-f[:,0],f[:,2]-f[:,0])/2
            if len(side)==4:n+=np.cross(f[:,2]-f[:,0],f[:,3]-f[:,0])/2
            parts.append(np.einsum("ij,ij->i",n,f.mean(axis=1)-center)/3)
        parts=np.array(parts).T;vol=parts.sum(axis=1)
        total+=float(vol.sum());min_cone=min(min_cone,float(parts.min()));min_cell=min(min_cell,float(vol.min()))
        for j in np.flatnonzero((parts<=0).any(axis=1)):
            cell=int(idx[j]);failures.append(dict(cellZeroBased=cell,nodeIDs=list(map(int,cells[cell])),
                                                   volumeD3=float(vol[j]),faceConeVolumesD3=parts[j].tolist()))
    return dict(volumeD3=total,minCellVolumeD3=min_cell,minFaceConeD3=min_cone,nonpositiveCones=failures)


def cuts(xyz,marks,wall):
    pairs=sorted({tuple(sorted((i,j))) for f in marks["sym0"] for i,j in zip(f,f[1:]+f[:1])})
    segments=xyz[np.asarray(pairs),:2];rows=[]
    for strip in range(START,STOP):
        side=wall[strip+1]-wall[strip];t=side/np.linalg.norm(side);n=np.array([-t[1],t[0]])
        for fraction in (.25,.50,.75):
            origin=wall[strip]+fraction*side
            s=(segments-origin)@t;h=(segments-origin)@n
            eligible=(s[:,0]*s[:,1]<=0)&(abs(s[:,1]-s[:,0])>1e-14)
            u=-s[eligible,0]/(s[eligible,1]-s[eligible,0])
            v=h[eligible,0]+u*(h[eligible,1]-h[eligible,0])
            v=np.r_[LOW,np.sort(v[(v>LOW+1e-9)&(v<HIGH-1e-9)]),HIGH]
            v=v[np.r_[True,np.diff(v)>1e-9]]
            assert len(v)>2 and abs(np.diff(v).sum()-(HIGH-LOW))<1e-12
            rows.append(dict(strip=strip,fraction=fraction,cutsD=v.tolist()))
    return rows


def run(base,candidate,selection,out):
    sel=read(selection);assert sel["constructible"] and sel["candidateMeshSHA256"]==sha(candidate/"mesh.su2")
    assert sel["sourceMeshSHA256"]==sha(base/"mesh.su2")
    assert sha(base/"nodes.npz")==sha(candidate/"nodes.npz")
    old,xyz,om=load_mesh(base/"mesh.su2");new,coords,nm=load_mesh(candidate/"mesh.su2")
    np.testing.assert_array_equal(coords,xyz);assert len(old)==len(new) and om.keys()==nm.keys()
    for k in ("wall","far"):assert om[k]==nm[k]
    assert all(len(om[k])==len(nm[k]) for k in om)
    changes={i for i,(a,b) in enumerate(zip(old,new)) if a!=b}
    assert len(changes)==sel["replacedCells"]==8*len(sel["pairs"])
    for tag in ("sym0","sym1"):
        assert sum(a!=b for a,b in zip(om[tag],nm[tag]))==2*len(sel["pairs"])
    interior=topology(new,nm)
    quality=volume_report(coords,new);v0=read(base/"geometry.json")["volumeD3"]
    assert np.isclose(quality["volumeD3"],v0,rtol=2e-12,atol=0)
    wall=meridian()[0];a,b=cuts(xyz,om,wall),cuts(coords,nm,wall)
    out.mkdir(parents=True,exist_ok=False);m=metrics(candidate,out/"measured")
    with zipfile.ZipFile(HERE.parent/"v306/evidence/measurements.zip") as z:
        baseline=json.loads(z.read("baseline/metrics.json"))
    assert baseline["meshSHA256"]==sel["sourceMeshSHA256"] and m["meshSHA256"]==sel["candidateMeshSHA256"]
    widths0=[];widths1=[];rows=[];locked=[]
    for u,v in zip(a,b):
        x,y=np.array(u["cutsD"]),np.array(v["cutsD"]);d0,d1=np.diff(x),np.diff(y)
        assert (u["strip"],u["fraction"])==(v["strip"],v["fraction"])
        detail=dict(strip=u["strip"],fraction=u["fraction"],baseCount=len(d0),newCount=len(d1),countRatio=len(d1)/len(d0),
                    baseMaxD=float(d0.max()),newMaxD=float(d1.max()),widthRatio=float(d1.max()/d0.max()),baseCutsD=x.tolist(),newCutsD=y.tolist())
        if u["fraction"]==.5:
            np.testing.assert_allclose(d0,baseline["rays"][u["strip"]-START]["intervalsD"],atol=1e-11,rtol=0)
            np.testing.assert_allclose(d1,m["rays"][u["strip"]-START]["intervalsD"],atol=1e-11,rtol=0)
            detail["maxCutDisplacementD"]=float(np.max(abs(x-y))) if len(x)==len(y) else None
            locked.append(len(x)==len(y) and np.max(abs(x-y))<=1e-9)
        else:
            widths0.extend(d0);widths1.extend(d1)
        rows.append(detail)
    hold=[r for r in rows if r["fraction"]!=.5]
    p0=np.percentile(widths0,[50,90]);p1=np.percentile(widths1,[50,90])
    reduction=1-m["closestEdgeMeanDeg"]/baseline["closestEdgeMeanDeg"]
    assert abs(reduction-sel["predictedFraction"])<1e-11 and abs(m["sampleAreaD2"]-baseline["sampleAreaD2"])<1e-12
    checks=dict(positiveFaceCones=not quality["nonpositiveCones"],direction=bool(reduction>=.2),nearNormalEdges=m["edgeFractionWithin15Deg"]>baseline["edgeFractionWithin15Deg"],
                allOriginalCuts=all(locked),quarterMedian=bool(.85<=p1[0]/p0[0]<=1.15),quarterP90=bool(.85<=p1[1]/p0[1]<=1.15),
                quarterCounts=all(.8<=r["countRatio"]<=1.2 for r in hold),quarterMaxima=all(.75<=r["widthRatio"]<=1.25 for r in hold))
    save(out/"rays.json",dict(rows=rows))
    result=dict(passed=all(checks.values()),checks=checks,baseMeanDeg=baseline["closestEdgeMeanDeg"],newMeanDeg=m["closestEdgeMeanDeg"],
                reduction=reduction,baseFractionWithin15=baseline["edgeFractionWithin15Deg"],newFractionWithin15=m["edgeFractionWithin15Deg"],
                targetTriangles=m["targetTriangles"],baseP50P90D=p0.tolist(),newP50P90D=p1.tolist(),interiorFacesPaired=interior,
                unchangedCells=len(old)-len(changes),replacedCells=len(changes),cells=len(old),nodes=len(coords),quality=quality,
                failedCountRays=[(r["strip"],r["fraction"]) for r in hold if not .8<=r["countRatio"]<=1.2],
                failedWidthRays=[(r["strip"],r["fraction"]) for r in hold if not .75<=r["widthRatio"]<=1.25],
                sourceMeshSHA256=sel["sourceMeshSHA256"],candidateMeshSHA256=m["meshSHA256"],selectionSHA256=sha(selection),
                raysSHA256=sha(out/"rays.json"),metricsSHA256=sha(out/"measured/metrics.json"),auditSHA256=sha(Path(__file__).resolve()))
    save(out/"geometry.json",result)
    print(json.dumps(result,indent=2))
    return result["passed"]


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("candidate",type=Path);p.add_argument("selection",type=Path)
    p.add_argument("output",type=Path);p.add_argument("--base",type=Path,default=DEFAULT_BASE);a=p.parse_args()
    out=a.output.resolve();out.relative_to(OUTPUTS)
    if not run(a.base.resolve(),a.candidate.resolve(),a.selection.resolve(),out):raise SystemExit(2)
