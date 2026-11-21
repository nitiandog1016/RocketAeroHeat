"""从真实网边独立交叉核对十六条射线，归档前中段的原面号、几何读数和局部网图。"""
import argparse
import csv
import io
from pathlib import Path
import shutil
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from audit_alignment import BASE,HERE,OLD,ROOT,LOW,HIGH,START,STOP,load_mesh,read,save,sha
from thermistorBead import D
from thermistorLeadMesh import original,meridian


def independent_cuts(folder,stats):
    """按端点相对射线的有符号切向坐标线性插值，独立复算每一个穿越位置。"""
    _,xyz,marks=load_mesh(folder/"mesh.su2")
    pairs=sorted({tuple(sorted((a,b))) for face in marks["sym0"] for a,b in zip(face,face[1:]+face[:1])})
    wall,_=meridian();edges=xyz[np.asarray(pairs),:2];results=[]
    assert [r["strip"] for r in stats["rays"]]==list(range(START,STOP))
    for row in stats["rays"]:
        j=row["strip"];t=wall[j+1]-wall[j];t/=np.linalg.norm(t);n=np.array([-t[1],t[0]])
        origin=(wall[j]+wall[j+1])/2
        np.testing.assert_array_equal(origin,row["originD"]);np.testing.assert_array_equal(n,row["normal"])
        s=(edges-origin)@t;h=(edges-origin)@n
        mask=(s[:,0]*s[:,1]<=0)&(abs(s[:,1]-s[:,0])>1e-14)
        u=-s[mask,0]/(s[mask,1]-s[mask,0])
        hits=h[mask,0]+u*(h[mask,1]-h[mask,0])
        hsort=np.sort(hits[(hits>LOW+1e-9)&(hits<HIGH-1e-9)])
        cuts=np.r_[LOW,hsort,HIGH];cuts=cuts[np.r_[True,np.diff(cuts)>1e-9]]
        gaps=np.diff(cuts);np.testing.assert_allclose(gaps,row["intervalsD"],rtol=0,atol=1e-11)
        assert len(gaps)==row["count"] and abs(gaps.max()-row["maximumD"])<1e-11
        k=int(np.argmax(gaps));results.append(dict(strip=j,xM=float(origin[0]*D),maximumDepthRangeD=cuts[k:k+2].tolist()))
    return edges,results


def plot(edges,stats,details,path):
    delta=np.array([b["maximumD"]/a["maximumD"]-1 for a,b in zip(stats[0]["rays"],stats[1]["rays"])])
    bad=int(np.argmax(abs(delta)));row=stats[0]["rays"][bad]
    origin=np.asarray(row["originD"]);n=np.asarray(row["normal"]);t=np.array([n[1],-n[0]])
    plt.rcParams.update({"font.family":"Microsoft YaHei","font.size":9,"axes.unicode_minus":False})
    fig,(left,mid,right)=plt.subplots(1,3,figsize=(11.2,4.5),gridspec_kw={"width_ratios":[.85,.85,1.15]})
    fig.subplots_adjust(left=.07,right=.98,top=.80,bottom=.19,wspace=.39)
    for i,(ax,lines,title) in enumerate(zip((left,mid),edges,("第三档基线","切向重排候选"))):
        coords=np.stack(((lines-origin)@t,(lines-origin)@n),axis=-1)*D*1e6
        lo=LOW*D*1e6;hi=HIGH*D*1e6
        visible=(coords[:,:,0].min(axis=1)<19)&(coords[:,:,0].max(axis=1)>-19)&(coords[:,:,1].min(axis=1)<hi)&(coords[:,:,1].max(axis=1)>lo)
        ax.add_collection(LineCollection(coords[visible],colors="#718795",linewidths=.65))
        ax.axvline(0,color="#404852",ls="--",lw=.8)
        ends=np.array(details[i][bad]["maximumDepthRangeD"])*D*1e6
        ax.plot([0,0],ends,color="#a23c20",lw=3,marker="o",ms=3)
        ax.set(xlim=(-19,19),ylim=(lo,hi),aspect="equal",xlabel="局部切向 [µm]",title=title+f"\n最宽间距 {np.diff(ends)[0]:.2f} µm")
        if i==0:ax.set_ylabel("距壁法向坐标 [µm]")
    inds=np.arange(START,STOP);right.axhspan(-25,25,color="#e4eced")
    right.axhline(25,color="#697a7a",lw=.8,ls="--");right.axhline(-25,color="#697a7a",lw=.8,ls="--")
    right.axhline(0,color=".6",lw=.7);right.plot(inds,100*delta,"o-",lw=1.2,ms=3,color="#2d617b")
    right.plot(inds[bad],100*delta[bad],"o",color="#a23c20",ms=5)
    right.annotate(f"{100*delta[bad]:+.2f}%",(inds[bad],100*delta[bad]),xytext=(4,6),textcoords="offset points")
    right.set(xlabel="原母线条带编号（从 0 起）",ylabel="最宽法向间距变化 [%]",title="十六条固定法线",ylim=(-32,37))
    right.set_xticks([54,58,62,66,69])
    fig.suptitle(f"热敏珠前中段网格　·　左两图为 {row['strip']} 号条带同一条法线",fontsize=14)
    fig.savefig(path,dpi=165,bbox_inches="tight");plt.close(fig)


def publish(base,candidate,reference,audit,out):
    result=read(audit/"geometry.json");placement=read(candidate/"placement.json")
    samples=(reference,audit);old,new=[read(d/"metrics.json") for d in samples]
    assert result["baseMesh"]==sha(base/"mesh.su2")==old["meshSHA256"]
    assert result["candidateMesh"]==sha(candidate/"mesh.su2")==new["meshSHA256"]==placement["meshSHA256"]
    assert placement["nodeMapSHA256"]==sha(candidate/"nodes.npz")
    assert result["baseMetricsSHA256"]==sha(reference/"metrics.json") and result["candidateMetricsSHA256"]==sha(audit/"metrics.json")
    assert result["planSHA256"]==sha(HERE/"PLAN.md") and result["auditSHA256"]==sha(HERE/"audit_alignment.py")
    assert placement["generatorSHA256"]==sha(ROOT/"scripts/thermistorFrontAlignment.py")
    for folder,s in zip(samples,(old,new)):
        assert sha(folder/"measurements.npz")==s["measurementSHA256"]
        with np.load(folder/"measurements.npz") as m:
            assert len(m["edgeNodeIDs"])==s["targetEdges"] and len(m["triangleNodeIDs"])==s["targetTriangles"]
            assert np.mean(m["edgeAngleDeg"]<=15)==s["edgeFractionWithin15Deg"]
            assert np.average(m["closestEdgeAngleDeg"],weights=m["triangleAreaD2"])==s["closestEdgeMeanDeg"]
    edges,details=zip(*(independent_cuts(f,s) for f,s in ((base,old),(candidate,new))))
    flags=result["checks"]
    assert not result["passed"] and {k for k,v in flags.items() if not v}=={"direction","rayMaxima"}
    failed=[r["strip"] for r,c in zip(old["rays"],result["maxIntervalRatios"]) if not .75<=c<=1.25]
    assert failed==[62]
    body,links,areas=original();strip=links[:,:,1].min(axis=1);mask=(strip>=START)&(strip<STOP)
    native=read(base/"body.json");assert body==native and len(areas)==10272 and mask.sum()==1536
    # 压缩包带真实面号和全部子午边读数。工程根目录只保留小型复核资料。
    out.mkdir(parents=True,exist_ok=False)
    with zipfile.ZipFile(out/"measurements.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for name,folder in zip(("baseline","candidate"),samples):
            for file in ("measurements.npz","metrics.json"):z.write(folder/file,name+"/"+file)
        stream=io.StringIO(newline="");w=csv.writer(stream);w.writerow(("face_id","strip_0","area_m2"))
        for i in np.flatnonzero(mask):w.writerow((i+1,strip[i],areas[i]))
        z.writestr("original_front_mid_faces.csv",stream.getvalue())
    rows=[]
    for j,(a,b) in enumerate(zip(old["rays"],new["rays"])):
        rows.append(dict(strip=a["strip"],xM=details[0][j]["xM"],originalFaces=int((strip==a["strip"]).sum()),
                         originalAreaM2=float(areas[strip==a["strip"]].sum()),baselineCount=a["count"],candidateCount=b["count"],
                         baselineMaxD=a["maximumD"],candidateMaxD=b["maximumD"],widthRatio=b["maximumD"]/a["maximumD"],
                         baselineMaxDepthD=details[0][j]["maximumDepthRangeD"],candidateMaxDepthD=details[1][j]["maximumDepthRangeD"]))
    save(out/"rays.json",dict(rows=rows,failedMaximumAtStrip=failed,originalFaces=int(mask.sum()),originalAreaM2=float(areas[mask].sum()),
         allFaces=len(areas),allAreaM2=float(areas.sum()),sourceBodySHA256=sha(base/"body.json"),geometryAccepted=False,
         measurementArchiveSHA256=sha(out/"measurements.zip")))
    for p in (candidate/"placement.json",audit/"geometry.json"):shutil.copyfile(p,out/p.name)
    plot(edges,(old,new),details,out/"alignment.png")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("candidate","reference","audit","output"):p.add_argument(name,type=Path)
    p.add_argument("--base",type=Path,default=OLD)
    args=p.parse_args();base=args.base.resolve();base.relative_to(ROOT/"outputs")
    paths=[getattr(args,k).resolve() for k in ("candidate","reference","audit","output")]
    for path in paths:path.relative_to(BASE)
    publish(base,*paths)
