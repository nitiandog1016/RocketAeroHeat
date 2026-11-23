"""从原生逐边数据、独立释放反力及流动边界通量重建双支管银球热量账。"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from run import BASE, CASE, HERE, PREP, ROOT, read, sha
from thermistorWireSurface import residual


def terminal(folder):
    proc=read(folder/"debugger_process.json")
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")
    return proc


def table(folder,name):
    return np.loadtxt(folder/name,delimiter=",",skiprows=1,ndmin=2)


def history(stage):
    terminal(stage)
    rec=read(stage/"stage.json")
    assert sha(stage/"history.csv")==rec["history"] and rec["trace"]==sha(stage/"stages.csv")
    h=table(stage,"history.csv")
    np.testing.assert_array_equal(h[:,2],np.arange(200))
    with (stage/"stages.csv").open(encoding="ascii") as f:rows=list(csv.DictReader(f))
    assert len(rows)==400
    for i,r in enumerate(rows):
        assert int(r["iteration"])==i//2 and r["stage"]==("actual" if i%2 else "low")
        assert int(r["above"])==int(r["below"])==0 and float(r["wallDrift"])<2e-14
        assert abs(float(r["balance"])) < 2e-12*max(1,float(r["scale"])) and float(r["rowError"])<2e-13
    raw=np.array([float(r["rawRMS"]) for r in rows[1::2]])
    limited=np.array([float(r["schemeRMS"]) for r in rows[1::2]])
    np.testing.assert_allclose(np.log10(limited),h[:,8],rtol=0,atol=2e-8)
    np.testing.assert_array_equal(rec["window"],h[-20:,3:9].max(axis=0))
    for it in (100,199):
        p=stage/f"restart_{it:06d}.csv";assert sha(p)==rec["frames"][p.name]
        scalar=np.loadtxt(p,delimiter=",",skiprows=1,usecols=9)
        np.testing.assert_allclose([scalar.min(),scalar.max()],[float(rows[2*it+1][k]) for k in ("min","max")],rtol=0,atol=1e-15)
    return h,raw,rows


def budget(folder,release,wet,scale):
    terminal(folder);terminal(release)
    n=np.fromfile(folder/"nodes.bin",dtype="<f8").reshape(-1,16)
    n=n[np.argsort(n[:,0])]
    np.testing.assert_array_equal(n[:,0],np.arange(len(n)))
    old,low,final=n[:,1:4].T
    np.testing.assert_array_equal(old,table(folder,"initial_flow.csv")[:,9])
    np.testing.assert_array_equal(final,table(folder,"restart.csv")[:,9])
    wall=n[:,4].astype(bool);np.testing.assert_array_equal(np.flatnonzero(wall),wet)
    assert np.isfinite(n).all() and min(low.min(),final.min())>=-1e-12
    upper=old[wall].max();assert max(low.max(),final.max())<=upper+1e-12
    e=np.fromfile(folder/"edges.bin",dtype="<f8").reshape(-1,6)
    assert np.isfinite(e).all()
    i,j=e[:,:2].astype(int).T;mass,high,corr,alpha=e[:,2:].T
    def accumulate(left,right):
        return np.bincount(i,left,minlength=len(n))+np.bincount(j,right,minlength=len(n))
    def net(face):
        return accumulate(face,-face)
    div=-old*(net(mass)+n[:,11])
    raw=net(high)+div+n[:,12]
    np.testing.assert_allclose(div,n[:,13],atol=2e-14,rtol=1e-12)
    np.testing.assert_allclose(raw,n[:,8],atol=2e-14,rtol=1e-12)
    np.testing.assert_allclose(raw[~wall],n[~wall,7],atol=2e-14,rtol=1e-12)
    # 依据前态、低阶预测、局部体积和伪步长重新决定每条边的共享系数。
    outward=accumulate(np.maximum(0,corr),np.maximum(0,-corr))
    inward=accumulate(np.maximum(0,-corr),np.maximum(0,corr))
    coeff=[]
    for amount,room in ((outward,np.minimum(old,low)),(inward,upper-np.maximum(old,low))):
        a=np.ones(len(n));mask=(amount>0)&~wall
        a[mask]=np.minimum(1,(n[mask,5]/n[mask,6])*np.maximum(0,room[mask])/amount[mask])
        coeff.append(a)
    give,take=coeff
    checked=np.where(corr>=0,np.minimum(give[i],take[j]),np.minimum(take[i],give[j]))
    np.testing.assert_allclose(checked,alpha,atol=1e-12,rtol=1e-12)
    applied=net(checked*corr)
    np.testing.assert_allclose(applied,n[:,10],atol=1e-16,rtol=1e-12)
    unlimited_release=residual(release/"restart.csv")
    assert sha(folder/"initial_flow.csv")==sha(release/"initial_flow.csv")
    releasedError=max(abs(unlimited_release[wet]+raw[wet]))*scale
    assert releasedError<1e-10
    corrected=raw-net(corr)+applied
    wallW=float(np.sum(-unlimited_release[wet]-net(corr)[wet]+applied[wet])*scale)
    np.testing.assert_allclose(wallW,corrected[wet].sum()*scale,rtol=0,atol=1e-10)
    far=np.loadtxt(folder/"flow_far.csv",delimiter=",");ids=far[:,0].astype(int);fm=far[:,1]
    assert len(ids)==len(np.unique(ids))==1346
    assert not wall[ids].any()
    # 给定远场标量为零。独立上风读数分别使用两种实际质量通量。
    np.testing.assert_allclose(n[ids,12],np.maximum(n[ids,11],0)*old[ids],rtol=1e-13,atol=1e-20)
    adv=float(n[:,12].sum()*scale);flowAdv=float(np.maximum(fm,0)@old[ids]*scale)
    divergence=float(div.sum()*scale);free=float(corrected[~wall].sum()*scale)
    err=wallW+free-adv-divergence
    assert abs(err)<1e-10
    farMask=np.zeros(len(n),bool);farMask[ids]=True
    return dict(wallW=wallW,scalarFarW=adv,flowFluxComparisonW=flowAdv,physicalGapW=adv-wallW,
                physicalGapFraction=(adv-wallW)/wallW,freeSpatialW=free,massDivergenceW=divergence,
                divergenceWallW=float(div[wall].sum()*scale),divergenceFarW=float(div[farMask].sum()*scale),
                divergenceInteriorW=float(div[~(wall|farMask)].sum()*scale),assemblyGapW=err,
                releaseMaxNodeErrorW=releasedError,farSignDifferences=int(np.sum(np.sign(n[ids,11])!=np.sign(fm))),
                edges=len(e),farNodes=len(ids),limitedEdges=int(np.sum((alpha<1)&(corr!=0))),
                sourceFrame=sha(folder/"initial_flow.csv"),observerSHA256=sha(folder/"restart.csv"),
                binaryFiles={x:sha(folder/x) for x in ("nodes.bin","edges.bin","flow_far.csv")})


def verify():
    s=read(PREP/"prepared.json")["condition"];affine=BASE/"build01/affine.log"
    with affine.open() as f:results=list(csv.reader(f))
    for n,(count,q,exact,alpha) in zip((3,9,27),results[:3]):
        assert int(count)==n and float(alpha)==1 and abs(float(q)-.05/s["pe"])<1e-14 and abs(float(exact)-float(q))<1e-14
    np.testing.assert_allclose([float(x) for x in results[3][1:]],[.01,.06,.015],rtol=0,atol=1e-14)
    stages=[BASE/"run01/steady01",BASE/"run01/steady02",BASE/"repeat03"]
    rejected=BASE/"run01/steady03";r=read(rejected/"debugger_process.json")
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and [Path(x).name for x,ok in r["inputsStable"].items() if not ok]==["PLAN.md"]
    repeatDiff={}
    for name in ("history.csv","restart_000100.csv","restart_000199.csv"):
        a=table(rejected,name);b=table(stages[-1],name)
        if name=="history.csv":a,b=a[:,2:9],b[:,2:9]
        repeatDiff[name]=float(np.max(abs(a-b)))
        assert repeatDiff[name]<1e-12
    summary=[];trend=[];seed=sha(PREP/"initial.csv")
    for index,stage in enumerate(stages):
        h,raw,rows=history(stage);rec=read(stage/"stage.json")
        assert rec["start"]==seed==sha(stage/"initial_flow.csv");seed=rec["frames"]["restart_000199.csv"]
        for k in [*range(0,200,5),199]:trend.append([200*index+k+1,*h[k,3:9],np.log10(raw[k])])
        summary.append(dict(frame=seed,window=rec["window"],rawMax20=raw[-20:].max(),seconds=rec["seconds"],
                            minAll=float(min(float(r["min"]) for r in rows)),maxAll=float(max(float(r["max"]) for r in rows)),
                            files=rec["frames"],peakBytes=rec["peakBytes"]))
    assert max(summary[-1]["window"])<-8 and summary[-1]["rawMax20"]<1e-8
    base=table(BASE/"run01/steady01","restart_000199.csv");orig=table(ROOT/"outputs/v326_validation/run01/steady01","restart_000199.csv")
    np.testing.assert_array_equal(base[:,:9],orig[:,:9])
    samples=BASE/"thermal01";receipts=read(samples/"samples.json")
    for name,rec in receipts.items():
        folder=samples/name;terminal(folder);assert sha(folder/"restart.csv")==rec["end"] and sha(folder/"initial_flow.csv")==rec["start"]
        assert rec["start"]==(receipts["plus20"]["end"] if name in ("late","releaseLate") else seed)
    a=table(samples/"early","restart.csv");b=table(samples/"control","restart.csv")
    observerError=float(abs(a-b).max());assert observerError<2e-12
    h20=table(samples/"plus20","history.csv");assert len(h20)==20 and np.max(h20[:,3:9])<-8
    with np.load(PREP/"mapping.npz") as m:wet=m["wet"]
    scale=s["rho"]*s["speedMps"]*s["cp"]*s["diameterM"]**2*s["scalarUnitK"]
    early=budget(samples/"early",samples/"releaseEarly",wet,scale)
    late=budget(samples/"late",samples/"releaseLate",wet,scale)
    change=abs(late["wallW"]/early["wallW"]-1)
    balancePass=max(abs(e["physicalGapFraction"]) for e in (early,late))<.005
    return dict(caseSHA256=sha(CASE),meshSHA256=sha(PREP/"mesh.su2"),coldSHA256=sha(PREP/"initial.csv"),
                coefficientToW=scale,stages=summary,trend=trend,repeatDifferences=repeatDiff,
                updates=600,nodes=len(a),boundViolations=0,extra20Fraction=change,thermalStable=bool(change<.002),
                balancePass=balancePass,early=early,late=late,observerMaxDifference=observerError,
                originalFaces=32480,originalAreaM2=read(ROOT/"cases/short_sphere/supportedFlux65.json")["originalAreaM2"],
                acceptedFaces=0,acceptedAreaM2=0.)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--write",action="store_true");args=parser.parse_args()
    result=verify();file=ROOT/"cases/short_sphere/supportedCapacity65.json"
    if args.write:
        with file.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    assert read(file)==result
    print(sha(file),result["early"]["physicalGapFraction"],result["extra20Fraction"])
