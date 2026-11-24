"""独立从两套原面归档比较同工况、同物理时间与同一步长的数值外域误差。"""
import argparse
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile
import json

from outer_grid import ROOT, HERE, read, sha
import numpy as np
sys.path.insert(0,str(HERE.parent/"v333"))
from verify_time import table, verify

CASE=ROOT/"cases/short_sphere"


def snapshot(data,branch,pair):
    path=data/("supportedTime65"+branch.title()+".zip")
    assert sha(path)==pair[branch+"Archive"]
    with ZipFile(path) as archive:
        s=json.loads(archive.read("summary.json"))
        w=table(archive,"wall.csv");stations=table(archive,"stations.csv");direction=table(archive,"direction.csv")
    assert s["frameSHA256"]==pair[branch+"Frame"] and s["meshSHA256"]==pair["meshSHA256"]
    np.testing.assert_array_equal(w["faceID"],np.arange(len(w))+1)
    np.testing.assert_array_equal(direction["faceID"],w["faceID"])
    return s,w,stations,direction["status"]


def compare(data):
    base=read(CASE/"supportedTime65.json");outer=read(data/"supportedTime65.json")
    assert sha(CASE/"supportedTime65.json")=="42c79caa7d36d7e142407945919252bc826c526d8481f9a6b53dc47a71fdb499"
    for k in ("condition","faces","areaM2","caseSHA256","bodySHA256","transferSHA256","angles","experimentNu","elapsedS","elapsedDU"):
        assert base[k]==outer[k],k
    for k in ("meshSHA256","seed","steadyFrame"):assert base[k]!=outer[k]
    result=dict(faces=base["faces"],areaM2=base["areaM2"],elapsedS=base["elapsedS"],baseMesh=base["meshSHA256"],
                outerMesh=outer["meshSHA256"],jointExperimentalAreaM2=0,angleDeg=base["angles"],measuredNu=base["experimentNu"])
    for branch in ("full","half"):
        a,x,sa,da=snapshot(CASE,branch,base);b,y,sb,db=snapshot(data,branch,outer)
        assert a["dtSeconds"]==b["dtSeconds"] and a["elapsedSeconds"]==b["elapsedSeconds"]
        for k in ("faceID","regionID","areaM2","wallInputK"):np.testing.assert_array_equal(x[k],y[k])
        for k in ("originalFaceID","reportedNu","angleDeg","wallK","leftNode","rightNode","rightWeight"):
            np.testing.assert_array_equal(sa[k],sb[k])
        area=x["areaM2"];region=x["regionID"]
        flux=x["fluxInWm2"];q=y["fluxInWm2"];regions=[]
        for k in range(1,6):
            sel=region==k
            w14=-float(area[sel]@flux[sel]);w20=-float(area[sel]@q[sel])
            np.testing.assert_allclose([w14,w20],[a["regions"][k-1]["outW"],b["regions"][k-1]["outW"]],rtol=0,atol=1e-10)
            regions.append(dict(region=k,faces=int(sel.sum()),areaM2=float(area[sel].sum()),baseW=w14,outerW=w20,deltaW=w20-w14))
        row=dict(baseFrame=a["frameSHA256"],outerFrame=b["frameSHA256"],dtS=a["dtSeconds"],regions=regions,
                 baseOutW=a["energy"]["wallW"],outerOutW=b["energy"]["wallW"],experimentalEquivalentW=a["observedFullSphereW"],
                 baseSilverEquivalentW=a["silverEquivalentFullSphereW"],outerSilverEquivalentW=b["silverEquivalentFullSphereW"],
                 baseRMSE=a["stationRMSE"],outerRMSE=b["stationRMSE"],baseNu=sa["computedNu"].tolist(),outerNu=sb["computedNu"].tolist(),
                 baseRelativeError=a["equivalentRelativeError"],outerRelativeError=b["equivalentRelativeError"],
                 basePeakFace=a["peakOutFace"],outerPeakFace=b["peakOutFace"],basePeakWm2=a["peakOutWm2"],outerPeakWm2=b["peakOutWm2"],
                 baseReverseFaces=int((da==-1).sum()),outerReverseFaces=int((db==-1).sum()),
                 reverseChangedAreaM2=float(area[(da==-1)!=(db==-1)].sum()))
        for key in ("fluxInWm2","gaugePa","firstK"):
            delta=y[key]-x[key];peak=int(np.argmax(abs(delta)))
            row[key]=dict(rms=float(np.sqrt(np.average(delta**2,weights=area))),maxAbs=float(abs(delta[peak])),face=peak+1)
        result[branch]=row
    return result


def package(work,out):
    out.mkdir(parents=True,exist_ok=False);pack=out/"outer20";pack.mkdir()
    ack=read(work/"completed.json");pair=work/"export/pair.json";assert sha(pair)==ack["pairSHA256"]
    shutil.copyfile(pair,pack/"supportedTime65.json")
    for branch in ("full","half"):
        shutil.copyfile(work/"export"/branch/"surface.zip",pack/("supportedTime65"+branch.title()+".zip"))
    c=read(work/"steady/converged.json");assert c["frame"]==read(work/"steady_complete.json")["frame"]
    source=work/"steady"/c["stage"]/"restart_000199.csv"
    checks=verify(pack,work/"time",base_frame=c["frame"],steady=source)
    for name,value in (("independent",checks),("comparison",compare(pack))):
        (out/(name+".json")).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf8")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("out",type=Path)
    args=p.parse_args();package(args.work.resolve(),args.out.resolve())
