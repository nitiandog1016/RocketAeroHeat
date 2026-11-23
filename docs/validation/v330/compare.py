"""在同一组原面上对照粗细网，分开计算面积加权差异、热峰区域和试验偏差。"""
import argparse
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

from check_grid import HERE, ROOT, sha

DATA=ROOT/"cases/short_sphere"


def load_surface(path):
    with ZipFile(path) as z:
        report=json.loads(z.read("summary.json"))
        w=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
        stations=np.genfromtxt(io.BytesIO(z.read("stations.csv")),names=True,delimiter=",")
    assert report["thermal"]["passed"] and report["faces"]==len(w)==32480
    assert abs(w["areaM2"]@w["fluxInWm2"]+report["thermal"]["early"]["wallW"])<1e-10
    return report,w,stations


def differences(values,areas):
    order=np.argsort(values,kind="stable")
    weight=np.cumsum(areas[order]);total=float(weight[-1])
    idx=np.searchsorted(weight,total*np.array([.5,.95,.99]))
    return dict(areaMean=float(areas@values/areas.sum()),p50=float(values[order[idx[0]]]),
                p95=float(values[order[idx[1]]]),p99=float(values[order[idx[2]]]),maximum=float(max(values)))


def compare(fine,grid):
    coarse=DATA/"supportedSurface65.zip"
    assert sha(coarse)=="2e5b065e3b7bd7f15dc46a22bdc6f5a041bf57ebb9a158c31cb46426953c772b"
    c,w,cs=load_surface(coarse);f,v,fs=load_surface(fine)
    for key in ("condition","bodySHA256","caseSHA256","transferSHA256","areaM2"):
        assert c[key]==f[key]
    assert f["thermal"]["meshSHA256"]==grid["meshSHA256"]!=c["thermal"]["meshSHA256"]
    assert f["frameSHA256"]!=c["frameSHA256"]
    for name in ("faceID","areaM2","regionID","wallInputK"):np.testing.assert_array_equal(w[name],v[name])
    np.testing.assert_array_equal(cs["angleDeg"],fs["angleDeg"]);np.testing.assert_array_equal(cs["reportedNu"],fs["reportedNu"])
    a=w["areaM2"];q=-w["fluxInWm2"];fq=-v["fluxInWm2"];base=q>=.99*q.max();changed=fq>=.99*fq.max()
    regions=[]
    for cr,fr in zip(c["regions"],f["regions"]):
        k=cr["id"];mask=w["regionID"]==k;assert fr["id"]==k and cr["faces"]==fr["faces"]==mask.sum()
        regions.append(dict(id=k,faces=cr["faces"],areaM2=cr["areaM2"],coarseW=cr["outW"],fineW=fr["outW"],
            relativePowerChange=fr["outW"]/cr["outW"]-1,coarsePeakWm2=q[mask].max(),finePeakWm2=fq[mask].max(),
            heatAbsDifference=differences(abs(q[mask]-fq[mask]),a[mask]),
            pressureAbsDifference=differences(abs(w["gaugePa"][mask]-v["gaugePa"][mask]),a[mask])))
    def station(s):return s["computedNu"].tolist()
    def limits(x,y):return [float(min(x.min(),y.min())),float(max(x.max(),y.max()))]
    return dict(condition=f["condition"],bodySHA256=f["bodySHA256"],caseSHA256=f["caseSHA256"],transferSHA256=f["transferSHA256"],
        coarseArchiveSHA256=sha(coarse),fineArchiveSHA256=sha(fine),coarseMesh=c["thermal"]["meshSHA256"],fineMesh=f["thermal"]["meshSHA256"],
        coarseFrame=c["frameSHA256"],fineFrame=f["frameSHA256"],faces=f["faces"],areaM2=f["areaM2"],jointVerifiedAreaM2=0,
        coarseW=c["thermal"]["early"]["wallW"],fineW=f["thermal"]["early"]["wallW"],
        heatAbsDifference=differences(abs(q-fq),a),pressureAbsDifference=differences(abs(w["gaugePa"]-v["gaugePa"]),a),
        regions=regions,anglesDeg=cs["angleDeg"].tolist(),measuredNu=cs["reportedNu"].tolist(),coarseNu=station(cs),fineNu=station(fs),
        coarseStationRMSE=c["stationRMSE"],fineStationRMSE=f["stationRMSE"],coarsePowerError=c["equivalentRelativeError"],finePowerError=f["equivalentRelativeError"],
        coarsePeak=dict(faceID=int(q.argmax()+1),peakWm2=q.max(),nearFaces=w["faceID"][base].astype(int).tolist(),nearAreaM2=a[base].sum()),
        finePeak=dict(faceID=int(fq.argmax()+1),peakWm2=fq.max(),nearFaces=v["faceID"][changed].astype(int).tolist(),nearAreaM2=a[changed].sum()),
        sharedNearAreaM2=a[base&changed].sum(),unionNearAreaM2=a[base|changed].sum(),
        sharedFluxRangeKWm2=limits(-q/1000,-fq/1000),sharedGaugeRangePa=limits(w["gaugePa"],v["gaugePa"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("fine",type=Path);p.add_argument("grid",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();result=compare(a.fine,json.loads(a.grid.read_bytes()))
    with a.out.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print(result["coarseW"],result["fineW"],result["finePowerError"])
