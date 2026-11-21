"""在相同原面上分别按给定壁温和给定总功率比较两种外域。"""
import argparse
import zipfile
from pathlib import Path

import numpy as np

from check_state import BASE,NEW,ROOT,read,save,sha,terminal


def load(folder):
    rec=terminal(folder);s=read(folder/"surface/summary.json");flow=read(folder/"surface/flow.json")
    assert s["sourceFrameSHA256"]==flow["sourceFrameSHA256"]==sha(folder/"restart.csv")
    assert s["sourceReactionSHA256"]==sha(folder/"reaction0/restart.csv")
    assert s["meshSHA256"]==sha(folder/"mesh.su2") and s["caseSHA256"]==sha(folder/"case.json")
    assert s["wallSHA256"]==sha(folder/"surface/wall.csv") and flow["tableSHA256"]==sha(folder/"surface/flow.csv")
    assert s["energy"]["sampledThetaExact"] and s["relativeChange20"]<.002 and abs(s["energy"]["relativeDifference"])<.005
    assert max(s["residualWindow"].values())< -8 and max(s["plus20Window"].values())< -8
    w=np.genfromtxt(str(folder/"surface/wall.csv"),names=True,delimiter=",")
    reverse=np.genfromtxt(str(folder/"surface/flow.csv"),names=True,delimiter=",")
    assert np.array_equal(reverse["face_id"],w["faceID"]+1)
    assert abs(-w["areaM2"]@w["heatFluxWm2"]-s["totalOutwardW"])<1e-13
    return s,w,reverse["reverse_both"].astype(bool)


def metrics(values,area):
    j=int(np.argmax(abs(values)))
    return dict(maxAbs=float(abs(values[j])),signedAtMax=float(values[j]),faceID=j+1,rms=float(np.sqrt(np.average(values**2,weights=area))))


def compare(out):
    a,b,r0=load(BASE);s,w,r=load(NEW)
    assert a==read(ROOT/"docs/validation/v300/evidence/summary.json")
    assert a["bodySHA256"]==s["bodySHA256"]==sha(NEW/"body.json")
    c=read(NEW/"case.json");assert c["outerRadiusOverD"]==30
    assert dict(c,outerRadiusOverD=20)==read(BASE/"case.json")
    assert np.array_equal(w["faceID"],np.arange(10272)) and np.array_equal(w["faceID"],b["faceID"])
    groups=w["regionID"];area=w["areaM2"]
    assert np.array_equal(area,b["areaM2"]) and np.array_equal(groups,b["regionID"])
    out.mkdir(parents=True,exist_ok=False)
    summary=dict(baseFrame=a["sourceFrameSHA256"],expandedFrame=s["sourceFrameSHA256"],baseMesh=a["meshSHA256"],expandedMesh=s["meshSHA256"],
                 reverseChangedFaces=int(np.count_nonzero(r0!=r)),reverseChangedAreaM2=float(area[r0!=r].sum()),
                 absolutePressureChangePa=metrics(w["pressurePa"]-b["pressurePa"],area),
                 changeInConductanceFraction=s["totalOutwardW"]/a["totalOutwardW"]-1)
    for mode in ("givenWall","givenPower"):
        scales=[1.,1.] if mode=="givenWall" else [.004/t["totalOutwardW"] for t in (a,s)]
        q0,q=scales[0]*b["heatFluxWm2"],scales[1]*w["heatFluxWm2"]
        t0,t=[293.15+scale*(field["nearWallCellTemperatureK"]-293.15) for field,scale in zip((b,w),scales)]
        tw0,tw=[293.15+scale*12 for scale in scales]
        watts0,watts=float(-area@q0),float(-area@q)
        if mode=="givenPower":assert max(abs(watts0-.004),abs(watts-.004))<1e-13
        regions=[]
        for g,label in enumerate(read(NEW/"body.json")["regionNames"],1):
            mask=groups==g;j0=int(np.flatnonzero(mask)[np.argmin(q0[mask])]);j=int(np.flatnonzero(mask)[np.argmin(q[mask])])
            old=float(-area[mask]@q0[mask]);new=float(-area[mask]@q[mask])
            regions.append(dict(name=label,faces=int(mask.sum()),areaM2=float(area[mask].sum()),oldW=old,newW=new,deltaW=new-old,
                                peak20Face=j0+1,peak30Face=j+1,peak20OutWm2=float(-q0[j0]),peak30OutWm2=float(-q[j])))
        summary[mode]=dict(wall20K=tw0,wall30K=tw,deltaWallK=tw-tw0,outward20W=watts0,outward30W=watts,regions=regions,
            deltaOutwardW=watts-watts0,fluxChangeWm2=metrics(q-q0,area),maxRelativeFluxChange=float(np.max(abs(q/q0-1))),
            firstTemperatureChangeK=metrics(t-t0,area),peak20Face=int(np.argmin(q0))+1,peak30Face=int(np.argmin(q))+1)
        np.testing.assert_allclose(sum(k["deltaW"] for k in regions),watts-watts0,atol=1e-13,rtol=0)
    save(out/"comparison.json",summary)
    # 30D 原面积分与温度表是工作台的独立对照，沿用求解器产出的原始归档。
    casefile=NEW/"surface/surface.zip"
    with zipfile.ZipFile(casefile) as z:
        assert z.testzip() is None
        for name in ("body.json","case.json"):assert z.read(name)==(NEW/name).read_bytes()
        for name in ("wall.csv","summary.json"):assert z.read(name)==(NEW/"surface"/name).read_bytes()
    save(out/"outer30.json",dict(archiveSHA256=sha(casefile),baseMeshSHA256=a["meshSHA256"],baseCaseSHA256=a["caseSHA256"],
         bodySHA256=s["bodySHA256"],meshSHA256=s["meshSHA256"],caseSHA256=s["caseSHA256"],sourceFrameSHA256=s["sourceFrameSHA256"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("out",type=Path)
    dest=p.parse_args().out.resolve();dest.relative_to(ROOT/"outputs");compare(dest)
