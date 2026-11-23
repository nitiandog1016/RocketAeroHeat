"""从同步终帧和壁反力独立复算修正工况的全场，并打包可在侧栏读取的数值核对。"""
import argparse
import csv
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v314")]
from audit_angular import independent_overlap
from thermistorBead import read,sha


def package(base,work,graph,out):
    request=read(work/"request.json");f=request["model"]["factor"];c=read(work/"prepared/case.json")
    reading=read(graph/"reading.json");assert abs(reading["reynoldsL"]-request["model"]["reynoldsL"])<1e-10
    source=base/read(base/"completed.json")["stage"];stage=work/read(work/"completed.json")["stage"]
    ref=read(base/"surface/summary.json");s=read(work/"surface/summary.json")
    assert ref["sourceFrameSHA256"]==request["wireFrame"]==sha(source/"restart.csv")
    assert s["wallSHA256"]==sha(work/"surface/wall.csv") and s["sourceFrameSHA256"]==sha(stage/"restart.csv")
    assert sha(source/"plus20/initial_flow.csv")==request["wireFrame"] and sha(source/"plus20/restart.csv")==ref["plus20FrameSHA256"]
    assert sha(source/"reaction20/restart.csv")==ref["reaction20SHA256"]
    assert sha(stage/"reaction0/initial_flow.csv")==s["sourceFrameSHA256"] and sha(stage/"reaction0/restart.csv")==s["reactionSHA256"]
    assert sha(stage/"wire_boundary.dat")==s["boundarySHA256"]==sha(work/"prepared/profile.dat")
    assert ref["meshSHA256"]==s["meshSHA256"] and ref["nodeMapSHA256"]==s["nodeMapSHA256"]
    for dest in (stage,stage/"reaction0",stage/"plus20",stage/"reaction20",stage/"step_control",stage/"step_observed",source/"plus20",source/"reaction20"):
        r=read(dest/"debugger_process.json");q=read(dest/"native.json")
        assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
        assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
        assert sha(dest/"case.cfg")==r["inputSHA256"][str(dest/"case.cfg")]
        assert sha(dest/"history.csv")==q["outputs"]["history.csv"] and sha(dest/"restart.csv")==q["outputs"]["restart.csv"]
    assert read(stage/"native.json")["steps"]==read(source/"plus20/native.json")["steps"]==20
    assert max(s["window"].values())<-8 and max(s["plus20Window"].values())<-8 and abs(s["energy"]["relativeDifference"])<.005 and s["relativeChange20"]<.002
    with np.load(work/"prepared/nodes.npz") as z:xyz,wall,first,pi,ci=(z[k] for k in ("xyz","wallMap","originalFirstMap","overlapParent","overlapChild"))
    body=read(work/"prepared/body.json");tri=np.array(body["V"])[np.array(body["F"])-1];group=np.array(body["regionID"])
    small=xyz[wall]*c["diameterM"]
    def areas(t):return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    area,child_area=areas(tri),areas(small);wet=np.unique(wall)
    weights=independent_overlap(tri,small,pi,ci)*child_area[ci,None];lump=np.bincount(wall.ravel(),np.repeat(child_area/3,3),minlength=len(xyz))
    def project(nodes):return np.bincount(pi,(weights*nodes[wall[ci]]).sum(axis=1),minlength=len(area))/area
    def fields(folder,probe):
        state=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1,usecols=(4,9))
        with (probe/"restart.csv").open(encoding="ascii") as stream:
            keys=next(csv.reader(stream));rows=np.loadtxt(stream,delimiter=",",usecols=keys.index("Residual_Species_0"))
        flux=np.zeros(len(xyz));flux[wet]=rows[wet]*c["rho"]*3*c["cp"]*24*c["diameterM"]**2/lump[wet]
        return project(flux),project(c["rho"]*9*state[:,0]),24*state[first,1].mean(axis=1),24*project(state[:,1])
    qr,pr,tr,wr=fields(source/"plus20",source/"reaction20")
    q,p,t,w=fields(stage,stage/"reaction0")
    tab=np.genfromtxt(work/"surface/wall.csv",names=True,delimiter=",")
    np.testing.assert_array_equal(tab["faceID"],np.arange(len(area)))
    for name,field in (("heatFluxWm2",q),("pressureGaugePa",p),("nearWallCellTemperatureK",c["airK"]+t),("wallInputK",c["airK"]+w)):
        np.testing.assert_allclose(tab[name],field,rtol=0,atol=1e-8)
    assert np.all(qr!=0)
    power=lambda val:float(-area@val)
    error=dict(maxPressurePa=float(max(abs(p-pr))),maxRelativeHeat=float(max(abs(q/(f*qr)-1))),
               maxFluxWm2=float(max(abs(q-f*qr))),maxFirstK=float(max(abs(t-f*tr))),maxWallK=float(max(abs(w-f*wr))),
               relativePower=abs(power(q)/(f*power(qr))-1),referencePlus20=ref["plus20FrameSHA256"],
               correctedSteps=20,referenceSteps=20)
    assert error["maxPressurePa"]<.001 and error["maxRelativeHeat"]<.0005 and error["maxFirstK"]<.001 and error["relativePower"]<.0005,error
    for r in s["regions"]:
        use=group==body["regionNames"].index(r["name"])+1
        assert abs(r["outwardW"]+area[use]@q[use])<1e-12
    old=np.genfromtxt(base/"surface/wall.csv",names=True,delimiter=",")
    drift=dict(maxPressurePa=float(max(abs(tab["pressurePa"]-old["pressurePa"]))),maxHeatWm2=float(max(abs(q-f*old["heatFluxWm2"]))),
               maxFirstK=float(max(abs(t-f*(old["nearWallCellTemperatureK"]-c["airK"])))))
    out.mkdir(parents=True,exist_ok=True);archive=out/"resistance12.zip"
    with ZipFile(archive,"x") as bundle:
        for n in ("wall.csv","summary.json"):
            entry=ZipInfo(n);entry.compress_type=ZIP_DEFLATED;bundle.writestr(entry,(work/"surface"/n).read_bytes())
    data=dict(archiveSHA256=sha(archive),summarySHA256=sha(work/"surface/summary.json"),frameSHA256=s["sourceFrameSHA256"],
              wallSHA256=s["wallSHA256"],boundarySHA256=s["boundarySHA256"],factor=f,referenceFrame=request["wireFrame"],
              figure12=reading,matchedCheck=error,unmatchedToMain=drift,
              summaryBeadHReference="resistance-equivalent internal overtemperature, 12 K",
              nativeApparentHWm2K=s["beadAndJointsW"]/(c["sourceAreaM2"]*12),
              nativeSurfaceHWm2K=s["beadAndJointsW"]/(c["sourceAreaM2"]*(request["surfaceK"]-request["ambientK"])))
    (out/"resistance12.json").write_bytes((json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode("utf8"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("base","work","graph","out"):p.add_argument(k,type=Path)
    a=p.parse_args();package(*(x.resolve() for x in (a.base,a.work,a.graph,a.out)))
