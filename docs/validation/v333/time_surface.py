"""将两个时间步分支的真实同刻端点传回同一湿壁，逐原面比较压力、热流和近壁流向。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from build_time import ROOT, HERE, PREP, PREP_KEY, BASE_FRAME, sha, read
sys.path[:0]=[str(HERE.parent/k) for k in ("v331","v330")]
from direction import geometry, frame_vectors, mapped_fields, report_grid, FLAGS
from map_surface import project_surface, write_surface, DATA, CASE
from compare import differences


def store(path,doc):
    with path.open("x",encoding="utf8",newline="\n") as f:json.dump(doc,f,indent=2,allow_nan=False);f.write("\n")


def mapped_step(job,folder,step,branch,transfer,geom,out,*,prepared=PREP):
    e=read(folder/"energy.json");prep=read(prepared/"prepared.json");s=prep["condition"];dt=e["dt"]
    assert e["step"]==step and abs(e["relativeBalance"])<.005 and abs(e["assemblyGapW"])<1e-9
    assert e["heatChange"]<.002 and e["releasePeakRelativeError"]<1e-5 and max(e["window"])<-8 and e["rawMax"]<1e-8
    assert all(sha(folder/name)==key for name,key in e["files"].items())
    frame=folder/f"restart_{step:05d}.csv";vtk=folder/f"surface_{step:05d}.vtk"
    assert sha(frame)==e["frame"] and sha(vtk)==e["vtk"]
    state=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=range(10))
    with np.load(folder/"wall_budget.npz") as m:node_in=m["nodeInW"]
    stats=project_surface(out,prepared,transfer,geom["transferSHA256"],state,node_in,e["wallW"])
    w=np.genfromtxt(out/"wall.csv",names=True,delimiter=",");stations=np.genfromtxt(out/"stations.csv",names=True,delimiter=",")
    keys=dict(frame=e["frame"],vtkSHA256=e["vtk"])
    vectors,wet,scales,_=frame_vectors(prepared,frame,vtk,geom["m"],s,keys)
    values=mapped_fields(geom,vectors,scales);flux=-w["fluxInWm2"];area=w["areaM2"]
    peak_mask=flux>=.99*max(flux)
    peak=dict(faceID=stats["peakOutFace"],peakWm2=flux.max(),nearFaces=w["faceID"][peak_mask].astype(int).tolist(),nearAreaM2=area[peak_mask].sum())
    flow=report_grid(branch,geom,values,wet,vectors,scales,w,stations,{branch+"Peak":peak,branch+"W":e["wallW"]},sectors=128)
    np.savetxt(out/"direction.csv",np.column_stack((w["faceID"],w["regionID"],values)),delimiter=",",
               header="faceID,regionID,shearPa,shearMagnitudePa,speedMps,speedMagnitudeMps,status,flags",comments="",
               fmt=["%d","%d"]+["%.17g"]*4+["%d","%d"])
    summary=dict(schema="short65-supported-bdf1-v1",branch=branch,energy=e,condition=s,caseSHA256=sha(CASE),
        bodySHA256=prep["supportedSHA256"],transferSHA256=geom["transferSHA256"],meshSHA256=prep["files"]["mesh.su2"],
        frameSHA256=e["frame"],releaseSHA256=read(folder/"release/native.json")["frame"],baselineFrame=job["baseFrame"],seedFrame=job["seed"],
        dtDU=dt,dtSeconds=dt*s["diameterM"]/s["speedMps"],elapsedSeconds=e["elapsedS"],peak=peak,direction=flow,flagNames=FLAGS,**stats)
    summary["files"]["direction.csv"]=sha(out/"direction.csv")
    write_surface(out,summary,("direction.csv",))
    return summary,w,values,stations


def export(work,out,*,prepared=PREP,prepared_sha=PREP_KEY,base_frame=BASE_FRAME,
           steady_surface=DATA/"supportedAngularSurface65.zip"):
    job=read(work/"job.json");assert job["baseFrame"]==base_frame and sha(work/"seed.csv")==job["seed"]
    assert job["preparedSHA256"]==sha(prepared/"prepared.json")==prepared_sha
    manifest=read(DATA/"supportedAngular65.json");transfer=DATA/"supportedAngularTransfer65.npz"
    assert sha(transfer)==manifest["transferSHA256"]
    geom=geometry(manifest,transfer);geom["transferSHA256"]=sha(transfer)
    first=read(work/"half/step1/energy.json")
    assert abs(first["relativeBalance"])<.005 and first["old"]==job["seed"]
    out.mkdir(parents=True,exist_ok=False)
    full,w,v,st=mapped_step(job,work/"full/step1",1,"full",transfer,geom,out/"full",prepared=prepared)
    half,z,vh,sth=mapped_step(job,work/"half/step2",2,"half",transfer,geom,out/"half",prepared=prepared)
    assert full["energy"]["old"]==job["seed"] and half["energy"]["old"]==first["frame"]
    assert full["elapsedSeconds"]==half["elapsedSeconds"]
    for k in ("faceID","regionID","areaM2","wallInputK"):np.testing.assert_array_equal(w[k],z[k])
    s=job["condition"];area=w["areaM2"];dyn=s["rho"]*s["speedMps"]**2/2
    def rms(x):return np.sqrt(area@(x*x)/area.sum())
    dp=w["gaugePa"]-z["gaugePa"];dq=w["fluxInWm2"]-z["fluxInWm2"];scale=max(abs(z["fluxInWm2"]))
    power=abs(full["energy"]["wallW"]/half["energy"]["wallW"]-1)
    audit=dict(pressureRMSPa=rms(dp),pressureRMSOverDynamic=rms(dp)/dyn,pressureMaxPa=max(abs(dp)),heatRMSWm2=rms(dq),
        heatRMSOverPeak=rms(dq)/scale,heatMaxWm2=max(abs(dq)),heatMaxOverPeak=max(abs(dq))/scale,
        firstTemperatureMaxK=max(abs(w["firstK"]-z["firstK"])),wallPowerRelativeDifference=power)
    assert audit["pressureRMSOverDynamic"]<.01 and audit["heatRMSOverPeak"]<.01 and audit["heatMaxOverPeak"]<.05 and power<.01,audit
    with ZipFile(steady_surface) as archive:
        baseline=json.loads(archive.read("summary.json"))
        assert baseline["frameSHA256"]==base_frame and baseline["thermal"]["meshSHA256"]==full["meshSHA256"]
        sw=np.genfromtxt(archive.open("wall.csv"),names=True,delimiter=",")
    np.testing.assert_array_equal(w["faceID"],sw["faceID"])
    ds=z["fluxInWm2"]-sw["fluxInWm2"];ps=z["gaugePa"]-sw["gaugePa"]
    comparison=dict(schema="short65-supported-time-pair-v1",condition=s,elapsedS=full["elapsedSeconds"],elapsedDU=.25,seed=job["seed"],
        steadyFrame=base_frame,fullFrame=full["frameSHA256"],halfFrame=half["frameSHA256"],firstHalfFrame=first["frame"],
        fullArchive=sha(out/"full/surface.zip"),halfArchive=sha(out/"half/surface.zip"),meshSHA256=full["meshSHA256"],
        bodySHA256=full["bodySHA256"],caseSHA256=full["caseSHA256"],transferSHA256=sha(transfer),faces=len(area),areaM2=area.sum(),
        jointExperimentalAreaM2=0,stepComparison=audit,fullEnergy=full["energy"],halfEnergy=half["energy"],firstHalfEnergy=first,
        fullPeak=full["peak"],halfPeak=half["peak"],angles=st["angleDeg"].tolist(),experimentNu=st["reportedNu"].tolist(),
        fullNu=st["computedNu"].tolist(),halfNu=sth["computedNu"].tolist(),fullReverse=full["direction"]["total"]["reverse"],
        halfReverse=half["direction"]["total"]["reverse"],reverseIntersectionM2=area[(v[:,4]==-1)&(vh[:,4]==-1)].sum(),
        fullStations=full["direction"]["stations"],halfStations=half["direction"]["stations"],flagNames=FLAGS,
        fullTransitions=full["direction"]["transitions"],halfTransitions=half["direction"]["transitions"],
        maxHeatChangeFromSteadyWm2=max(abs(ds)),maxPressureChangeFromSteadyPa=max(abs(ps)),
        flowRangeMax={name:np.nanmax(np.abs(np.r_[v[:,j],vh[:,j]])) for j,name in enumerate(
            ("ShearPa","ShearMagnitudePa","SpeedMps","SpeedMagnitudeMps","Status"))},
        heatRangeKWm2=[min(w["fluxInWm2"].min(),z["fluxInWm2"].min())/1000,max(w["fluxInWm2"].max(),z["fluxInWm2"].max())/1000],
        gaugeRangePa=[min(w["gaugePa"].min(),z["gaugePa"].min()),max(w["gaugePa"].max(),z["gaugePa"].max())])
    store(out/"pair.json",comparison);return comparison


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("work",type=Path);p.add_argument("out",type=Path);a=p.parse_args();print(export(a.work,a.out)["stepComparison"])
