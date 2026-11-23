"""核对银球基流同态收支，按完整球网积分并在论文实际极角取样。"""
import argparse
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from shortSphere65 import CASE, kelvin, conductivity, read, sha
from thermistorWireSurface import residual


def recover(work,out):
    out.mkdir(parents=True,exist_ok=False);prep=work/"prepared";p=read(prep/"prepared.json");s=p["condition"]
    c=read(CASE);end=read(work/"endpoints.json");stage=work/end["stage"]
    assert sha(CASE)==end["case"]==p["caseSHA256"] and all(sha(prep/n)==v for n,v in p["files"].items())
    phases=[stage,stage/"reaction",stage/"plus20",stage/"reaction20",stage/"control",stage/"observer"]
    seals=[]
    for dest in phases:
        a=read(dest/"native.json");r=read(dest/"debugger_process.json")
        assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
        assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
        assert a["frame"]==sha(dest/"restart.csv") and a["input"]==sha(dest/"initial_flow.csv")
        assert sha(dest/"history.csv")==a["history"] and sha(dest/"case.cfg")==r["inputSHA256"][str(dest/"case.cfg")]
        seals.append(a)
    base=seals[0]["frame"];assert base==end["frame"] and seals[0]["passed"] and seals[2]["passed"]
    assert all(seals[i]["input"]==base for i in [1,2,4,5]) and seals[3]["input"]==seals[2]["frame"]
    state=np.loadtxt(stage/"restart.csv",delimiter=",",skiprows=1)
    control,observed=(np.loadtxt(d/"restart.csv",delimiter=",",skiprows=1) for d in phases[-2:])
    same_step=float(abs(control-observed).max());assert same_step<2e-7
    with np.load(prep/"mapping.npz") as z:
        wet,first,v,f,centers,area,lump,ang,line,theta=(z[k] for k in ("wet","first","V","F","centers","area","lump","angles","line","theta"))
    np.testing.assert_allclose(state[wet,9],theta,rtol=0,atol=2e-14)
    scale=s["rho"]*s["speedMps"]*s["cp"]*s["scalarUnitK"]*s["diameterM"]**2
    raw=residual(phases[1]/"restart.csv");later=residual(phases[3]/"restart.csv")
    stages={}
    for name in ["start","beforeFar","afterFar","beforeWall","afterWall","end"]:
        rows=np.loadtxt(phases[-1]/(name+".csv"),delimiter=",",skiprows=1);rows=rows[rows[:,0].argsort()]
        np.testing.assert_array_equal(rows[:,0],np.arange(len(state)));np.testing.assert_array_equal(rows[:,1],state[:,9])
        stages[name]=float(rows[:,2].sum()*scale)
        if name=="beforeWall":node_error=float(abs(rows[wet,2]+raw[wet]).max()*scale)
    released=stages["beforeWall"]-stages["afterWall"];far=stages["afterFar"]-stages["beforeFar"]
    power=-raw[wet].sum()*scale
    change=float(abs(later[wet].sum()/raw[wet].sum()-1));balance=(released-far)/released
    gates=dict(wallOutW=power,farOutW=far,balance=balance,extra20Change=change,sameStepMax=same_step,
               wallRowErrorW=node_error,stagesW=stages,window=seals[0]["window"],extraWindow=seals[2]["window"])
    (out/"balance.json").write_text(json.dumps(gates,indent=2),encoding="utf8")
    assert min(power,far)>0 and node_error<1e-11 and abs(power-released)<1e-10 and abs(balance)<.005 and change<.002,gates
    # 三角面常值为线性场的面积平均；四边形中心以四个节点的平均值定义。
    def faces(x):return np.r_[x,x[centers].mean(axis=1)][f].mean(axis=1)
    q=scale*raw[wet]/lump
    gauge=s["rho"]*s["speedMps"]**2*state[wet,4]
    wall=s["airK"]+s["scalarUnitK"]*theta;near=s["airK"]+s["scalarUnitK"]*state[first,9]
    flux=faces(q);pressure=faces(gauge);tw=faces(wall);tc=faces(near)
    assert abs(flux@area+power)<1e-11 and np.all(s["airPa"]+pressure>0)
    coords=v[f].mean(axis=1);angle=np.degrees(np.arctan2(np.linalg.norm(coords[:,1:],axis=1),coords[:,0]))
    body=dict(V=v.tolist(),F=(f+1).tolist())
    (out/"body.json").write_text(json.dumps(body,separators=(",",":")),encoding="utf8")
    fields=["faceID","areaM2","angleDeg","absolutePa","gaugePa","fluxInWm2","wallInputK","firstK"]
    np.savetxt(out/"wall.csv",np.column_stack((np.arange(1,len(f)+1),area,angle,s["airPa"]+pressure,pressure,flux,tw,tc)),
               delimiter=",",header=",".join(fields),comments="",fmt=["%d"]+["%.17g"]*7)
    # 论文七站的极角从迎流点量起。在 +Z 子午线上读节点，不使用最近三角面的中心冒充驻点。
    sample=np.asarray(c["angleDeg"]);axis=ang[line];assert np.all(np.diff(axis)>0) and abs(axis[0])+abs(axis[-1]-180)<1e-10
    local_q=np.interp(sample,axis,q[line]);wall_sample=np.interp(sample,axis,wall[line])
    k_sample=conductivity(wall_sample);nu=-local_q*s["diameterM"]/((wall_sample-s["airK"])*k_sample)
    exp=np.array(c["localNuSurface"]);err=nu-exp
    refwall=kelvin(c["surfaceF"]);equiv=exp*conductivity(refwall)*(refwall-s["airK"])/s["diameterM"]
    np.savetxt(out/"stations.csv",np.column_stack((sample,exp,nu,err,err/exp,refwall,wall_sample,local_q,equiv)),
               header="angleDeg,reportedNu,computedNu,deltaNu,relativeError,reportedWallK,gridWallK,computedFluxInWm2,referenceEquivalentOutWm2",
               comments="",delimiter=",",fmt="%.17g")
    mean_wall=area@tw/area.sum();d=s["diameterM"];t0=s["airK"];BTU=1055.05585262
    observed_w=c["observedBtuS"]*BTU;integral_w=c["gradientIntegratedBtuS"]*BTU
    sphere_area=np.pi*d*d;equivalent_w=power*sphere_area/area.sum()
    summary=dict(condition=s,caseSHA256=sha(CASE),frameSHA256=base,reactionSHA256=seals[1]["frame"],meshSHA256=p["files"]["mesh.su2"],
                 observerFrameSHA256=seals[-1]["frame"],observerProgramSHA256=sha(Path(read(phases[-1]/"debugger_process.json")["command"][0])),
                 areaM2=float(area.sum()),analyticAreaM2=sphere_area,faces=len(f),
                 wholeOutW=power,equivalentSphereW=equivalent_w,observedW=observed_w,gradientIntegratedW=integral_w,
                 equivalentPowerRelativeError=equivalent_w/observed_w-1,meanWallK=mean_wall,reportedMeanWallK=float(kelvin(c["meanSurfaceF"])),
                 meanNuFree=power*d/(area.sum()*(mean_wall-t0)*conductivity(t0)),meanNuSurface=power*d/(area.sum()*(mean_wall-t0)*conductivity(mean_wall)),
                 stationBias=float(err.mean()),stationRMSE=float(np.sqrt((err**2).mean())),stationMaxAbs=float(abs(err).max()),
                 zeroVerifiedAreaM2=0.,balance=gates,files={k:sha(out/k) for k in ["wall.csv","stations.csv","body.json"]})
    (out/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf8")
    with ZipFile(out/"surface.zip","x",compression=ZIP_DEFLATED,compresslevel=9) as z:
        for name in ["body.json","wall.csv","stations.csv","summary.json"]:z.write(out/name,name)
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();s=recover(a.work.resolve(),a.out.resolve());print(s["wholeOutW"],s["stationRMSE"])
