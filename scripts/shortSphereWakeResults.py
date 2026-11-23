"""按真实物理时刻比较两种步长，保留背风测站、尾迹速度及可在工作台读取的球面结果。"""
import argparse
import csv
import io
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from shortSphereWake import BASE, CASE, DT, ROOT, SAMPLES, read, sha, write_json
from shortSphere65 import conductivity


def samples(work,branch,time):
    step=round(time/DT[branch]);folder=work/branch/f"step{step:03d}"/"readout"
    seal=read(folder/"accepted.json")
    assert seal["timeDU"]==time and all(sha(folder/n)==h for n,h in seal["files"].items())
    with np.load(folder/"fields.npz") as z:fields={k:z[k] for k in z.files}
    return seal,fields


def wake_trace(work,branch,end,xyz,volume):
    # 同一 (x,r) 环内先转换到柱坐标，扣除周向均值；轴节点另纳入横向均速。
    x,y,z=xyz.T;rad=np.hypot(y,z);use=(x>=-3)&(x<=-.6)&(rad<=1)
    ring=use&(rad>1e-10);weights=volume[ring];cos=y[ring]/rad[ring];sin=z[ring]/rad[ring]
    _,group=np.unique(np.round(np.column_stack((x[ring],rad[ring])),11),axis=0,return_inverse=True)
    weight_sum=np.bincount(group,weights);dt=DT[branch];rows=[];j=read(work/"job.json")
    last=j["seedSHA256"]
    for step in range(round(end/dt)+1):
        if step:
            folder=work/branch/f"step{step:03d}";path=folder/f"restart_{step:05d}.csv";rec=read(folder/"native.json")
            assert rec["old"]==last and sha(path)==rec["frame"] and max(rec["window"].values())<-8
            last=rec["frame"]
        else:path=work/"seed.csv";assert sha(path)==last
        vel=np.loadtxt(path,delimiter=",",skiprows=1,usecols=(5,6,7))
        u,v,w=vel[ring].T;variance=0.
        for field in [u,v*cos+w*sin,-v*sin+w*cos]:
            mean=np.bincount(group,weights*field)/weight_sum
            variance+=weights@(field-mean[group])**2
        rows.append(dict(step=step,timeDU=step*dt,transverseMeanOverU=float(np.average(vel[use,2],weights=volume[use])),
                         nonAxisymmetricRmsOverU=float(np.sqrt(variance/weights.sum())),frame=last))
    return rows


def compare(work,through,out):
    out.mkdir(parents=True,exist_ok=False);assert through in SAMPLES
    j=read(work/"job.json");prep=Path(j["prepared"]);s=j["condition"];c=read(CASE)
    assert sha(prep/"prepared.json")==j["prepSHA256"]
    p=read(prep/"prepared.json");assert sha(prep/"mapping.npz")==p["files"]["mapping.npz"]
    with np.load(prep/"mapping.npz") as z:xyz,area=z["xyz"],z["area"]
    match=[];history={b:[] for b in DT}
    for time in SAMPLES:
        if time>through:break
        ends={b:samples(work,b,time) for b in DT}
        coarse,fine=(ends[b][1] for b in ("coarse","fine"))
        dP=fine["gauge"]-coarse["gauge"];dq=fine["q"]-coarse["q"];dT=fine["first"]-coarse["first"]
        def rms(v):return float(np.sqrt(np.average(v**2,weights=area)))
        pRms=rms(dP)/(.5*s["rho"]*s["speedMps"]**2);qref=max(abs(coarse["q"]))
        qRms=rms(dq)/qref;qMax=max(abs(dq))/qref;powerDiff=abs(ends["fine"][0]["wallW"]/ends["coarse"][0]["wallW"]-1)
        match.append(dict(timeDU=time,pressureRmsPa=rms(dP),pressureMaxPa=float(max(abs(dP))),pressureRmsOverDynamic=pRms,
                          heatRmsWm2=rms(dq),heatMaxWm2=float(max(abs(dq))),heatRmsOverPeak=qRms,heatMaxOverPeak=qMax,
                          powerRelative=powerDiff,firstRmsK=rms(dT),firstMaxK=float(max(abs(dT))),
                          withinTolerance=bool(pRms<.01 and qRms<.01 and qMax<.05 and powerDiff<.01)))
        for b,(seal,fields) in ends.items():
            entry={k:v for k,v in seal.items() if k not in ("files","stagesW","previousStagesW")}
            entry["observerDirectory"]=seal.get("observerDirectory","observer")
            history[b].append(entry)
    trace={b:wake_trace(work,b,through,xyz,ends[b][1]["volume"]) for b in DT}
    assert abs(ends["fine"][1]["volume"].sum()/p["geometry"]["volumeD3"]-1)<1e-11
    for b in DT:
        path=out/(b+"_wake.csv")
        with path.open("x",encoding="utf8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(trace[b][0]));writer.writeheader();writer.writerows(trace[b])
    with ZipFile(CASE.parent/"surface65.zip") as z:
        body=z.read("body.json");ref= json.loads(z.read("summary.json"));raw=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True)
        stations=np.genfromtxt(io.BytesIO(z.read("stations.csv")),delimiter=",",names=True)
    assert ref["frameSHA256"]==BASE and ref["meshSHA256"]==p["files"]["mesh.su2"]
    np.testing.assert_allclose(raw["areaM2"],area,rtol=1e-12,atol=0)
    summary=dict(timeDU=through,seed=j["seedSHA256"],base=BASE,case=sha(CASE),mesh=ref["meshSHA256"],
                 comparison=match,samples=history,trace=trace,statisticallyAveraged=False)
    write_json(out/"comparison.json",summary)
    if not all(m["withinTolerance"] for m in match):
        raise ValueError("同刻差异超过登记门槛，逐时刻差异已保留。")
    # 瞬时读数单独绑定参考面积、热功率与温度定义。压缩包不复制全域重启场。
    for b,(seal,fields) in ends.items():
        folder=out/b;folder.mkdir();w=raw.copy();w["gaugePa"]=fields["gauge"];w["absolutePa"]=fields["gauge"]+s["airPa"]
        w["fluxInWm2"]=fields["q"];w["wallInputK"]=fields["wall"];w["firstK"]=fields["first"]
        st=stations.copy();nu=np.asarray(seal["localNu"]);error=nu-st["reportedNu"]
        st["computedNu"]=nu;st["deltaNu"]=error;st["relativeError"]=error/st["reportedNu"];st["computedFluxInWm2"]=seal["stationFluxIn"]
        power=seal["wallW"];assert abs(w["fluxInWm2"]@area+power)<1e-11
        def table(data,path):
            np.savetxt(path,np.column_stack([data[k] for k in data.dtype.names]),delimiter=",",comments="",header=",".join(data.dtype.names),fmt="%.17g")
        table(w,folder/"wall.csv");table(st,folder/"stations.csv");(folder/"body.json").write_bytes(body)
        subset={k:ref[k] for k in ("condition","caseSHA256","meshSHA256","areaM2","analyticAreaM2","faces","observedW","gradientIntegratedW","reportedMeanWallK","zeroVerifiedAreaM2")}
        tw=np.average(fields["wall"],weights=area);equiv=power*subset["analyticAreaM2"]/area.sum()
        path=work/b/f"step{seal['physicalStep']:03d}"/"readout"
        release=read(path/"release_native/native.json");assert release["old"]==seal["frame"]
        subset.update(frameSHA256=seal["frame"],reactionSHA256=release["frame"],previousFrameSHA256=seal["old"],seedSHA256=j["seedSHA256"],
                      baseSHA256=BASE,dtSeconds=DT[b]*s["diameterM"]/s["speedMps"],timeDU=through,dtDU=DT[b],physicalStep=seal["physicalStep"],
                      elapsedPhysicalTimeS=seal["seconds"],wholeOutW=power,equivalentSphereW=equiv,equivalentPowerRelativeError=equiv/subset["observedW"]-1,
                      meanWallK=tw,meanNuFree=power*s["diameterM"]/(area.sum()*(tw-s["airK"])*conductivity(s["airK"])),
                      meanNuSurface=power*s["diameterM"]/(area.sum()*(tw-s["airK"])*conductivity(tw)),
                      stationBias=float(error.mean()),stationRMSE=float(np.sqrt(np.mean(error**2))),stationMaxAbs=float(max(abs(error))),
                      balance=dict(balance=seal["relativeBalance"],wallOutW=power,farOutW=seal["farW"],airStorageW=seal["storageW"],
                                   independentStorageW=seal["storageFromFramesW"],assemblyChange=seal["relativeWallChange"],window=seal["window"]),
                      wakeComparison=match[-1],history=history[b],wake=trace[b],baselineNu=stations["computedNu"].tolist(),
                      files={k:sha(folder/k) for k in ("wall.csv","stations.csv","body.json")})
        write_json(folder/"summary.json",subset);stem="wake_"+b
        with ZipFile(out/(stem+".zip"),"x",compression=ZIP_DEFLATED,compresslevel=9) as z:
            for name in ("body.json","wall.csv","stations.csv","summary.json"):z.write(folder/name,name)
        write_json(out/(stem+".json"),dict(archiveSHA256=sha(out/(stem+".zip")),summarySHA256=sha(folder/"summary.json"),caseSHA256=sha(CASE)))
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("out",type=Path);p.add_argument("--through",type=float,required=True)
    a=p.parse_args();compare(a.work.resolve(),a.through,a.out.resolve())
