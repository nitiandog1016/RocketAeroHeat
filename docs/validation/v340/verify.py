"""从独立释放的壁行、原生内部边和保存体场复算量热边界，保留可离线审计的小包。"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile, ZIP_DEFLATED

for var in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):os.environ[var]="1"
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CASE=ROOT/"cases/short_sphere";PROOF=HERE/"evidence"
sys.path[:0]=[str(HERE.parent/"v331"),str(ROOT/"scripts")]
from verify_direction import normals_and_weights
from shortSphere65 import conductivity, read, sha
from thermistorWireSurface import residual


def store(path,obj):
    with path.open("x",encoding="utf8",newline="\n") as f:json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False);f.write("\n")


def publish(job,stem="calorimetry65"):
    ack=read(job/"completed.json");assert sha(job/"power65.zip")==ack["archiveSHA256"]
    # 发布回执保留来源摘要；本机大场路径及原始回执保存在证据包。
    source=ack["source"]
    ack["source"]={k:v for k,v in source.items() if k not in ("prepared","frame")}
    archive=CASE/(stem+".zip");manifest=CASE/(stem+".json")
    if manifest.exists():assert read(manifest)==ack
    if archive.exists():assert sha(archive)==ack["archiveSHA256"]
    else:shutil.copyfile(job/"power65.zip",archive)
    if not manifest.exists():store(manifest,ack)


def terminal(p,exe):
    d=read(p/"debugger_process.json")
    assert d["status"]=="terminal" and d["debuggerExitCode"]==0 and all(d["inputsStable"].values())
    assert not(d["stopReason"] or d["survivingProcessIDs"] or d["supervisorError"])
    assert "Exit Success (SU2_CFD)" in (p/"debugger.log").read_text(errors="replace")
    assert d["inputSHA256"][d["command"][0]]==exe and sha(p/"case.cfg")==d["inputSHA256"][str(p/"case.cfg")]
    return d


def pack(job,*,film=False,proof_dir=PROOF):
    transfer=CASE/"supportedAngularTransfer65.npz"
    with np.load(transfer) as t:xyz=t["nativeWallXYZD"];firstxyz=t["firstXYZD"]
    ack=read(job/"completed.json");included={};large={}
    source=read(job/"source.json");base=Path(source["prepared"]);parent=Path(source["frame"])
    assert sha(parent)==source["frameSHA256"] and sha(base/"prepared.json")==source["preparedSHA256"]
    original=np.loadtxt(parent,delimiter=",",skiprows=1,usecols=range(10))
    with np.load(base/"mapping.npz") as mapping:wet_ids=mapping["wet"];old_wall=mapping["theta"]
    previous_excess=None
    with ZipFile(proof_dir/"native.zip","x",compression=ZIP_DEFLATED,compresslevel=6) as z:
        def put(name,data):
            included[name]=hashlib.sha256(data).hexdigest();z.writestr(name,data)
        def add(f):put(f.relative_to(job).as_posix(),f.read_bytes())
        add(job/"source.json");add(job/"completed.json")
        if film:
            add(job/"acceptedTrials.json")
            for f in (ROOT/"scripts/shortSupportFilmProperties.py",ROOT/"scripts/shortSupportFilm65.py",
                      ROOT/"scripts/shortSupportPower65.py",HERE.parent/"v329/run_shared.py"):
                put("source/"+f.name,f.read_bytes())
        for lab in (ack["order"] if film else ("high","low","balance")):
            trial=job/lab;work=trial/"run";p=ack["trials"][lab] if film else ack[lab];build=read(trial/"driver/build.json");exe=build["exe"]
            assert sha(trial/"driver/uniform.exe")==exe and abs(build["expectedTheta"]*90-p["excessK"])<1e-12
            for f in ("trial.json","prepared/prepared.json","driver/build.json","driver/CScalarSolver.inl","driver/species_process.json",
                      "driver/link_process.json","run/converged.json","run/thermal.json"):add(trial/f)
            thermal=read(work/"thermal.json");prior=sha(trial/"prepared/initial.csv")
            spec=read(trial/"prepared/prepared.json")
            assert prior==p["initialSHA256"] and spec["parentFrame"]==sha(parent)
            seed=np.loadtxt(trial/"prepared/initial.csv",delimiter=",",skiprows=1)
            np.testing.assert_array_equal(seed[:,:8],original[:,:8])
            new_wall=p["excessK"]/90
            if film:
                assert spec["initialSourceExcessK"]==previous_excess and p["parentFrame"]==sha(parent)
                if previous_excess is None:
                    theta=np.zeros(len(seed));enthalpy=np.full(len(seed),1-298.15);enthalpy[wet_ids]+=new_wall
                else:
                    ratio=p["excessK"]/previous_excess
                    theta=original[:,9]*ratio;enthalpy=(1-298.15)+(original[:,8]-(1-298.15))*ratio
                np.testing.assert_array_equal(seed[:,8],enthalpy)
            else:theta=original[:,9]*(new_wall/max(old_wall))
            theta[wet_ids]=new_wall
            np.testing.assert_array_equal(seed[:,9],theta)
            for boundary,scalar in (("profile.dat",new_wall),("release.dat",0.)):
                rows=np.loadtxt(trial/"prepared"/boundary);np.testing.assert_array_equal(rows[:,0],wet_ids)
                np.testing.assert_array_equal(rows[:,1],np.full(len(rows),1+new_wall))
                np.testing.assert_array_equal(rows[:,2],np.full(len(rows),scalar))
            dirs=sorted(work.glob("steady*"))
            for stage in dirs+[work/name for name in ("plus20","early","control","late","releaseEarly","releaseLate")]:
                record=read(stage/"record.json");release=stage.name.startswith("release")
                receipt=terminal(stage,ack["source"]["releaseSHA256"] if release else exe)
                if film:assert receipt["inputSHA256"][str(HERE.parent/"v329/run_shared.py")]==included["source/run_shared.py"]
                for name in ("record.json","case.cfg","debugger_process.json","history.csv","debugger.log","stages.csv","far_pairs.csv"):
                    path=stage/name
                    if path.is_file():add(path)
                assert sha(stage/"initial_flow.csv")==record["start"] and sha(stage/"restart.csv")==record["end"]
                for name,digest in record.get("frames",{}).items():
                    file=stage/name;assert sha(file)==digest;large[file.relative_to(job).as_posix()]=digest
                if stage in dirs:
                    assert record["start"]==prior;prior=record["end"]
                    trace=np.genfromtxt(stage/"stages.csv",names=True,delimiter=",",dtype=None,encoding="ascii")
                    assert len(trace)==400 and not trace["below"].any() and not trace["above"].any()
                    assert max(trace["rowError"])<2e-13 and np.all(abs(trace["balance"])<2e-12*np.maximum(1,trace["scale"]))
            assert prior==p["earlyFrame"]==thermal["early"]["sourceFrame"] and thermal["passed"]
            end=read(dirs[-1]/"record.json");assert max(end["window"])<-8 and end["rawMax20"]<1e-8
            for label,release in (("early","releaseEarly"),("late","releaseLate")):
                stage=work/label;e=thermal[label];frame=stage/"initial_flow.csv";edges=stage/"edges.bin";far=stage/"flow_far.csv"
                for f,h in ((frame,e["sourceFrame"]),(edges,e["binaryFiles"]["edges.bin"]),(far,e["binaryFiles"]["flow_far.csv"])):
                    assert sha(f)==h;large[f.relative_to(job).as_posix()]=h
                state=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=range(10));tree=cKDTree(state[:,1:4])
                dw,wet=tree.query(xyz);df,first=tree.query(firstxyz)
                assert max(dw.max(),df.max())<1e-12 and len(np.unique(wet))==11906
                edge=np.fromfile(edges,"<f8").reshape(-1,6);i,j=edge[:,:2].astype(int).T
                qty=edge[:,4]*(edge[:,5]-1)
                correction=np.bincount(i,qty,minlength=len(state))-np.bincount(j,qty,minlength=len(state))
                freed=residual(work/release/"restart.csv")[wet]
                sides=np.loadtxt(far,delimiter=",");ids=sides[:,0].astype(int)
                buf=io.BytesIO()
                np.savez_compressed(buf,wet=wet,release=freed,correction=correction[wet],p=state[wet,4],theta=state[wet,9],
                                    first=state[first,9],farMass=sides[:,1],farTheta=state[ids,9])
                put(lab+"/"+label+".npz",buf.getvalue())
            if film:
                previous_excess=p["excessK"];parent=work/"early/initial_flow.csv"
                original=np.loadtxt(parent,delimiter=",",skiprows=1,usecols=range(10))
        failed=job.parent/"power01/high/run/steady01"
        if not film and failed.is_dir():
            for k in ("debugger_process.json","debugger.log","history.csv","case.cfg"):
                put("first_attempt/"+k,(failed/k).read_bytes())
            put("first_attempt/attempt_driver.py",(job.parent/"power01/attempt_driver.py").read_bytes())
            assert sha(job.parent/"power01/high/prepared/initial.csv")==ack["high"]["initialSHA256"]
        for f in (proof_dir.parent/"PLAN.md",HERE/"NATIVE.md"):
            put("contracts/"+f.name,f.read_bytes())
        z.writestr("manifest.json",json.dumps(dict(files=included,largeLocal=large),indent=2))


def verify(*,film=False,proof_dir=PROOF,stem="calorimetry65"):
    p=read(CASE/(stem+".json"));case=read(CASE/"test65.json");archive=CASE/(stem+".zip")
    assert p["schema"]==("short65-sphere-power-film-v1" if film else "short65-sphere-power-v1")
    final=p["balanceName"] if film else "balance"
    assert p["archiveSHA256"]==sha(archive) and p["source"]["caseSHA256"]==sha(CASE/"test65.json")
    target=case["observedBtuS"]*1055.05585262/case["supportAreaFactor"];assert p["targetSphereW"]==target
    with ZipFile(archive) as z:
        s=json.loads(z.read("summary.json"));w=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")
        stations=np.genfromtxt(z.open("stations.csv"),names=True,delimiter=",")
        for f,h in s["files"].items():assert hashlib.sha256(z.read(f)).hexdigest()==h
    transfer=CASE/"supportedAngularTransfer65.npz";assert sha(transfer)==s["transferSHA256"]==p["source"]["transferSHA256"]
    with np.load(transfer) as grid:g=dict(grid)
    area=g["targetArea"];region=g["region"];src=g["nativeArea"];corners=g["nativeCorners"];xyz=g["nativeWallXYZD"]
    _,frac=normals_and_weights(xyz,corners)
    B=csr_matrix(((src[:,None]*frac).ravel(),(np.repeat(np.arange(len(src)),4),corners.ravel())),shape=(len(src),len(xyz)))
    lump=B.sum(axis=0).A1
    Q=csr_matrix((g["sourceFraction"],(g["target"],g["source"])),shape=(len(area),len(src)))
    I=csr_matrix((g["targetFraction"],(g["target"],g["source"])),shape=Q.shape)
    def mean(v):return I@((B@v)/src)
    np.testing.assert_array_equal(w["faceID"],np.arange(len(area))+1);np.testing.assert_array_equal(w["regionID"],region)
    np.testing.assert_array_equal(w["areaM2"],area);assert len(area)==32480 and s["jointExperimentalAreaM2"]==0
    cond=s["condition"];air=cond["airK"];proof={};per_degree={}
    with ZipFile(proof_dir/"native.zip") as z:
        manifest=json.loads(z.read("manifest.json"))
        for f,h in manifest["files"].items():assert hashlib.sha256(z.read(f)).hexdigest()==h
        full=json.loads(z.read("completed.json"))
        assert all(full[k]==p[k] for k in p if k!="source")
        assert all(full["source"][k]==v for k,v in p["source"].items())
        for lab in (p["order"] if film else ("high","low","balance")):
            trial=p["trials"][lab] if film else p[lab];thermal=json.loads(z.read(lab+"/run/thermal.json"));rows={}
            conv=json.loads(z.read(lab+"/run/converged.json"));previous=trial["initialSHA256"]
            for index in range(1,len(conv["trend"])+1):
                stage=f"{lab}/run/steady{index:02d}/";record=json.loads(z.read(stage+"record.json"))
                proc=json.loads(z.read(stage+"debugger_process.json"))
                assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
                assert not(proc["stopReason"] or proc["survivingProcessIDs"] or proc["supervisorError"])
                history=np.loadtxt(z.open(stage+"history.csv"),delimiter=",",skiprows=1)
                trace=np.genfromtxt(z.open(stage+"stages.csv"),names=True,delimiter=",",dtype=None,encoding="ascii")
                assert record["start"]==previous and len(history)==200 and len(trace)==400
                previous=record["end"];window=history[-20:,3:9].max(axis=0);raw=trace["rawRMS"][1::2][-20:].max()
                np.testing.assert_array_equal(window,record["window"]);assert raw==record["rawMax20"]
            assert max(window)<-8 and raw<1e-8 and previous==conv["frame"]==trial["earlyFrame"]
            scale=cond["rho"]*cond["cp"]*cond["speedMps"]*cond["diameterM"]**2*cond["scalarUnitK"]
            assert scale==thermal["scale"] and thermal["thermalChange"]<.002
            for label in ("early","late"):
                with np.load(io.BytesIO(z.read(lab+"/"+label+".npz"))) as v:
                    np.testing.assert_allclose(v["theta"],trial["excessK"]/cond["scalarUnitK"],atol=2e-14,rtol=0)
                    load=(v["release"]-v["correction"])*scale;q=Q@(B@(load/lump))/area
                    zones=[-q[region==k]@area[region==k] for k in range(1,6)]
                    total=sum(zones);farW=np.maximum(v["farMass"],0)@v["farTheta"]*scale
                    np.testing.assert_allclose(zones,trial[label+"RegionW"],atol=1e-10,rtol=0)
                    assert abs(total-thermal[label]["wallW"])<1e-10 and abs(farW-thermal[label]["scalarFarW"])<1e-10
                    assert abs(total-farW)/total<.005 and np.argmin(q)+1==trial[label+"Peak"]
                    rows[label]=dict(regionW=zones,totalW=total,farW=farW,peakFace=int(np.argmin(q)+1),peakWm2=-min(q))
                    if lab==final:
                        assert abs(zones[0]/target-1)<.0002
                    if label=="early":
                        per_degree[lab]=(q/trial["excessK"],mean(v["first"])*cond["scalarUnitK"]/trial["excessK"])
                    if lab!=final or label!="early":continue
                    gauge=mean(cond["rho"]*cond["speedMps"]**2*v["p"]);wt=air+cond["scalarUnitK"]*v["theta"]
                    expected=dict(absolutePa=cond["airPa"]+gauge,gaugePa=gauge,fluxInWm2=q,wallInputK=mean(wt),
                                  firstK=mean(air+cond["scalarUnitK"]*v["first"]))
                    errors={k:max(abs(w[k]-x)) for k,x in expected.items()};assert max(errors.values())<1e-8
                    rows["faceErrors"]=errors
                    look={j:i for i,j in enumerate(v["wet"])}
                    for st in stations:
                        ix=[look[int(st[k])] for k in ("leftNode","rightNode")];mix=st["rightWeight"];weights=np.array([1-mix,mix])
                        temp=weights@wt[ix];heat=weights@(load/lump)[ix]
                        nu=-heat*cond["diameterM"]/((temp-air)*conductivity(temp))
                        assert abs(nu-st["computedNu"])<1e-9
            proof[lab]=rows
    q0=proof["low"]["early"]["regionW"][0];q1=proof["high"]["early"]["regionW"][0]
    assert q0<target<q1
    candidate=p["temperatureK"]
    if not film:
        assert abs(q1/(2*q0)-1)<.002
        candidate=air+20+20*(target-q0)/(q1-q0)
        assert abs(candidate-p["temperatureK"])<1e-8
    err=[max(abs(a-b)) for a,b in zip(per_degree["high"],per_degree["low"])]
    if not film:assert err[0]<.002*max(abs(per_degree["high"][0])) and err[1]<.002
    return dict(trials=proof,targetSphereW=target,temperatureK=candidate,experimentalMeanK=(case["meanSurfaceF"]-32)*5/9+273.15,
                normalizedMaxErrors=err,faces=len(w),areaM2=sum(area),peakFace=int(np.argmin(w["fluxInWm2"])+1),
                sphereFaces=int(sum(region==1)),nu=stations["computedNu"].tolist(),experimentalNu=case["localNuSurface"],
                stationRMSE=float(np.sqrt(np.mean((stations["computedNu"]-stations["reportedNu"])**2))))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--publish",type=Path);parser.add_argument("--pack",type=Path)
    parser.add_argument("--record",action="store_true");a=parser.parse_args();PROOF.mkdir(exist_ok=True)
    if a.publish:publish(a.publish.resolve())
    if a.pack:pack(a.pack.resolve())
    if not a.publish or a.pack:
        result=verify()
        if a.record:store(PROOF/"comparison.json",result)
        print("Verified 32,480 faces, three equilibrium trials, six energy balances and seven heat-transfer stations.")
