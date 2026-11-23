"""双支管试件从登记冷态自动推进并对同一收敛端点执行热量核对。"""
import argparse
import os
from pathlib import Path

import numpy as np

from build_shared import HERE, PREV, KIT, ROOT, DRIVER, build, read, sha, write_json
from run import CASE, PREP, configuration, observe
from verify import budget, table


def advance(folder,start,exe,steps,release=False,capture=False,*,prepared=PREP,plan=HERE/"PLAN.md"):
    folder.mkdir()
    s=read(prepared/"prepared.json")
    names={"mesh.su2":prepared/"mesh.su2","wire_boundary.dat":prepared/("release.dat" if release else "profile.dat"),"initial_flow.csv":start}
    for name,path in names.items():os.link(path,folder/name)
    assert all(sha(folder/n)==s["files"][names[n].name] for n in ("mesh.su2","wire_boundary.dat"))
    cfg=folder/"case.cfg";cfg.write_text(configuration(s["condition"],steps,release,"BOUNDED_SCALAR"),encoding="ascii")
    if steps==200:assert sha(cfg)=="f1e0a615abb10c46560e1c446d5d8b61db9eac8d968aab748534fc1ca41a01cf"
    locks=[exe,cfg,*[folder/k for k in names],prepared/"prepared.json",CASE,plan,Path(__file__).resolve()]
    if capture:
        flag=folder/"capture.flag";flag.write_bytes(b"Short 65 endpoint\n");locks.append(flag)
    os.environ["PATH"]=str(KIT/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
    proc=observe([str(exe),"case.cfg","--threads","4"],folder,{str(p):sha(p) for p in locks},
                 300 if steps==200 else 150,450*1024**2)
    assert proc["status"]=="terminal" and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
    assert proc["debuggerExitCode"]==0 and "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace"),proc["debuggerExitCode"]
    h=table(folder,"history.csv")
    assert len(h)==steps and np.isfinite(h).all()
    np.testing.assert_array_equal(h[:,2],np.arange(steps))
    record=dict(start=sha(start),exe=sha(exe),end=sha(folder/"restart.csv"),history=sha(folder/"history.csv"),
                window=h[-20:,3:9].max(axis=0).tolist(),seconds=proc["elapsedSeconds"],peakBytes=proc["peakFamilyRSSBytes"])
    assert record["end"]==sha(folder/f"restart_{steps-1:06d}.csv")
    if not release:
        pairs=np.loadtxt(folder/"far_pairs.csv",delimiter=",",ndmin=2)
        np.testing.assert_array_equal(pairs[:,0],np.arange(steps))
        assert np.all(pairs[:,1]==1346) and np.all(pairs[:,2:4]>0) and np.all(pairs[:,4]<=1e-15)
        trace=np.genfromtxt(folder/"stages.csv",names=True,delimiter=",",dtype=None,encoding="ascii")
        assert len(trace)==steps*2 and np.all(trace["stage"][::2]=="low") and np.all(trace["stage"][1::2]=="actual")
        assert np.all(trace["below"]==0) and np.all(trace["above"]==0) and np.max(trace["wallDrift"])<2e-14
        record.update(rawMax20=float(max(trace["rawRMS"][1::2][-20:])),pairs=sha(folder/"far_pairs.csv"),trace=sha(folder/"stages.csv"))
        with np.load(prepared/"mapping.npz") as m:wet,wall=m["wet"],m["theta"]
        frames={}
        for index in sorted({steps-1,*range(100,steps,100)}):
            name=f"restart_{index:06d}.csv";state=table(folder,name);fields=s["condition"]
            np.testing.assert_array_equal(state[:,0],np.arange(s["geometry"]["nodes"]));assert np.isfinite(state).all()
            np.testing.assert_allclose(state[wet,9],wall,atol=2e-14,rtol=0)
            assert (fields["airPa"]+fields["rho"]*fields["speedMps"]**2*state[:,4]).min()>0
            assert state[:,9].min()>=-1e-12 and state[:,9].max()<=wall.max()+1e-12
            np.testing.assert_allclose([state[:,9].min(),state[:,9].max()],[trace[2*index+1][k] for k in ("min","max")],atol=1e-15,rtol=0)
            frames[name]=sha(folder/name)
        record["frames"]=frames
    write_json(folder/"record.json",record)
    return record


def solve(work,driver,*,prepared=PREP,plan=HERE/"PLAN.md"):
    work.mkdir(parents=True,exist_ok=False)
    record=read(driver/"build.json");exe=driver/"shared.exe"
    assert sha(exe)==record["exe"] and all(sha(driver/p)==v for p,v in record["files"].items())
    spec=read(prepared/"prepared.json")
    assert sha(prepared/"prepared.json")==record.get("preparedSHA256","ca9c622e056924cb0482e06e5f3561b402d7c40df079da93dd62cd61a98587b0")
    seed=prepared/"initial.csv"
    assert sha(seed)==spec["files"]["initial.csv"] and spec["caseSHA256"]==sha(CASE)
    trend=[]
    for i in range(1,13):
        stage=work/f"steady{i:02d}";r=advance(stage,seed,exe,200,prepared=prepared,plan=plan)
        worst=max(max(r["window"]),np.log10(r["rawMax20"]));trend.append(worst)
        print(stage.name,worst,flush=True);seed=stage/"restart.csv"
        if worst < -8:break
        if len(trend)>=3 and trend[-1]>=trend[-2]>=trend[-3]:raise RuntimeError("连续两段不改善")
    else:raise RuntimeError("达到登记阶段上限")
    write_json(work/"converged.json",dict(stage=stage.name,frame=sha(seed),trend=trend))
    release=DRIVER/"wire.exe";assert sha(release)=="9721e25124a866652e8d004b531854ec88d2f856e99ca7b80c0ab34d2bbebb55"
    for name,start,steps,let,snap in (("plus20",seed,20,False,False),("control",seed,1,False,False),
            ("early",seed,1,False,True),("releaseEarly",seed,1,True,False),
            ("late",work/"plus20/restart.csv",1,False,True),("releaseLate",work/"plus20/restart.csv",1,True,False)):
        advance(work/name,start,release if let else exe,steps,let,snap,prepared=prepared,plan=plan)
    check(work,prepared)


def check(work,prepared=PREP):
    conv=read(work/"converged.json");seed=conv["frame"]
    for name in ("early","late","control","plus20","releaseEarly","releaseLate"):
        d=work/name;r=read(d/"record.json")
        assert r["start"]==sha(d/"initial_flow.csv")==(read(work/"plus20/record.json")["end"] if name in ("late","releaseLate") else seed)
        assert r["end"]==sha(d/"restart.csv") and r["history"]==sha(d/"history.csv")
    controlError=float(np.max(abs(table(work/"control","restart.csv")-table(work/"early","restart.csv"))))
    assert controlError<2e-12 and max(read(work/"plus20/record.json")["window"])<-8
    s=read(prepared/"prepared.json")["condition"];scale=s["rho"]*s["cp"]*s["speedMps"]*s["diameterM"]**2*s["scalarUnitK"]
    with np.load(prepared/"mapping.npz") as m:wet=m["wet"]
    ends={}
    for name,rel in (("early","releaseEarly"),("late","releaseLate")):
        obs=work/name;ends[name]=budget(obs,work/rel,wet,scale)
        flow=np.loadtxt(obs/"flow_far.csv",delimiter=",");n=np.fromfile(obs/"nodes.bin","<f8").reshape(-1,16)
        n=n[n[:,0].argsort()];np.testing.assert_array_equal(n[flow[:,0].astype(int),11],flow[:,1])
    change=abs(ends["late"]["wallW"]/ends["early"]["wallW"]-1)
    passed=change<.002 and all(abs(e["physicalGapFraction"])<.005 for e in ends.values())
    result=dict(caseSHA256=sha(CASE),meshSHA256=sha(prepared/"mesh.su2"),**ends,scale=scale,thermalChange=change,controlError=controlError,
                passed=bool(passed),updates=len(conv["trend"])*200)
    write_json(work/"thermal.json",result)
    if not passed:raise ValueError("端点热量门槛未通过")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path);p.add_argument("driver",type=Path)
    a=p.parse_args();solve(a.output.resolve(),a.driver.resolve())
