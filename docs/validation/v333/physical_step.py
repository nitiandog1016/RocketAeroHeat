"""Short 双支管银球的独立物理步与释放壁行，按真实前一时间层接续。"""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np
from build_time import BASE, BASE_FRAME, BUILD, HERE, KIT, PREP, PREP_KEY, ROOT, read, sha, write_json
from shortSphere65 import CASE, configuration
from thermistorMeasuredCold import HEADER, RESIDUALS
from runSU2NativeDebugger import observe

DT={"full":.25,"half":.125}
UPDATES=360


def setup(work,driver,*,prepared=PREP,prepared_sha=PREP_KEY,steady=BASE/"run01/steady05",base_frame=BASE_FRAME,steady_build=BUILD):
    work.mkdir(parents=True,exist_ok=False)
    p=read(prepared/"prepared.json");assert sha(prepared/"prepared.json")==prepared_sha
    assert all(sha(prepared/name)==h for name,h in p["files"].items()) and sha(CASE)==p["caseSHA256"]
    frame=steady/"restart_000199.csv";end=read(steady/"record.json")
    rec=read(steady/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not(rec["supervisorError"] or rec["stopReason"] or rec["survivingProcessIDs"])
    assert sha(frame)==end["end"]==base_frame and max(end["window"])<-8 and end["rawMax20"]<1e-8
    heat=read(steady.parent/"thermal.json")
    assert heat["passed"] and heat["early"]["sourceFrame"]==base_frame and heat["meshSHA256"]==p["files"]["mesh.su2"]
    src=read(driver/"build.json");assert src["baseline"]==sha(steady_build/"build.json")
    assert all(sha(driver/name)==h for name,h in src["files"].items())
    state=np.loadtxt(frame,skiprows=1,delimiter=",")
    with np.load(prepared/"mapping.npz") as m:xyz,wet=m["xyz"],m["wet"]
    np.testing.assert_array_equal(xyz,state[:,1:4]);x,y,z=xyz.T
    dx=x+1.15;support=1-(dx/.55)**2-(y*y+z*z)/.45**2
    env=np.maximum(support,0.)**2
    delta=np.column_stack((6*z/.45**2*env,np.zeros(len(x)),-6*dx/.55**2*env))
    delta*=.01/np.linalg.norm(delta,axis=1).max()
    wall_or_far=np.r_[wet,np.flatnonzero(np.linalg.norm(xyz,axis=1)>13.9)]
    assert (support[wall_or_far]<0).all() and not delta[support<=0].any()
    prior=state.copy();state[:,5:8]+=delta
    path=work/"seed.csv";np.savetxt(path,state,delimiter=",",fmt=["%d"]+["%.17g"]*9,header=HEADER,comments="")
    np.testing.assert_array_equal(np.loadtxt(path,delimiter=",",skiprows=1),state)
    np.testing.assert_array_equal(prior[support<=0],state[support<=0])
    np.testing.assert_array_equal(prior[:,[0,1,2,3,4,8,9]],state[:,[0,1,2,3,4,8,9]])
    j=dict(baseFrame=base_frame,seed=sha(path),preparedSHA256=prepared_sha,driver=str(driver),buildSHA256=sha(driver/"build.json"),
           condition=p["condition"],nodes=len(xyz),wetNodes=len(wet),changed=int(np.any(prior!=state,axis=1).sum()),
           peakDeltaOverU=float(np.linalg.norm(state[:,5:8]-prior[:,5:8],axis=1).max()))
    write_json(work/"job.json",j);return j


def config(j,dt,step,release):
    limit=1 if release else UPDATES+1
    old=configuration(j["condition"],limit,release,"BOUNDED_SCALAR")
    drop={"ITER","WRT_RESTART_OVERWRITE","WRT_SURFACE_OVERWRITE","WRT_RESTART_COMPACT"}
    change=dict(CONV_STARTITER=str(limit+2),CONV_RESIDUAL_MINVAL="-9.5",SOLUTION_FILENAME="seed",OUTPUT_WRT_FREQ="(1, 1)")
    lines=[]
    for line in old.splitlines():
        k=line.partition("=")[0].strip()
        if k not in drop:lines.append(k+"= "+change.pop(k) if k in change else line)
    assert not change
    lines+=["TIME_DOMAIN= YES","TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER",f"TIME_STEP= {dt}","MAX_TIME= 10.0",
            f"RESTART_ITER= {step}",f"TIME_ITER= {step+1}",f"INNER_ITER= {limit}","HISTORY_WRT_FREQ_INNER= 1","WRT_RESTART_COMPACT= NO"]
    return "\n".join(lines)+"\n"


def advance(work,branch,step,old,release=False,*,prepared=PREP,plan=HERE/"PLAN.md"):
    assert (branch,step) in (("full",1),("half",1),("half",2))
    j=read(work/"job.json");driver=Path(j["driver"]);p=read(prepared/"prepared.json")
    assert sha(driver/"build.json")==j["buildSHA256"] and sha(prepared/"prepared.json")==j["preparedSHA256"]
    physical=work/branch/f"step{step}";dt=DT[branch]
    folder=physical/"release" if release else physical
    folder.mkdir(parents=True,exist_ok=False)
    if release:
        assert read(physical/"native.json")["frame"]==sha(old)
    else:
        expected=j["seed"] if step==1 else read(work/branch/f"step{step-1}"/"native.json")["frame"]
        assert sha(old)==expected
    exe=driver/("release.exe" if release else "physical.exe")
    assert sha(exe)==read(driver/"build.json")["files"][exe.name]
    for name,dest in (("mesh.su2","mesh.su2"),("release.dat" if release else "profile.dat","wire_boundary.dat")):
        assert sha(prepared/name)==p["files"][name];os.link(prepared/name,folder/dest)
    seed=folder/f"seed_{step-1:05d}.csv";os.link(old,seed)
    cfg=folder/"case.cfg";cfg.write_text(config(j,dt,step,release),encoding="ascii")
    locks=[exe,work/"job.json",plan,Path(__file__).resolve(),cfg,seed,folder/"wire_boundary.dat",folder/"mesh.su2",driver/"build.json"]
    binpath=str(KIT/"toolchain/w64devkit/bin")
    if binpath not in os.environ["PATH"].split(os.pathsep):os.environ["PATH"]=binpath+os.pathsep+os.environ["PATH"]
    r=observe([str(exe),"case.cfg","--threads","4"],folder,{str(f):sha(f) for f in locks},150 if release else 480,(350 if release else 550)*1024**2)
    return inspect_step(work,branch,step,old,release,prepared=prepared)


def inspect_step(work,branch,step,old,release,*,prepared=PREP):
    physical=work/branch/f"step{step}";folder=physical/"release" if release else physical
    j=read(work/"job.json");dt=DT[branch];seed=folder/f"seed_{step-1:05d}.csv"
    r=read(folder/"debugger_process.json");driver=Path(j["driver"]);p=read(prepared/"prepared.json")
    exe=driver/("release.exe" if release else "physical.exe")
    assert r["command"]==[str(exe),"case.cfg","--threads","4"]
    locks=r["inputSHA256"];assert sha(exe)==locks[str(exe)]==read(driver/"build.json")["files"][exe.name]
    for name,dest in (("mesh.su2","mesh.su2"),("release.dat" if release else "profile.dat","wire_boundary.dat")):
        assert sha(folder/dest)==p["files"][name]==locks[str(folder/dest)]
    assert sha(seed)==sha(old)==locks[str(seed)] and sha(work/"job.json")==locks[str(work/"job.json")]
    assert sha(folder/"case.cfg")==locks[str(folder/"case.cfg")] and (folder/"case.cfg").read_text("ascii")==config(j,dt,step,release)
    assert r["status"]=="terminal" and all(r["inputsStable"].values()) and not(r["survivingProcessIDs"] or r["supervisorError"] or r["stopReason"])
    assert r["debuggerExitCode"]==0 and "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace"),r["debuggerExitCode"]
    file=folder/f"history_{step:05d}.csv";frame=folder/f"restart_{step:05d}.csv"
    with file.open() as stream:
        keys=[v.strip() for v in next(csv.reader(stream,skipinitialspace=True))];h=np.loadtxt(stream,delimiter=",",ndmin=2)
    size=1 if release else UPDATES+1
    assert len(h)==size and np.all(h[:,keys.index("Time_Iter")]==step)
    if release:
        # 释放壁行时流动更新被显式跳过，五个零流动残差的 log10 为负无穷。
        frozen=[keys.index(k) for k in RESIDUALS[:-1]]
        assert np.isneginf(h[:,frozen]).all() and np.isfinite(np.delete(h,frozen,axis=1)).all()
    else:assert np.isfinite(h).all()
    np.testing.assert_array_equal(h[:,keys.index("Inner_Iter")],np.arange(size))
    info=dict(dt=dt,step=step,elapsedDU=dt*step,elapsedS=dt*step*j["condition"]["diameterM"]/j["condition"]["speedMps"],
              old=sha(seed),frame=sha(frame),history=sha(file),nativeElapsedSeconds=r["elapsedSeconds"],peakRSS=r["peakFamilyRSSBytes"])
    if not release:
        state=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=range(10))
        with np.load(prepared/"mapping.npz") as m:
            wet=m["wet"];theta=m["theta"];np.testing.assert_array_equal(state[:,1:4],m["xyz"])
        assert np.isfinite(state).all() and state[:,9].min()>=-1e-12 and state[:,9].max()<=theta.max()+1e-12
        assert np.max(abs(state[wet,5:8]))<1e-12 and min(j["condition"]["airPa"]+j["condition"]["rho"]*j["condition"]["speedMps"]**2*state[:,4])>0
        np.testing.assert_allclose(state[wet,9],theta,atol=2e-14,rtol=0)
        trace=np.genfromtxt(folder/"stages.csv",delimiter=",",names=True,dtype=None,encoding="ascii")
        np.testing.assert_array_equal(trace["iteration"],np.repeat(np.arange(size),2))
        assert np.all(trace["stage"][::2]=="low") and np.all(trace["stage"][1::2]=="actual")
        assert not trace["below"].any() and not trace["above"].any() and max(trace["wallDrift"])<=2e-14
        assert max(trace["rowError"])<2e-13 and np.all(abs(trace["balance"])<2e-12*np.maximum(1,trace["scale"]))
        pairs=np.loadtxt(folder/"far_pairs.csv",delimiter=",",ndmin=2)
        np.testing.assert_array_equal(pairs[:,0],np.arange(size))
        assert np.all(pairs[:,1]==2690) and np.all(pairs[:,2:4]>0) and np.all(pairs[:,4]<=1e-15)
        actual=trace[1::2];window=h[-21:-1,[keys.index(k) for k in RESIDUALS]].max(axis=0)
        info.update(updates=UPDATES,endpointAssemblies=1,window=window.tolist(),rawMax=max(actual["rawRMS"][-21:-1]),
                    heatChange=max(abs(actual["wallTotal"][-3:]/actual["wallTotal"][-1]-1)),maxRowError=max(trace["rowError"]),
                    vtk=sha(folder/f"surface_{step:05d}.vtk"))
        saved=folder/"observed.json"
        if saved.exists():assert read(saved)==info
        else:write_json(saved,info)
        assert max(info["window"])<-8 and info["rawMax"]<1e-8 and info["heatChange"]<.002,info
    saved=folder/"native.json"
    if saved.exists():assert read(saved)==info
    else:write_json(saved,info)
    return frame


def run_steps(work,*,prepared=PREP,plan=HERE/"PLAN.md",expected_edges=921290):
    for branch,steps in (("full",1),("half",2)):
        old=work/"seed.csv"
        for step in range(1,steps+1):
            folder=work/branch/f"step{step}"
            frame=inspect_step(work,branch,step,old,False,prepared=prepared) if folder.exists() else advance(work,branch,step,old,prepared=prepared,plan=plan)
            if (folder/"release").exists():inspect_step(work,branch,step,frame,True,prepared=prepared)
            else:advance(work,branch,step,frame,True,prepared=prepared,plan=plan)
            from energy import endpoint
            result=endpoint(work,branch,step,prepared=prepared,expected_edges=expected_edges)
            print(branch,step,result["wallW"],result["storageW"],flush=True);old=frame


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("work",type=Path);parser.add_argument("driver",type=Path);parser.add_argument("--resume",action="store_true")
    a=parser.parse_args();work=a.work.resolve();driver=a.driver.resolve();work.relative_to(ROOT/"outputs")
    if not a.resume:setup(work,driver)
    else:
        job=read(work/"job.json")
        assert job["buildSHA256"]==sha(driver/"build.json") and job["seed"]==sha(work/"seed.csv") and job["baseFrame"]==BASE_FRAME
    run_steps(work)
