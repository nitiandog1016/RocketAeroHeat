"""Short Test 65：在非轴对称扰动下推进同一光滑球，保存真实 BDF1 时间层。"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from runSU2NativeDebugger import observe
from shortSphere65 import CASE, DRIVER, ROOT, configuration, read, sha, write_json
from thermistorMeasuredCold import HEADER, RESIDUALS
from buildThermistorWireDriver import invoke

DT = {"coarse": .25, "fine": .125}
SAMPLES = (.25, 1., 2., 4., 6., 8.)
BASE = "794193ab0027408bea43686d9145e1e2e0b6d133f76ac37bab102641eca079d9"
TEMPLATE = ROOT/"cases/thermistor_b/startup6/CIntegration.cpp"
TEMPLATE_SHA = "8f769fc6aa061835806e5bb034e5a190984a820923510e2f0e4e76e5f669e274"


def setup(work,base,driver):
    work.mkdir(parents=True,exist_ok=False)
    prep=base/"prepared";ep=read(base/"endpoints.json");p=read(prep/"prepared.json");steady=base/ep["stage"]
    seal=read(steady/"native.json");native=read(steady/"debugger_process.json")
    assert ep["frame"]==BASE==sha(steady/"restart.csv")==seal["frame"] and seal["passed"]
    assert native["status"]=="terminal" and native["debuggerExitCode"]==0 and all(native["inputsStable"].values())
    assert not(native["stopReason"] or native["supervisorError"] or native["survivingProcessIDs"])
    key=read(CASE.parent/"surface65.json");archive=CASE.parent/"surface65.zip"
    assert sha(archive)==key["archiveSHA256"] and sha(CASE)==p["caseSHA256"]==key["caseSHA256"]
    with ZipFile(archive) as z:published=json.loads(z.read("summary.json"))
    assert published["frameSHA256"]==BASE and p["files"]["mesh.su2"]==published["meshSHA256"]
    assert all(sha(prep/n)==h for n,h in p["files"].items()) and sha(driver/"wire.exe")==DRIVER["wire"]
    state=np.loadtxt(steady/"restart.csv",delimiter=",",skiprows=1)
    with np.load(prep/"mapping.npz") as z:xyz,wet=z["xyz"],z["wet"]
    np.testing.assert_array_equal(state[:,1:4],xyz)
    # 矢势仅有 y 分量，三次紧支撑函数的旋度给出连续、无散度的初始扰动。
    x,y,z=xyz.T;dx=x+1.15;envelope=np.maximum(1-(dx/.55)**2-(y*y+z*z)/.45**2,0.)**2
    delta=np.column_stack((6*z/.45**2*envelope,np.zeros(len(x)),-6*dx/.55**2*envelope))
    delta*=.01/np.linalg.norm(delta,axis=1).max()
    assert np.any(delta) and np.all(delta[wet]==0) and np.all(delta[np.linalg.norm(xyz,axis=1)>13.9]==0)
    state[:,5:8]+=delta;seed=work/"seed.csv"
    np.savetxt(seed,state,delimiter=",",comments="",header=HEADER,fmt=["%d"]+["%.17g"]*9)
    write_json(work/"job.json",dict(base=str(base),baseFrame=BASE,prepared=str(prep),driver=str(driver),
        prepSHA256=sha(prep/"prepared.json"),seedSHA256=sha(seed),affectedNodes=int(np.count_nonzero(np.linalg.norm(delta,axis=1))),
        peakPerturbationOverU=float(np.linalg.norm(delta,axis=1).max()),condition=p["condition"]))


def observer_dir(work, workers):
    if workers not in (1, 4):
        raise ValueError("银球端点回放仅登记了单线程及原生四线程。")
    return work/("observer" if workers == 1 else "observer4")


def build_observer(work, workers=1):
    """只为 37,290 节点的本试件重定向已有固定端点读取器，与逐节点热界入口链接。"""
    j=read(work/"job.json");driver=Path(j["driver"]);out=observer_dir(work,workers);out.mkdir()
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build")).resolve()
    rec=read(driver/"build.json")
    assert rec["programs"]["wire"]==sha(driver/"wire.exe")==DRIVER["wire"] and sha(TEMPLATE)==TEMPLATE_SHA
    assert rec["driverSourceSHA256"]==sha(ROOT/"scripts/buildThermistorWireDriver.py")
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    old=Path(next(f for f in cache["cachedBuildInputs"] if f.endswith("CIntegration.cpp")))
    prefix=old.parents[3]
    locked={}
    # 构建缓存含源码与编译目录，两者共同以工具包根目录为起点。
    for name,h in cache["cachedBuildInputs"].items():locked[kit/Path(name).relative_to(prefix.parent)]=h
    assert all(sha(p)==h for p,h in locked.items())
    source=TEMPLATE.read_text("utf8")
    old_step="(config->GetTimeIter()==1 || config->GetTimeIter()==2)"
    assert source.count(old_step)==2 and source.count("==172482")==2
    source=source.replace("==172482","==37290").replace(old_step,"(config->GetTimeIter()>=1)")
    if workers == 4:
        # 取样前等所有线程结束装配，取样时只由主线程读写，随后共同继续或冻结。
        replacements = {
            "omp_get_num_threads()==1": ("omp_get_num_threads()==4", 2),
            "if (!snapshot) return;": ("if (!snapshot) return;\n#pragma omp barrier\n#pragma omp master\n    {", 1),
            'std::fclose(l);': ('std::fclose(l);\n    }\n#pragma omp barrier', 1),
            'auto* log=std::fopen': ('#pragma omp master\n    {\n    auto* log=std::fopen', 1),
            'std::fclose(log);': ('std::fclose(log);\n    }\n#pragma omp barrier', 1),
        }
        for before,(after,count) in replacements.items():
            assert source.count(before)==count
            source=source.replace(before,after)
    src=out/"CIntegration.cpp";src.write_text(source,encoding="utf8")
    for stage in ("CIntegration","wire_observer"):
        recipe=read(driver/(stage+"_process.json"));args=recipe["args"]
        assert recipe["exitCode"]==0 and not recipe["livePids"] and sha(Path(args[0]))==rec["compilerSHA256"]
        if stage=="CIntegration":
            args=[str(out/Path(x).name) if Path(x).parent==driver else x for x in args]
        else:
            args=[str(out/"wake_readout.exe") if Path(x)==driver/"wire_observer.exe" else
                  str(out/"CIntegration.obj") if Path(x)==driver/"CIntegration.obj" else x for x in args]
        invoke(stage,args,Path(recipe["cwd"]),out,180,3,bin_directory=Path(args[0]).parent)
    assert all(sha(p)==h for p,h in locked.items())
    write_json(out/"build.json",dict(exe=sha(out/"wake_readout.exe"),source=sha(src),template=TEMPLATE_SHA,
               mainObject=sha(driver/"SU2_CFD.obj"),originalDriver=DRIVER["wire"],compiler=rec["compilerSHA256"]))


def cfg(s,dt,step,limit,mode):
    text=configuration(s,limit,mode=="release")
    remove={"ITER","WRT_RESTART_OVERWRITE","WRT_SURFACE_OVERWRITE","WRT_RESTART_COMPACT"}
    change={"CONV_RESIDUAL_MINVAL":"-9.5","SOLUTION_FILENAME":"seed","OUTPUT_WRT_FREQ":"(1, 1)",
            "VOLUME_OUTPUT":"(COORDINATES, SOLUTION, RESIDUAL)","CONV_STARTITER":str(50 if mode=="advance" else limit+2)}
    result=[]
    for line in text.splitlines():
        key=line.split("=")[0].strip()
        if key not in remove:result.append(key+"= "+change.pop(key) if key in change else line)
    assert not change
    result+=["TIME_DOMAIN= YES","TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER",f"TIME_STEP= {dt:.17g}",
             "MAX_TIME= 10.0",f"RESTART_ITER= {step}",f"TIME_ITER= {step+1}",f"INNER_ITER= {limit}",
             "HISTORY_WRT_FREQ_INNER= 1","WRT_RESTART_COMPACT= NO"]
    return "\n".join(result)+"\n"


def native(work,branch,step,old,folder,limit,mode,observer_workers=1):
    j=read(work/"job.json");p=Path(j["prepared"]);assert sha(p/"prepared.json")==j["prepSHA256"]
    files=read(p/"prepared.json")["files"];dt=DT[branch];s=j["condition"]
    folder.mkdir(parents=True,exist_ok=False);seed=folder/f"seed_{step-1:05d}.csv"
    boundary="release.dat" if mode=="release" else "profile.dat"
    for name,target in [("mesh.su2","mesh.su2"),(boundary,"wire_boundary.dat")]:
        assert sha(p/name)==files[name];os.link(p/name,folder/target)
    os.link(old,seed);(folder/"case.cfg").write_text(cfg(s,dt,step,limit,mode),encoding="ascii")
    if mode in ("advance","release"):
        exe=Path(j["driver"])/"wire.exe";digest=DRIVER["wire"];threads=4 if mode=="advance" else 1
    else:
        reader=observer_dir(work,observer_workers)
        exe=reader/"wake_readout.exe";digest=read(reader/"build.json")["exe"];threads=observer_workers
    assert sha(exe)==digest
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build")).resolve()
    toolpath=str(kit/"toolchain/w64devkit/bin")
    if toolpath not in os.environ["PATH"].split(os.pathsep):os.environ["PATH"]=toolpath+os.pathsep+os.environ["PATH"]
    locks={str(q):sha(q) for q in [exe,seed,CASE,work/"job.json",Path(__file__).resolve(),*(folder/n for n in ("case.cfg","wire_boundary.dat","mesh.su2"))]}
    observe([str(exe),"case.cfg","--threads",str(threads)],folder,locks,120 if mode=="advance" else 180,(150 if mode=="advance" else 200)*1024**2)
    receipt=read(folder/"debugger_process.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not(receipt["stopReason"] or receipt["supervisorError"] or receipt["survivingProcessIDs"])
    with (folder/f"history_{step:05d}.csv").open() as f:
        keys=[k.strip() for k in next(csv.reader(f,skipinitialspace=True))];hist=np.loadtxt(f,delimiter=",",ndmin=2)
    assert np.isfinite(hist).all() and np.all(hist[:,keys.index("Time_Iter")]==step)
    np.testing.assert_array_equal(hist[:,keys.index("Inner_Iter")],np.arange(len(hist)))
    frame=folder/f"restart_{step:05d}.csv";state=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=range(10))
    with np.load(p/"mapping.npz") as z:xyz,wet,theta=z["xyz"],z["wet"],z["theta"]
    log=(folder/"debugger.log").read_text("ascii",errors="replace")
    assert f"Applied thermistor wall nodes: {len(wet)}" in log and "Exit Success (SU2_CFD)" in log
    np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-12)
    assert np.isfinite(state).all() and np.min(s["airPa"]+s["rho"]*s["speedMps"]**2*state[:,4])>0
    if mode!="release":
        assert state[:,9].min()>=-1e-12 and state[:,9].max()<=max(theta)+1e-12
        np.testing.assert_allclose(state[wet,9],theta,atol=2e-14,rtol=0)
    assert np.max(abs(state[wet,5:8]))<1e-10
    residuals=hist[:,[keys.index(k) for k in RESIDUALS]]
    window=dict(zip(RESIDUALS,np.max(residuals[-20:],axis=0).tolist()))
    record=dict(step=step,dt=dt,inner=len(hist),old=sha(old),frame=sha(frame),history=sha(folder/f"history_{step:05d}.csv"),window=window,
                seconds=receipt["elapsedSeconds"],peakRSS=receipt["peakFamilyRSSBytes"],mode=mode,positive=True)
    write_json(folder/"native.json",record)
    if mode=="advance":assert len(hist)>=20 and max(window.values()) < -8,window
    else:assert len(hist)==limit
    return state,residuals,record


def advance(work,branch,through):
    """每个物理步均独立保存；只有已核对的同一父帧才能接续。"""
    dt=DT[branch];end=int(round(through/dt));assert 0<end*dt==through<=8
    j=read(work/"job.json");old=work/"seed.csv";assert sha(old)==j["seedSHA256"]
    for step in range(1,end+1):
        folder=work/branch/f"step{step:03d}";frame=folder/f"restart_{step:05d}.csv"
        if folder.exists():
            rec=read(folder/"native.json");r=read(folder/"debugger_process.json")
            assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values()) and not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
            assert rec["old"]==sha(old) and rec["frame"]==sha(frame) and rec["inner"]>=20 and max(rec["window"].values())<-8
            assert rec["step"]==step and rec["dt"]==dt and rec["mode"]=="advance"
            assert rec["history"]==sha(folder/f"history_{step:05d}.csv")
            assert (folder/"case.cfg").read_text("ascii")==cfg(j["condition"],dt,step,360,"advance")
            prep=Path(j["prepared"]);assert sha(prep/"prepared.json")==j["prepSHA256"]
            files=read(prep/"prepared.json")["files"];locks=r["inputSHA256"]
            for name,digest in [(f"seed_{step-1:05d}.csv",rec["old"]),("wire_boundary.dat",files["profile.dat"]),("mesh.su2",files["mesh.su2"])]:
                assert sha(folder/name)==locks[str(folder/name)]==digest
            assert sha(folder/"case.cfg")==locks[str(folder/"case.cfg")]
            assert locks[str(work/"job.json")]==sha(work/"job.json") and locks[str(CASE)]==sha(CASE)
            assert r["command"][1:]==["case.cfg","--threads","4"] and locks[r["command"][0]]==DRIVER["wire"]
        else:
            _,_,rec=native(work,branch,step,old,folder,360,"advance")
            print(branch,step,rec["inner"],max(rec["window"].values()),flush=True)
        old=frame


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=("setup","build","advance"))
    parser.add_argument("work",type=Path);parser.add_argument("--base",type=Path,default=ROOT/"outputs/v321_validation/run01")
    parser.add_argument("--driver",type=Path,default=ROOT/"outputs/v317_validation/build01")
    parser.add_argument("--branch",choices=DT);parser.add_argument("--through",type=float);a=parser.parse_args()
    w=a.work.resolve();w.relative_to((ROOT/"outputs").resolve())
    if a.action=="setup":setup(w,a.base.resolve(),a.driver.resolve())
    elif a.action=="build":build_observer(w)
    else:advance(w,a.branch,a.through)
