"""从同一冷态回放已拒绝的质量通量方案，定位首次越界，读取对应边、边界和矩阵行。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];BASE=ROOT/"outputs/v326_validation"
sys.path.insert(0,str(ROOT/"scripts"))
from shortSphere65 import configuration,read,sha
from runSU2NativeDebugger import observe,atomic_write_json


def replay(directory,build,scout=None):
    original=BASE/"run01/steady01";prep=ROOT/"outputs/v325_validation/prep02/prepared"
    meta=read(prep/"prepared.json");s=meta["condition"]
    cfg=(original/"case.cfg").read_text("ascii")
    assert cfg==configuration(s,200,False,"BOUNDED_SCALAR") and (BASE/"run01/failure.json").is_file()
    r=read(original/"debugger_process.json")
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"])
    directory.mkdir(parents=True,exist_ok=False)
    locks={}
    for name,h in [("case.cfg",sha(original/"case.cfg")),("mesh.su2",meta["files"]["mesh.su2"]),
                   ("wire_boundary.dat",meta["files"]["profile.dat"]),("initial_flow.csv",meta["files"]["initial.csv"])]:
        src=original/name
        assert sha(src)==h==r["inputSHA256"][str(src)]
        os.link(src,directory/name);locks[str(directory/name)]=h
    step,point=-1,-1
    if scout:
        rec=read(scout/"replay.json")
        assert sha(scout/"extrema.csv")==rec["files"]["extrema.csv"]
        step,point=rec["firstIteration"],rec["worstGlobalID"]
    with np.load(prep/"mapping.npz") as m:ceiling=float(m["theta"].max())
    (directory/"probe.in").write_text(f"{step} {point} {ceiling:.17g}\n",encoding="ascii")
    exe=build/"probe.exe";br=read(build/"build.json")
    assert sha(exe)==br["exe"] and all(sha(build/n)==h for n,h in br["generated"].items())
    kit=ROOT/"tools/runtime/su2-thermal-build"
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    old,=[Path(p) for p in cache["cachedBuildInputs"] if Path(p).as_posix().endswith("/source/SU2_CFD/src/integration/CIntegration.cpp")]
    old_kit=old.parents[4]
    assert all(sha(kit/Path(p).relative_to(old_kit))==h for p,h in cache["cachedBuildInputs"].items())
    assert sha(kit/"su2.zip")==br["sourceArchive"]==cache["sourceArchiveSHA256"]
    for f in [exe,build/"build.json",directory/"probe.in",HERE/"PLAN.md",Path(__file__).resolve()]:locks[str(f)]=sha(f)
    os.environ["PATH"]=str(kit/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
    r=observe([str(exe),"case.cfg","--threads","4"],directory,locks,300,450*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"])
    assert "Exit Success (SU2_CFD)" in (directory/"debugger.log").read_text(errors="replace")
    hist=lambda d:np.loadtxt(d/"history.csv",skiprows=1,delimiter=",")
    ref=hist(original);run=hist(directory);max_log=float(abs(ref[:,3:9]-run[:,3:9]).max())
    assert len(run)==200 and max_log<.01,max_log
    trace=np.loadtxt(directory/"extrema.csv",delimiter=",")
    np.testing.assert_array_equal(trace[:,0],np.arange(200))
    flags=(trace[:,4]+trace[:,5])>0;first=trace[np.flatnonzero(flags)[0]]
    differences={}
    for i in [100,199]:
        n=f"restart_{i:06d}.csv";a=np.loadtxt(original/n,skiprows=1,delimiter=",");b=np.loadtxt(directory/n,skiprows=1,delimiter=",")
        differences[n]=float(abs(a-b).max())
        assert differences[n]<1e-8 and abs(trace[i,2]-a[:,9].min())<1e-13
        assert int(trace[i,4])==int(np.sum(a[:,9]<-1e-12)) and int(trace[i,5])==int(np.sum(a[:,9]>ceiling+1e-12))
    names=["extrema.csv","first_update.csv"]
    if scout:
        assert (step,point)==tuple(first[:2].astype(int))
        names+=["flow_edges.csv","scalar_edges.csv","flow_far.csv","scalar_far.csv","row.csv","matrix.csv"]
    atomic_write_json(directory/"replay.json",dict(firstIteration=int(first[0]),worstGlobalID=int(first[1]),maxLogDifference=max_log,
                      frameMaxDifference=differences,files={n:sha(directory/n) for n in names},build=br["exe"]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("name");p.add_argument("--build",default="build01");p.add_argument("--scout")
    a=p.parse_args();assert all(x.isidentifier() for x in [a.name,a.build,*([a.scout] if a.scout else [])])
    replay(BASE/a.name,BASE/a.build,BASE/a.scout if a.scout else None)
