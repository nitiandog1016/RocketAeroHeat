"""监督已登记的双支管银球阶段，从原冷态开始，保留拒绝现场或正常保存帧。"""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from shortSphere65 import CASE, configuration, read, sha, write_json
from runSU2NativeDebugger import observe

BASE=ROOT/"outputs/v328_validation"
PREP=ROOT/"outputs/v325_validation/prep02/prepared"
REF=ROOT/"outputs/v326_validation/run01/steady01"


def stage(out,build,start,spec):
    out.mkdir()
    br=read(build/"build.json");exe=build/"silver.exe"
    assert sha(exe)==br["exe"] and sha(HERE/"build.py")==br["recipe"]
    for name in ("silver.hpp","capacity.hpp","affine.cpp"):
        assert br["sources"][name]==sha(HERE/name)==sha(build/name)
    files={"initial_flow.csv":start,"mesh.su2":PREP/"mesh.su2","wire_boundary.dat":PREP/"profile.dat","case.cfg":REF/"case.cfg"}
    assert sha(files["case.cfg"])=="f1e0a615abb10c46560e1c446d5d8b61db9eac8d968aab748534fc1ca41a01cf"
    assert files["case.cfg"].read_text("ascii")==configuration(spec["condition"],200,False,"BOUNDED_SCALAR")
    assert sha(files["mesh.su2"])==spec["files"]["mesh.su2"] and sha(files["wire_boundary.dat"])==spec["files"]["profile.dat"]
    for name,path in files.items():
        os.link(path,out/name)
    locked=[out/n for n in files]+[exe,build/"build.json",PREP/"prepared.json",HERE/"PLAN.md",Path(__file__).resolve(),CASE]
    kit=ROOT/"tools/runtime/su2-thermal-build/toolchain/w64devkit/bin"
    os.environ["PATH"]=str(kit)+os.pathsep+os.environ["PATH"]
    proc=observe([str(exe),"case.cfg","--threads","4"],out,{str(p):sha(p) for p in locked},300,450*1024**2)
    assert proc["status"]=="terminal" and all(proc["inputsStable"].values())
    assert not(proc["survivingProcessIDs"] or proc["stopReason"] or proc["supervisorError"]),proc
    with (out/"stages.csv").open(encoding="ascii") as stream:rows=list(csv.DictReader(stream))
    code=proc["debuggerExitCode"]
    record=dict(exit=code,start=sha(start),seconds=proc["elapsedSeconds"],peakBytes=proc["peakFamilyRSSBytes"],trace=sha(out/"stages.csv"))
    if code==42:
        assert (out/"failure_nodes.csv").is_file() and (out/"failure_edges.csv").is_file()
        record.update(status="rejected",last=rows[-1],completedUpdates=sum(r["stage"]=="actual" for r in rows),
                      failureFiles={p.name:sha(p) for p in out.glob("failure_*.csv")})
    else:
        assert code==0 and "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")
        h=np.loadtxt(out/"history.csv",delimiter=",",skiprows=1,ndmin=2)
        assert len(h)==200 and len(rows)==400 and np.isfinite(h).all()
        frames={}
        with np.load(PREP/"mapping.npz") as m:wet,theta=m["wet"],m["theta"]
        for it in (100,199):
            file=out/f"restart_{it:06d}.csv";a=np.loadtxt(file,delimiter=",",skiprows=1)
            assert np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(len(a))) and len(a)==101122
            np.testing.assert_allclose(a[wet,9],theta,rtol=0,atol=2e-14)
            assert abs(a[:,9].min()-float(rows[2*it+1]["min"]))<1e-15 and a[:,9].min()>=-1e-12 and a[:,9].max()<=theta.max()+1e-12
            s=spec["condition"];assert (s["airPa"]+s["rho"]*s["speedMps"]**2*a[:,4]).min()>0
            frames[file.name]=sha(file)
        assert sha(out/"restart.csv")==frames["restart_000199.csv"]
        window=h[-20:,3:9].max(axis=0)
        raw=max(float(r["rawRMS"]) for r in rows[-40:] if r["stage"]=="actual")
        record.update(status="normal",window=window.tolist(),rawRMS=raw,frames=frames,history=sha(out/"history.csv"),
                      completedUpdates=200,passed=bool(max(window)<-8 and raw<1e-8))
    write_json(out/"stage.json",record)
    return record


def run(out,build):
    out.mkdir(parents=True,exist_ok=False)
    assert sha(PREP/"prepared.json")=="ca9c622e056924cb0482e06e5f3561b402d7c40df079da93dd62cd61a98587b0"
    spec=read(PREP/"prepared.json");start=PREP/"initial.csv"
    assert sha(start)==spec["files"]["initial.csv"] and spec["caseSHA256"]==sha(CASE)
    records=[];trend=[];reason="stage-budget"
    for i in range(1,13):
        folder=out/f"steady{i:02d}";r=stage(folder,build,start,spec);records.append(r)
        print(folder.name,r["status"],r.get("window"),flush=True)
        if r["status"]=="rejected":reason="method-gate";break
        if r["passed"]:reason="residuals-passed";break
        trend.append(max(max(r["window"]),np.log10(max(1e-300,r["rawRMS"]))))
        if len(trend)>=3 and trend[-1]>=trend[-2]>=trend[-3]:reason="two-nonimprovements";break
        start=folder/"restart.csv"
    write_json(out/"decision.json",dict(reason=reason,stages=len(records),qualifiedAreaM2=0,
               stageRecords=[sha(out/f"steady{i:02d}/stage.json") for i in range(1,len(records)+1)]))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("name");p.add_argument("--build",default="build01")
    a=p.parse_args();assert a.name.isidentifier() and a.build.isidentifier()
    run(BASE/a.name,BASE/a.build)
