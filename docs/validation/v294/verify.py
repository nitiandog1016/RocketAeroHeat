"""核对同一热敏珠的两步状态及固定端点热量账；每个阶段使用新的独立目录。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from advance import BUILD,GRID,HERE,LIMIT,PREVIOUS,ROOT,RUN,grid_info,observe,prefix,read,settings,sha
from check_state import audit
from read_wall import energy,sample
from thermistorStartupSurface import write
from thermistorStartupCache import bound_source
from thermistorStartupObserver import ensure_observer


def check(index,stage,area,out):
    label="first" if index==1 else "second"
    folder=area/(label+"01")
    if stage=="state":
        audit(index,GRID,folder,out,settings,LIMIT,grid_info,bound_source)
        report=read(out/"state.json")
        if index==1:prefix(folder/"history_00001.csv")
        report["verifiedPrefixRows"]=560 if index==1 else 0
        write(out/"verified.json",report)
    else:
        proof=area/(label+"_state")/"verified.json"
        if index==1:assert read(proof)["verifiedPrefixRows"]==560
        energy(index,GRID,folder,area/(label+"_observe"),area/(label+"_release"),proof,out,settings,LIMIT,BUILD,grid_info,bound_source)


def execute(stage,index,area,worker):
    label="first" if index==1 else "second";src=area/(label+"01")
    out=area/(label+"_"+stage);proof=area/(label+"_state")/"verified.json"
    if worker:
        assert stage in ("state","energy");check(index,stage,area,out);return
    grid_info()
    ensure_observer()
    out.mkdir(exist_ok=False)
    locks=[HERE/"PLAN.md",HERE/"advance.py",Path(__file__).resolve(),PREVIOUS/"run_step.py",PREVIOUS/"check_state.py",
           PREVIOUS/"read_wall.py",PREVIOUS/"evidence/first_history.csv.gz",PREVIOUS/"evidence/readout_build.json",
           GRID/"grid.json",GRID/"nodes.npz",ROOT/"scripts/thermistorStartupSurface.py",ROOT/"scripts/beadSurfaceResult.py",
           ROOT/"scripts/thermistorStartupCache.py",ROOT/"scripts/thermistorStartupObserver.py",bound_source("geometry.json")]
    if stage in ("observe","release"):
        sample(index,stage,area,out,GRID,src,proof,settings,LIMIT,BUILD,1100,locks,grid_info)
        return
    locks.extend(src/n for n in ("debugger_process.json","case.cfg",f"restart_{index:05}.csv",
                                f"seed_{index-1:05}.csv",f"history_{index:05}.csv"))
    if index==2:locks.append(area/"first_energy/accepted.json")
    if stage=="energy":
        if index==2:locks.extend([bound_source("body.json"),bound_source("nodes.npz")])
        locks.append(proof)
        locks.extend(area/(label+"_observe")/n for n in ("before_wall.csv","after_wall.csv","stages.csv","readout_updates.csv",
                                "debugger_process.json",f"restart_{index:05}.csv"))
        locks.extend(area/(label+"_release")/n for n in ("debugger_process.json",f"restart_{index:05}.csv"))
    rec=observe([sys.executable,str(Path(__file__).resolve()),stage,str(index),"--root",str(area),"--worker"],
                out,{str(p):sha(p) for p in locks},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
    if stage=="state" and not read(out/"verified.json")["passed"]:
        raise SystemExit("半步未通过物理验收，终场和诊断已保留。")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("state","observe","release","energy"));cli.add_argument("step",type=int,choices=(1,2))
    cli.add_argument("--root",type=Path,default=RUN);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);execute(a.stage,a.step,area,a.worker)
