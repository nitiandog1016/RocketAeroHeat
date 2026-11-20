"""复核热敏珠保存帧及旧新时间层，并对实际已接受的半步核实同阶段热量账。"""
import argparse
import os
from pathlib import Path
import sys

for k in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[k]="1"

from run_step import HERE, LIMIT, PREVIOUS, ROOT, RUN, input_file, observe, prefix, read, settings, sha
from check_half import RMS, columns, energy, state
from observe_half import build, sample
from thermistorStartupSurface import write


def check(stage,index,folder,out,area):
    label="first" if index==1 else "second"
    history=folder/f"history_{index:05}.csv"
    if stage=="state":
        state(index,folder,out,settings,LIMIT)
        record=read(out/"state.json")
        if index==1:
            prefix(history)
            assert record["oldFrameSHA256"]==sha(input_file("seed_00000.csv"))
        else:
            parent=read(area/"first_energy/accepted.json")
            assert parent["accepted"] and record["oldFrameSHA256"]==parent["energy"]["acceptedFrameSHA256"]
        good=(columns(history,RMS)<-8).all(axis=1).tolist()
        earliest=next((i for i in range(19,len(good)) if all(good[i-19:i+1])),None)
        record.update(prefixRows=240 if index==1 else 0,firstPassingWindowEndIndex=earliest)
        write(out/"verified.json",record)
    else:
        proof=area/(label+"_state")/"verified.json"
        if index==1:
            assert read(proof)["prefixRows"]==240
            prefix(history)
        energy(index,folder,out,proof,area/(label+"_observe"),area/(label+"_release"),settings,LIMIT)
        record=read(out/"energy.json")
        record["prefixRows"]=read(proof)["prefixRows"]
        if index==2:
            parent=read(area/"first_energy/accepted.json")
            assert record["energy"]["oldFrameSHA256"]==parent["energy"]["acceptedFrameSHA256"]
            assert record["elapsedPhysicalTimeS"]==2*parent["elapsedPhysicalTimeS"]
        write(out/"accepted.json",record)


def execute(stage,index,area,worker):
    area.relative_to(RUN)
    label="first" if index==1 else "second"
    folder=area/(label+"01")
    out=area/("build01" if stage=="build" else label+"_"+stage)
    if worker:
        assert stage in ("state","energy")
        check(stage,index,folder,out,area)
        return
    out.mkdir(parents=True,exist_ok=False)
    common=[HERE/"PLAN.md",HERE/"run_step.py",Path(__file__).resolve(),PREVIOUS/"check_half.py",
            PREVIOUS/"observe_half.py",PREVIOUS/"run_half.py",ROOT/"scripts/thermistorStartupSurface.py",
            ROOT/"scripts/beadSurfaceResult.py",ROOT/"scripts/runSU2NativeDebugger.py"]
    if stage=="build":
        build(out)
        return
    proof=area/(label+"_state")/"verified.json"
    if stage in ("observe","release"):
        sample(stage,folder,proof,out,area/"build01",settings,LIMIT,common)
        return
    files=common+[input_file("nodes.npz"),input_file("geometry.json"),input_file("body.json"),
                  PREVIOUS/"evidence/history.csv.gz",HERE.parent/"v289/check_fields.py",HERE.parent/"v287/audit_step.py",
                  folder/"debugger_process.json",folder/"case.cfg",folder/f"seed_{index-1:05}.csv",
                  folder/f"restart_{index:05}.csv",folder/f"history_{index:05}.csv"]
    if index==2:
        files.append(area/"first_energy/accepted.json")
    if stage=="energy":
        files.append(proof)
        files.extend(area/(label+"_observe")/n for n in ("before_wall.csv","after_wall.csv","stages.csv","debugger_process.json"))
        files.extend(area/(label+"_release")/n for n in (f"restart_{index:05}.csv","debugger_process.json"))
    cmd=[sys.executable,str(Path(__file__).resolve()),stage,str(index),"--root",str(area),"--worker"]
    rec=observe(cmd,out,{str(p):sha(p) for p in files},120,100*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("build","state","observe","release","energy"))
    cli.add_argument("step",type=int,choices=(1,2))
    cli.add_argument("--root",type=Path,default=RUN)
    cli.add_argument("--worker",action="store_true")
    a=cli.parse_args()
    execute(a.stage,a.step,a.root.resolve(),a.worker)
