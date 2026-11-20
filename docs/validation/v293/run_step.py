"""在通过嵌套网格与初场验收的热敏珠试件上推进一个物理半步。"""
import argparse
import os
from pathlib import Path
import sys

from prepare import HERE, ROOT, RUN, SOURCE, bound_source, read, sha, observe
sys.path.insert(0,str(HERE.parent/"v290"))
from observe_heat import BIN,ENGINE_SHA,OFFICIAL,replace

LIMIT=560


def configuration(index):
    base=ROOT/("outputs/v292_validation/first01" if index==1 else "outputs/v292_validation/pair/second01")/"case.cfg"
    manifest=read(HERE.parent/"v292/evidence/local_assets.json")
    assert sha(base)==manifest[base.relative_to(ROOT).as_posix()]["sha256"]
    return replace(base.read_text("ascii"),"INNER_ITER= 260\n",f"INNER_ITER= {LIMIT}\n")


def inputs(grid):
    info=read(grid/"grid.json")
    receipt=read(grid/"debugger_process.json")
    assert info["admitted"] and receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and not receipt["stopReason"]
    assert not receipt["supervisorError"] and not receipt["survivingProcessIDs"] and all(receipt["inputsStable"].values())
    for file,key in (("mesh.su2","gridSHA256"),("nodes.npz","nodeMapSHA256"),("seed_00000.csv","seedSHA256")):
        assert sha(grid/file)==info[key]
    return info


def run(index,grid,out,previous):
    info=inputs(grid)
    locks=[HERE/"PLAN.md",Path(__file__).resolve(),grid/"grid.json",grid/"debugger_process.json",OFFICIAL,
           ROOT/"scripts/runSU2NativeDebugger.py",HERE.parent/"v292/evidence/local_assets.json"]
    assert sha(OFFICIAL)==ENGINE_SHA and info["nodes"]==172482
    if index==1:
        assert previous is None
        old=grid/"seed_00000.csv"
    else:
        previous.relative_to(RUN)
        r=read(previous)
        assert r["step"]==1 and r["accepted"] and r["gridSHA256"]==info["gridSHA256"]
        old=Path(r["acceptedFrame"]);old.relative_to(RUN)
        assert sha(old)==r["energy"]["acceptedFrameSHA256"] and abs(r["energy"]["relativeImbalance"])<.005
        locks.append(previous)
    out.relative_to(RUN);out.mkdir(parents=True,exist_ok=False)
    seed=out/f"seed_{index-1:05}.csv"
    os.link(old,seed);os.link(grid/"mesh.su2",out/"mesh.su2")
    cfg=out/"case.cfg";cfg.write_text(configuration(index),encoding="ascii")
    locks.extend([cfg,seed,out/"mesh.su2"])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    rec=observe([str(OFFICIAL),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in locks},900,220*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["survivingProcessIDs"]
    assert "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("step",choices=[1,2],type=int);cli.add_argument("output",type=Path)
    cli.add_argument("--grid",type=Path,default=RUN/"grid01");cli.add_argument("--previous",type=Path)
    a=cli.parse_args();run(a.step,a.grid.resolve(),a.output.resolve(),a.previous.resolve() if a.previous else None)
