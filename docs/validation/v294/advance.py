"""在已验收的热敏珠 56 层体网上按 v294 限额推进，逐项绑定源帧、网格与原生历史。"""
import argparse
import gzip
import hashlib
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/"v293"
sys.path.insert(0,str(PREVIOUS))
from run_step import ROOT, OFFICIAL, ENGINE_SHA, BIN, observe, read, replace, sha
from thermistorStartupCache import BUILD, CASE, GRID, ensure_grid, grid_info
from thermistorSeedRecovery import identical

RUN=ROOT/"outputs/v294_validation"
LIMIT=800


def settings(index):
    cfg=identical(CASE/"halfstep.cfg","72c487e7378c25db6a9b57ebf51fdc0b5b8658df920a01dcce667a84224658be").read_text("ascii")
    if index==1:return cfg
    assert index==2
    return replace(replace(cfg,"RESTART_ITER= 1\n","RESTART_ITER= 2\n"),"TIME_ITER= 2\n","TIME_ITER= 3\n")


def prefix(history):
    record=read(PREVIOUS/"evidence/first_state.json")
    saved=gzip.decompress((PREVIOUS/"evidence/first_history.csv.gz").read_bytes())
    assert not record["passed"] and hashlib.sha256(saved).hexdigest()==record["historySHA256"]
    assert len(saved.splitlines())==561
    assert history.read_bytes().startswith(saved),"首半步前 560 行发生变化"


def run(index,area):
    g=ensure_grid()
    assert sha(OFFICIAL)==ENGINE_SHA
    files=[HERE/"PLAN.md",Path(__file__).resolve(),PREVIOUS/"run_step.py",GRID/"grid.json",CASE/"halfstep.cfg",
           ROOT/"scripts/thermistorStartupCache.py",ROOT/"scripts/thermistorSeedRecovery.py",
           PREVIOUS/"evidence/grid.json",PREVIOUS/"evidence/first_history.csv.gz",OFFICIAL,ROOT/"scripts/runSU2NativeDebugger.py"]
    if index==1:
        old=GRID/"seed_00000.csv"
        assert sha(old)==g["seedSHA256"]
    else:
        proof=area/"first_energy/accepted.json";r=read(proof)
        assert r["accepted"] and r["step"]==1 and r["gridSHA256"]==g["gridSHA256"] and r["verifiedPrefixRows"]==560
        old=Path(r["acceptedFrame"]);old.relative_to(area)
        assert sha(old)==r["energy"]["acceptedFrameSHA256"] and abs(r["energy"]["relativeImbalance"])<.005
        files.append(proof)
    out=area/("first01" if index==1 else "second01")
    out.mkdir(parents=True,exist_ok=False)
    seed=out/f"seed_{index-1:05}.csv";cfg=out/"case.cfg"
    os.link(old,seed);os.link(GRID/"mesh.su2",out/"mesh.su2");cfg.write_text(settings(index),encoding="ascii")
    files.extend([seed,cfg,out/"mesh.su2"])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    rec=observe([str(OFFICIAL),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in files},1100,220*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["survivingProcessIDs"]
    assert "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")
    if index==1:prefix(out/"history_00001.csv")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("step",type=int,choices=(1,2));cli.add_argument("--root",type=Path,default=RUN)
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);run(a.step,area)
