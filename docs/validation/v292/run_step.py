"""从同一热敏珠旧态运行 260 次复核，第二步只接受核对过的中间帧。"""
import argparse
import gzip
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/"v291"
sys.path.insert(0,str(PREVIOUS))
from run_half import BIN, ENGINE_SHA, MANIFEST, OFFICIAL, ROOT, RUN as PREVIOUS_RUN, config, input_file, observe, read, replace, sha

RUN=ROOT/"outputs/v292_validation"
LIMIT=260


def settings(index):
    return replace(config(index),"INNER_ITER= 240\n",f"INNER_ITER= {LIMIT}\n")


def prefix(history):
    previous=PREVIOUS/"evidence/history.csv.gz"
    prior=gzip.decompress(previous.read_bytes())
    assert len(prior.splitlines())==241 and history.read_bytes().startswith(prior),"前 240 行原生历史不同"
    return previous


def run(index,out,previous):
    ref=read(PREVIOUS/"evidence/first_half.json")
    assert not ref["passed"] and ref["innerRows"]==240 and sha(OFFICIAL)==ENGINE_SHA
    original=PREVIOUS_RUN/"first01/case.cfg"
    assert sha(original)==ref["configSHA256"] and original.read_text("ascii")==config(1)
    locks=[HERE/"PLAN.md",Path(__file__).resolve(),PREVIOUS/"run_half.py",MANIFEST,original,
           PREVIOUS/"evidence/first_half.json",PREVIOUS/"evidence/history.csv.gz",
           ROOT/"scripts/runSU2NativeDebugger.py",input_file("case.cfg"),OFFICIAL]
    if index==1:
        assert previous is None
        old=input_file("seed_00000.csv")
        assert sha(old)==ref["oldFrameSHA256"]
    else:
        assert previous is not None
        previous.relative_to(RUN)
        report=read(previous)
        assert report["accepted"] and report["step"]==1 and report["prefixRows"]==240
        assert abs(report["energy"]["relativeImbalance"])<.005
        old=Path(report["acceptedFrame"])
        old.relative_to(RUN)
        assert sha(old)==report["energy"]["acceptedFrameSHA256"]
        locks.append(previous)
    out.relative_to(RUN)
    out.mkdir(parents=True,exist_ok=False)
    seed=out/f"seed_{index-1:05}.csv"
    os.link(old,seed)
    os.link(input_file("mesh.su2"),out/"mesh.su2")
    (out/"case.cfg").write_text(settings(index),encoding="ascii")
    locks.extend([seed,out/"mesh.su2",out/"case.cfg"])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    rec=observe([str(OFFICIAL),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in locks},420,180*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["survivingProcessIDs"]
    assert "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")
    if index==1:
        prefix(out/"history_00001.csv")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("step",choices=[1,2],type=int)
    cli.add_argument("output",type=Path)
    cli.add_argument("--previous",type=Path)
    a=cli.parse_args()
    run(a.step,a.output.resolve(),a.previous.resolve() if a.previous else None)
