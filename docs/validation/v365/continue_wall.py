"""根据真实原生终态在独立目录重读首个更新或从末帧有限续算。"""
import argparse
import json
from pathlib import Path
import shutil
from run_wall import BASE, DOC, EXE, observe, read, sha


def advance(previous,folder,stage):
    previous=previous.resolve()
    folder=folder.resolve()
    for p in (previous,folder):
        p.relative_to(BASE)
    receipt=read(previous/"debugger_process.json")
    assert receipt["debuggerExitCode"]==0 and not receipt["survivingProcessIDs"] and not receipt["stopReason"]
    assert all(sha(Path(k))==v for k,v in receipt["inputSHA256"].items())
    folder.mkdir(parents=True,exist_ok=False)
    for name in ("mesh.su2","condition.json"):
        shutil.copyfile(previous/name,folder/name)
    lines=(previous/"case.cfg").read_text().splitlines()
    cfg=dict((key.strip(),v.strip()) for key,v in (line.split("=",1) for line in lines))
    assert cfg["ITER"]=="2500" and cfg["RESTART_SOL"]=="NO"
    if stage=="continue":
        seed=previous/"restart_002499.csv"
        assert sha(seed)==sha(previous/"restart.csv")
        shutil.copyfile(seed,folder/"seed.csv")
        cfg.update(RESTART_SOL="YES",READ_BINARY_RESTART="NO",SOLUTION_FILENAME="seed",
                   RAMP_MUSCL="NO",ITER="10000",OUTPUT_WRT_FREQ="(1000,1000)")
        seconds,cap=300,160
        initial=folder/"seed.csv"
    else:
        cfg.update(ITER="1",OUTPUT_WRT_FREQ="(1,1)")
        seconds,cap=15,24
        initial=previous/"initial.npz"
    path=folder/"case.cfg"
    path.write_text("".join(f"{key}= {v}\n" for key,v in cfg.items()),encoding="ascii")
    locks=dict(receipt["inputSHA256"])
    for p in (path,folder/"mesh.su2",folder/"condition.json",initial,DOC/"CONTINUE.md",Path(__file__)):
        locks[str(p.resolve())]=sha(p)
    result=observe([str(EXE),path.name,"--threads","2"],folder,locks,seconds,cap*1024**2)
    print(json.dumps({k:result[k] for k in ("elapsedSeconds","debuggerExitCode","stopReason","survivingProcessIDs")}),flush=True)
    assert result["debuggerExitCode"]==0 and not result["stopReason"] and "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("previous",type=Path)
    p.add_argument("folder",type=Path)
    p.add_argument("stage",choices=("continue","first"))
    a=p.parse_args()
    advance(a.previous,a.folder,a.stage)
