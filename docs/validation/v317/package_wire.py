"""核对长引线正常终帧后打包原面表和来源，保留体场在本机。"""
import argparse
from pathlib import Path
import sys
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha,write_json


def package(run,out):
    done=read(run/"completed.json");stage=run/done["stage"];summary=read(run/"surface/summary.json")
    assert done["frame"]==summary["sourceFrameSHA256"]==sha(stage/"restart.csv")
    assert done["summarySHA256"]==sha(run/"surface/summary.json") and summary["wallSHA256"]==done["wallSHA256"]==sha(run/"surface/wall.csv")
    assert summary["reactionSHA256"]==sha(stage/"reaction0/restart.csv") and summary["boundarySHA256"]==sha(run/"prepared/profile.dat")
    for name in ("reaction0","plus20","reaction20","step_control","step_observed"):
        r=read(stage/name/"debugger_process.json");assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and not (r["stopReason"] or r["survivingProcessIDs"])
    assert max(summary["window"].values())<-8 and max(summary["plus20Window"].values())<-8
    assert abs(summary["energy"]["relativeDifference"])<.005 and summary["relativeChange20"]<.002
    out.mkdir(parents=True,exist_ok=True);archive=out/"longWire12.zip"
    with ZipFile(archive,"x") as z:
        for n in ("wall.csv","summary.json"):
            entry=ZipInfo(n);entry.compress_type=ZIP_DEFLATED;z.writestr(entry,(run/"surface"/n).read_bytes())
    write_json(out/"longWire12.json",dict(archiveSHA256=sha(archive),summarySHA256=done["summarySHA256"],frameSHA256=done["frame"],
               wallSHA256=done["wallSHA256"],boundarySHA256=summary["boundarySHA256"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("run",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();package(a.run.resolve(),a.out.resolve())
