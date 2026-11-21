"""在相同带线体网上独立计算 6 K 温差，保留强标量的 12 K 量纲与同阶段壁面反力。"""
import argparse
import os
import re
import sys
import uuid
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v282")]
from beadSurfaceResult import thermal_frame
from observe_bead import check as check_energy
from observe_flux import invoke
from thermistorBead import ENGINE, ENGINE_SHA, read, sha, solve, write_json
from thermistorReaction import observe, terminal

BASE=ROOT/"outputs/v300_validation/wired01"


def prepare(folder):
    rec=terminal(BASE)
    reference=read(HERE.parent/"v300/evidence/summary.json")
    inputs=read(HERE.parent/"v300/evidence/inputs.json")["inputs"]
    assert all(sha(BASE/k)==value for k,value in inputs.items())
    assert sha(ENGINE)==ENGINE_SHA and sha(BASE/"restart.csv")==reference["sourceFrameSHA256"]
    observer=BASE/"ledger/thermistor_observe.exe"
    assert sha(observer)==reference["energy"]["observerSHA256"]
    folder.mkdir(parents=True,exist_ok=False)
    for name in ("mesh.su2","nodes.npz","body.json","geometry.json"):
        os.link(BASE/name,folder/name)
    cfg=(BASE/"case.cfg").read_text("ascii")
    for before,after in (("MARKER_WALL_SPECIES= (wall, VALUE, 1.0)","MARKER_WALL_SPECIES= (wall, VALUE, 0.5)"),
                         ("MARKER_ISOTHERMAL= (wall, 2)","MARKER_ISOTHERMAL= (wall, 1.5)")):
        assert cfg.count(before)==1
        cfg=cfg.replace(before,after)
    (folder/"case.cfg").write_text(cfg,"ascii")
    with np.load(folder/"nodes.npz") as grid:walls=np.unique(grid["wallMap"])
    initial=np.loadtxt(BASE/"initial_flow.csv",delimiter=",",skiprows=1)
    assert np.array_equal(np.flatnonzero(initial[:,9]),walls) and np.all(initial[walls,9]==1)
    initial[walls,9]=.5
    assert np.all(initial[walls,8]==2-298.15)
    initial[walls,8]=1.5-298.15
    with (BASE/"initial_flow.csv").open() as f:header=f.readline().strip()
    np.savetxt(folder/"initial_flow.csv",initial,delimiter=",",comments="",fmt=["%d"]+["%.17g"]*9,header=header)
    spec=read(BASE/"case.json");spec.update(wallK=299.15,scalarReferenceDeltaK=12.)
    write_json(folder/"case.json",spec)
    files=[HERE/"PLAN.md",Path(__file__).resolve(),ROOT/"scripts/thermistorReaction.py",HERE.parent/"v282/observe_bead.py",
           BASE/"case.cfg",BASE/"initial_flow.csv",*(folder/name for name in ("geometry.json","case.json","nodes.npz","body.json","initial_flow.csv"))]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               threads=1,maximumSeconds=720,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in files]))
    write_json(folder/"inputs.json",dict(referenceFrame=reference["sourceFrameSHA256"],referenceTable=reference["wallSHA256"],
            executable=rec["inputSHA256"][str(ENGINE)],observerSHA256=sha(observer),case=spec,
            hashes={name:sha(folder/name) for name in ("mesh.su2","case.cfg","initial_flow.csv","nodes.npz","body.json","case.json","geometry.json","engine_request.json")}))


def probe(folder):
    thermal_frame(folder)
    observe(folder,folder/"reaction0",release=True,wall_value=.5)
    observe(folder,folder/"plus20",release=False,wall_value=.5)
    observe(folder/"plus20",folder/"reaction20",release=True,wall_value=.5)
    ledger=folder/"ledger";ledger.mkdir()
    instrument=ledger/"thermistor_observe.exe"
    os.link(BASE/"ledger"/instrument.name,instrument)
    assert sha(instrument)==read(folder/"inputs.json")["observerSHA256"]
    cfg,n=re.subn(r"^ITER= \d+$","ITER= 1",(folder/"case.cfg").read_text("ascii"),flags=re.M)
    assert n==1
    for tag,exe,trace in (("control",ENGINE,False),("observed",instrument,True)):
        step=ledger/tag;step.mkdir()
        os.link(folder/"mesh.su2",step/"mesh.su2");os.link(folder/"restart.csv",step/"initial_flow.csv")
        (step/"case.cfg").write_text(cfg,"ascii")
        invoke("step",[exe,"case.cfg","--threads","1"],step,step,60,1,trace=trace)
    check_energy(folder,folder,ledger)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("action",choices=("prepare","solve","probes"));p.add_argument("folder",type=Path)
    a=p.parse_args();f=a.folder.resolve();f.relative_to(ROOT/"outputs")
    if a.action=="prepare":prepare(f)
    elif a.action=="solve":
        assert all(sha(f/k)==v for k,v in read(f/"inputs.json")["hashes"].items())
        solve(f);terminal(f)
    else:probe(f)
