"""在固定的 12 K 量纲下独立改变 SU2 原生强标量边界，从冷标量计算并读取壁反力。"""
import os
from pathlib import Path
import re
import sys
import uuid

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"docs/validation/v272")]
from thermistorBead import ENGINE,ROOT,read,sha,solve,write_json
from thermistorReaction import observe,terminal
from observe_flux import invoke

BASE=ROOT/"outputs/v282_validation/strong_outer01"
OBSERVER=ROOT/"outputs/v282_validation/scalar_observer01/thermistor_observe.exe"


def run(folder,value):
    terminal(BASE)
    assert sha(BASE/"initial_flow.csv")=="27f0f97031d10b0f1e39a60ccac289ce25f7716b276a0e88b33e6b6bb23e37f6"
    assert sha(BASE/"case.cfg")=="b61456e3fc30a5161d4ff74717fb100d1b8418ac2822e027d2975e447d33c6f7"
    assert sha(OBSERVER)=="62f95580241e968df972e1ddaec71727913aff3811e8eb781773cc2c90c25135"
    strong=folder/"strong";strong.mkdir(parents=True,exist_ok=False)
    os.link(BASE/"mesh.su2",strong/"mesh.su2")
    initial=np.loadtxt(BASE/"initial_flow.csv",delimiter=",",skiprows=1)
    with np.load(BASE/"nodes.npz") as grid:
        wall=np.unique(grid["ids"][:,:65])
    assert np.array_equal(np.flatnonzero(initial[:,-1]),wall)
    initial[wall,-1]=value
    head=(BASE/"initial_flow.csv").open()
    with head: header=head.readline().strip()
    np.savetxt(strong/"initial_flow.csv",initial,delimiter=",",comments="",header=header,fmt=["%d"]+["%.17g"]*9)
    line="MARKER_WALL_SPECIES= (wall, VALUE, 1.0)"
    cfg=(BASE/"case.cfg").read_text("ascii");assert cfg.count(line)==1
    (strong/"case.cfg").write_text(cfg.replace(line,f"MARKER_WALL_SPECIES= (wall, VALUE, {value})"),encoding="ascii")
    inputs=[BASE/"initial_flow.csv",BASE/"case.cfg",BASE/"nodes.npz",strong/"initial_flow.csv",Path(__file__).resolve()]
    write_json(strong/"engine_request.json",dict(workDirectory=str(strong),configFile="case.cfg",meshFile="mesh.su2",
        maximumSeconds=360,threads=1,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in inputs]))
    solve(strong);terminal(strong)
    observe(strong,folder/"reaction0",release=True,wall_value=value)
    observe(strong,folder/"plus20",release=False,wall_value=value)
    observe(folder/"plus20",folder/"reaction20",release=True,wall_value=value)
    for tag,exe,enable in (("observed",OBSERVER,True),("control",ENGINE,False)):
        step=folder/tag;step.mkdir()
        os.link(strong/"mesh.su2",step/"mesh.su2");os.link(strong/"restart.csv",step/"initial_flow.csv")
        cfg,n=re.subn(r"^ITER= \d+$","ITER= 1",(strong/"case.cfg").read_text("ascii"),flags=re.M);assert n==1
        (step/"case.cfg").write_text(cfg,encoding="ascii")
        invoke("step",[exe,"case.cfg","--threads","1"],step,step,60,1,trace=enable)


if __name__=="__main__":
    dest=Path(sys.argv[1]).resolve();dest.relative_to(ROOT/"outputs")
    for name,value in (("half",0.5),("zero",0.0)):
        run(dest/name,value)
