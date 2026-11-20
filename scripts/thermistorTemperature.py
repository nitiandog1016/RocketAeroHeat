"""用固定 SU2 原生被动标量求热敏珠的恒物性温度，壁面强制给定温度。"""
import argparse
import os
from pathlib import Path
import uuid

import numpy as np
from readThermistorBead import frame
from thermistorBead import PE, ROOT, read, sha, solve, write_json


def prepare(parent, out):
    previous = frame(parent)
    assert previous["normalTerminal"] and previous["passesResidualAndHeatWindow"]
    out.mkdir(parents=True, exist_ok=False)
    for name in ("mesh.su2","case.json","body.json","geometry.json","nodes.npz"):
        os.link(parent/name,out/name)
    initial = np.loadtxt(parent/"restart.csv",delimiter=",",skiprows=1)
    assert initial.shape[1] == 9 and np.isfinite(initial).all()
    temp = np.zeros(len(initial))
    ids = np.load(parent/"nodes.npz")["ids"]
    m = read(parent/"case.json")["meridianIntervals"]+1
    temp[ids[:,:m]] = 1
    with (parent/"restart.csv").open() as f:
        header = f.readline().strip()+',"Species_0"'
    np.savetxt(out/"initial_flow.csv", np.column_stack((initial,temp)), delimiter=",",comments="",
               header=header,fmt=["%d"]+["%.17g"]*9)
    cfg=(parent/"case.cfg").read_text("ascii")
    cfg=cfg.replace("RMS_VELOCITY-Z, RMS_ENTHALPY)","RMS_VELOCITY-Z, RMS_ENTHALPY, RMS_SPECIES_0)")
    cfg += f"""KIND_SCALAR_MODEL= SPECIES_TRANSPORT
SPECIES_INIT= 0.0
DIFFUSIVITY_MODEL= CONSTANT_DIFFUSIVITY
DIFFUSIVITY_CONSTANT= {1/PE:.17g}
MARKER_WALL_SPECIES= (wall, VALUE, 1.0)
MUSCL_SPECIES= NO
CONV_NUM_METHOD_SPECIES= SCALAR_UPWIND
SPECIES_CLIPPING= NO
"""
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    inputs=[parent/"restart.csv",out/"initial_flow.csv",Path(__file__).resolve(),ROOT/"scripts/thermistorBead.py",
            out/"case.json",out/"body.json",out/"geometry.json",out/"nodes.npz"]
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
        maximumSeconds=360,threads=1,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in inputs]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    out=args.new_output.resolve()
    out.relative_to(ROOT/"outputs")
    prepare(args.parent.resolve(),out)
    solve(out)
