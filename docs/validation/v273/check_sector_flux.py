"""复算低速三维扇区各真实边界的显热收支，验证对称面和内部面的贡献。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Tunnel import CP,D,DT,RHO,U,write_json
from readCM011TunnelSector import state_at_frame
from runEmbeddedSU2 import sha


def audit(parent,control,observed,frame_reader=state_at_frame,meridian_nodes=97,copies=60,
          symmetry=("sym0","sym1"),schema="cm011-sector-energy-v1"):
    frame=frame_reader(parent)
    assert frame["normalTerminal"] and frame["passesResidualAndHeatWindow"]
    source=sha(parent/"restart.csv")
    for folder in (control,observed):
        sample=json.loads((folder/"sample.json").read_bytes())
        step=json.loads((folder/"step_process.json").read_bytes())
        assert step["exitCode"]==0 and step["stopReason"] is None and step["livePids"]==[]
        assert source==sample["parentRestartSHA256"]==sha(folder/"initial_flow.csv")
        assert all(sha(Path(n))==v for n,v in sample["inputSHA256"].items())
        assert all(sha(folder/n)==h for n,h in sample["outputSHA256"].items())
    assert (control/"case.cfg").read_bytes()==(observed/"case.cfg").read_bytes()
    assert sha(control/"mesh.su2")==sha(observed/"mesh.su2")==sha(parent/"mesh.su2")
    def restart(path):
        a=np.loadtxt(path,delimiter=",",skiprows=1)
        assert np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(len(a)))
        return a
    initial=restart(parent/"restart.csv")
    traced=np.loadtxt(observed/"observed_input.csv",delimiter=",",skiprows=1)
    traced=traced[np.argsort(traced[:,0])]
    assert np.array_equal(traced[:,0],initial[:,0]) and np.array_equal(initial[:,4:],traced[:,1:])
    updated=restart(observed/"restart.csv")
    diff=float(np.max(abs(updated[:,4:]-restart(control/"restart.csv")[:,4:])))
    assert diff<2e-7
    ref=np.loadtxt(observed/"far_reference.csv",delimiter=",",skiprows=1)
    grid=np.load(parent/"nodes.npz")
    assert np.array_equal(np.sort(ref[:,0]),np.unique(grid["ids"][:,-meridian_nodes:]))
    href=ref[0,1]
    assert np.all(ref[:,1]==href) and abs(href+297.15)<1e-12
    with (observed/"sector_stages.csv").open() as f:
        rows=list(csv.DictReader(f))
    names=[r["stage"] for r in rows]
    assert names==["start","internal","before_far","after_far","before_wall","after_wall",*symmetry,"end"]
    values=np.array([[float(r["mass"]),float(r["enthalpy"])] for r in rows])
    assert np.isfinite(values).all()
    deltas=np.diff(values,axis=0,prepend=np.zeros((1,2)))
    mass_scale=copies*RHO*U*D**2
    power_scale=mass_scale*CP*DT
    increments={name:dict(massKgS=m*mass_scale,rawEnthalpyW=e*power_scale,referenceAdvectionW=m*href*power_scale,
                         excessHeatW=(e-m*href)*power_scale) for name,(m,e) in zip(names,deltas)}
    for key in ("start","internal","before_far","before_wall",*symmetry,"end"):
        assert abs(increments[key]["excessHeatW"])<1e-12 and abs(increments[key]["massKgS"])<1e-12,key
    far=increments["after_far"]["excessHeatW"]
    wall=-increments["after_wall"]["excessHeatW"]
    assert far>0 and wall>0
    m,e=values[-1]
    residual=(e-href*m)*power_scale
    assert abs(residual-(far-wall))<1e-12
    return dict(schema=schema,sourceRestartSHA256=source,savedIndex=frame["savedIndex"],
                sourceFrameSHA256=frame["frameSHA256"],copies=copies,referenceNativeEnthalpy=href,
                sampledInputExact=True,updatedStateDifference=diff,allBoundaryIncrements=increments,
                nativeWallReleaseW=wall,nativeFarExcessOutW=far,wallMinusFarW=wall-far,
                wallMinusFarFraction=(wall-far)/wall,weightedResidualW=residual,
                belowHalfPercent=bool(abs(wall-far)<.005*wall),physicallyValidated=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("parent","official","observed","new_report"): parser.add_argument(name,type=Path)
    args=parser.parse_args()
    write_json(args.new_report,audit(args.parent.resolve(),args.official.resolve(),args.observed.resolve()))
