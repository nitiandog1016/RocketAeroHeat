"""固定 0.05 m/s 风洞的独立恒物性动量时间步及原生终帧核对。"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe
from cm011Tunnel import laminar_config

EXE=ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
BASE=ROOT/"outputs/v352_validation/prepared01"
CONDITION=ROOT/"cases/simple_launcher/tunnel.json"
STEPS={"full":.125,"half":.0625}
HEADER=["PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z"]


def read(path):
    return json.loads(path.read_bytes())


def save(path,record):
    path.write_bytes((json.dumps(record,allow_nan=False,indent=2)+"\n").encode())


def config(p,c,dt,step):
    re=c["densityKgM3"]*.05*p["diameterM"]/c["dynamicViscosityPaS"]
    pe=re*c["dynamicViscosityPaS"]*c["specificHeatJkgK"]/c["conductivityWmK"]
    # 关闭原生焓时，HEATFLUX 标记只供无滑移壁识别；热壁边界在独立显热方程中设置。
    changes=dict(AXISYMMETRIC="NO",INC_ENERGY_EQUATION="NO",MUSCL_FLOW="NO",
        CONV_FIELD="(RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_VELOCITY-Z)",
        CONV_RESIDUAL_MINVAL="-9.0",SCREEN_OUTPUT="(TIME_ITER, INNER_ITER, RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_VELOCITY-Z)",
        SCREEN_WRT_FREQ_INNER="20",HISTORY_OUTPUT="(ITER, RMS_RES)",OUTPUT_FILES="(RESTART_ASCII)",
        OUTPUT_WRT_FREQ="1",VOLUME_OUTPUT="(COORDINATES, SOLUTION)")
    lines=[]
    for line in laminar_config(re,pe,600).splitlines():
        key=line.partition("=")[0].strip()
        if key not in ("ITER","MARKER_ISOTHERMAL","MARKER_SYM","WRT_SURFACE_OVERWRITE","WRT_RESTART_OVERWRITE"):
            lines.append(key+"= "+changes.pop(key) if key in changes else line)
    assert not changes
    lines+=["MARKER_HEATFLUX= (wall, 0)","RESTART_SOL= YES","READ_BINARY_RESTART= NO","SOLUTION_FILENAME= seed",
        "RESTART_FILENAME= restart","TIME_DOMAIN= YES","TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER",
        f"TIME_STEP= {dt}",f"RESTART_ITER= {step}",f"TIME_ITER= {step+1}","MAX_TIME= 10","INNER_ITER= 600",
        "HISTORY_WRT_FREQ_INNER= 1","WRT_RESTART_COMPACT= NO"]
    return "\n".join(lines)+"\n"


def inspect(folder, p, seed, branch, size, step, prepared):
    receipt=read(folder/"debugger_process.json")
    good=(receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and not
          (receipt["survivingProcessIDs"] or receipt["supervisorError"] or receipt["stopReason"]) and all(receipt["inputsStable"].values()))
    rec=dict(branch=branch,stepSize=size,step=step,dtDU=STEPS[size],elapsedPhysicalS=step*STEPS[size]*p["diameterM"]/.05,
             seedSHA256=sha(seed),meshSHA256=p["files"]["mesh.su2"],seconds=receipt["elapsedSeconds"],
             peakRSSBytes=receipt["peakFamilyRSSBytes"],exitCode=receipt["debuggerExitCode"],accepted=False)
    if good:
        assert receipt["inputSHA256"][str(seed)]==rec["seedSHA256"]
        assert receipt["inputSHA256"][str(folder/"case.cfg")]==sha(folder/"case.cfg")
        frame=folder/f"restart_{step:05d}.csv"
        with frame.open() as stream:
            assert next(csv.reader(stream))==HEADER
        history=folder/f"history_{step:05d}.csv"
        with history.open() as stream:
            reader=csv.DictReader(stream,skipinitialspace=True);reader.fieldnames=[k.strip() for k in reader.fieldnames]
            records=list(reader)
        keys=["rms[P]","rms[U]","rms[V]","rms[W]"]
        assert set(k for k in reader.fieldnames if k.startswith("rms["))==set(keys)
        assert all(int(row["Time_Iter"])==step for row in records)
        np.testing.assert_array_equal([int(row["Inner_Iter"]) for row in records],np.arange(len(records)))
        norms=np.array([[float(row[k]) for k in keys] for row in records])
        a=np.loadtxt(frame,delimiter=",",skiprows=1)
        np.testing.assert_array_equal(a[:,0],np.arange(p["nodes"]))
        with np.load(prepared/"stations.npz") as m:
            np.testing.assert_allclose(a[:,1:4],m["xyz"],rtol=0,atol=1e-13)
            wall_speed=float(np.max(np.abs(a[m["wall"],5:8])))
        min_p=float(np.min(101325+1.177*.05**2*a[:,4]))
        log=(folder/"debugger.log").read_text(errors="replace")
        assert "No energy equation." in log and "Incompressible Laminar Navier-Stokes' equations." in log and "Exit Success (SU2_CFD)" in log
        rec.update(frameSHA256=sha(frame),historySHA256=sha(history),inner=len(records),maxWallSpeedOverU=wall_speed,
                   minimumAbsolutePressurePa=min_p,last20Maximum=dict(zip(keys,norms[-20:].max(axis=0).tolist())),
                   accepted=bool(np.isfinite(a).all() and len(records)>=20 and np.isfinite(norms).all()
                                 and norms[-20:].max() < -8 and wall_speed<1e-9 and min_p>0))
    return rec


def run(out,branch,size,step,prepared=BASE):
    p=read(prepared/"prepared.json");c=read(CONDITION)
    assert c["velocityMps"]==[-.05,0,0] and c["wallTemperatureK"]==300.001 and p["sectors"]==48
    assert sha(EXE)=="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
    for name in ("mesh.su2","stations.npz",branch+".csv"):
        assert sha(prepared/name)==p["files"][name]
    folder=out/f"{branch}_{size}_{step}"
    seed=folder/f"seed_{step-1:05d}.csv"
    if folder.exists():
        assert (folder/"debugger_process.json").is_file()
        assert (folder/"case.cfg").read_text("ascii")==config(p,c,STEPS[size],step)
    else:
        folder.mkdir(parents=True)
        os.link(prepared/"mesh.su2",folder/"mesh.su2")
        if step==1:
            with (prepared/(branch+".csv")).open("rb") as old,seed.open("xb") as new:
                for line in old:new.write(line.rstrip().rsplit(b",",1)[0]+b"\n")
        else:
            prev=out/f"{branch}_{size}_{step-1}"
            r=read(prev/"step.json");old=prev/f"restart_{step-1:05d}.csv"
            assert r["accepted"] and r["frameSHA256"]==sha(old) and r["meshSHA256"]==p["files"]["mesh.su2"]
            os.link(old,seed)
        cfg=folder/"case.cfg";cfg.write_text(config(p,c,STEPS[size],step),encoding="ascii")
        paths=[EXE,cfg,seed,folder/"mesh.su2",prepared/"prepared.json",CONDITION,DOC/"PLAN.md",Path(__file__).resolve()]
        observe([str(EXE),"case.cfg","--threads","4"],folder,{str(f):sha(f) for f in paths},420,400*1024**2)
    rec=inspect(folder,p,seed,branch,size,step,prepared)
    if not (folder/"step.json").exists():save(folder/"step.json",rec)
    else:assert read(folder/"step.json")==rec
    print(branch,size,step,rec["accepted"],rec.get("last20Maximum"),flush=True)
    return rec["accepted"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    args=parser.parse_args();out=args.out.resolve();out.relative_to(ROOT/"outputs")
    for size,steps in (("full",1),("half",2)):
        for branch in ("control","disturbed"):
            for step in range(1,steps+1):
                if not run(out,branch,size,step):sys.exit(2)
