"""固定风洞工况的一步全周可行性试算，保留程序终态和未通过的残差。"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Tunnel import laminar_config
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe


def run(prepared,out,branch,scheme,control=None):
    prepared.relative_to(ROOT/"outputs");out.relative_to(ROOT/"outputs")
    p=json.loads((prepared/"prepared.json").read_bytes())
    for name,digest in p["files"].items():
        assert sha(prepared/name)==digest
    if branch=="disturbed":
        if control is None:
            raise ValueError("需指定同一时间步及空间离散下已通过的无扰动对照目录。")
        reference=json.loads((control/"step.json").read_bytes())
        if not (reference["pressureAccepted"] and reference["branch"]=="control" and reference["scheme"]==scheme
                and reference["meshSHA256"]==p["files"]["mesh.su2"] and reference["seedSHA256"]==p["files"]["control.csv"]
                and reference["nativeFrameSHA256"]==sha(control/"restart_00001.csv")
                and reference["historySHA256"]==sha(control/"history_00001.csv")):
            raise ValueError("无扰动步未通过或输入不一致，无法构成配对试算。")
    case=json.loads((ROOT/"cases/simple_launcher/tunnel.json").read_bytes())
    assert case["velocityMps"]==[-.05,0,0] and case["wallTemperatureK"]==300.001
    rho=case["densityKgM3"];u=.05;d=p["diameterM"]
    re=rho*u*d/case["dynamicViscosityPaS"]
    pe=rho*u*d*case["specificHeatJkgK"]/case["conductivityWmK"]
    out.mkdir(parents=True,exist_ok=False)
    for name,target in [("mesh.su2","mesh.su2"),(branch+".csv","seed_00000.csv")]:
        os.link(prepared/name,out/target)
    settings=dict(AXISYMMETRIC="NO",CONV_RESIDUAL_MINVAL="-9.5",
        CONV_FIELD="(RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY)",
        SCREEN_OUTPUT="(TIME_ITER, INNER_ITER, RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY)",
        SCREEN_WRT_FREQ_INNER="10",OUTPUT_FILES="(RESTART_ASCII)",VOLUME_OUTPUT="(COORDINATES, SOLUTION)",
        OUTPUT_WRT_FREQ="1",RESTART_SOL="YES",READ_BINARY_RESTART="NO",SOLUTION_FILENAME="seed",
        RESTART_FILENAME="restart",TIME_DOMAIN="YES",TIME_MARCHING="DUAL_TIME_STEPPING-1ST_ORDER",
        TIME_STEP="0.125",MAX_TIME="10.0",RESTART_ITER="1",TIME_ITER="2",INNER_ITER="360",
        HISTORY_WRT_FREQ_INNER="1",WRT_RESTART_COMPACT="NO",
        MUSCL_FLOW="NO" if scheme=="firstOrder" else "YES")
    lines=[]
    for line in laminar_config(re,pe,360).splitlines():
        key=line.partition("=")[0].strip()
        if key not in ("ITER","MARKER_SYM","WRT_SURFACE_OVERWRITE","WRT_RESTART_OVERWRITE"):
            lines.append(key+"= "+settings.pop(key) if key in settings else line)
    lines.extend(f"{k}= {v}" for k,v in settings.items())
    cfg=out/"case.cfg";cfg.write_text("\n".join(lines)+"\n",encoding="ascii")
    exe=ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
    assert sha(exe)=="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
    inputs=[exe,cfg,out/"mesh.su2",out/"seed_00000.csv",prepared/"prepared.json",HERE/"PLAN.md",
            Path(__file__).resolve(), ROOT/"scripts/runSU2NativeDebugger.py", ROOT/"cases/simple_launcher/tunnel.json"]
    r=observe([str(exe),"case.cfg","--threads","4"],out,{str(f):sha(f) for f in inputs},360,400*1024**2)
    record=dict(branch=branch,scheme=scheme,seedSHA256=p["files"][branch+".csv"],meshSHA256=p["files"]["mesh.su2"],
                processSeconds=r["elapsedSeconds"],peakRSSBytes=r["peakFamilyRSSBytes"],exitCode=r["debuggerExitCode"],
                stopReason=r["stopReason"],dtSeconds=.125*d/u,dtOverBaseConvectiveTime=.125*d/p["baseDiameterM"],
                pressureAccepted=False,heatAccepted=False)
    if r["debuggerExitCode"]==0 and not (r["supervisorError"] or r["survivingProcessIDs"] or r["stopReason"]):
        with (out/"history_00001.csv").open() as f:
            reader=csv.DictReader(f,skipinitialspace=True)
            reader.fieldnames=[k.strip() for k in reader.fieldnames]
            rows=list(reader)
        assert all(int(row["Time_Iter"])==1 for row in rows)
        np.testing.assert_array_equal([int(row["Inner_Iter"]) for row in rows],np.arange(len(rows)))
        keys=["rms[P]","rms[U]","rms[V]","rms[W]","rms[h]"]
        residual=np.array([[float(row[k]) for k in keys] for row in rows])
        assert np.isfinite(residual).all()
        frame=out/"restart_00001.csv";data=np.loadtxt(frame,delimiter=",",skiprows=1)
        with np.load(prepared/"stations.npz") as m:
            np.testing.assert_allclose(data[:,1:4],m["xyz"],rtol=0,atol=1e-13)
            np.testing.assert_array_equal(data[:,0],np.arange(p["nodes"]))
            finite=bool(np.isfinite(data).all())
            wall_max=float(np.max(np.abs(data[m["wall"],5:8])))
            min_pressure=float(np.min(case["ambientPressurePa"]+rho*u*u*data[:,4]))
            base=data[m["base"],4]*rho*u**2*1000
            radii=np.linalg.norm(m["xyz"][m["base"][0],1:],axis=1)
        assert "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")
        record.update(nativeFrameSHA256=sha(frame),historySHA256=sha(out/"history_00001.csv"),inner=len(rows),
                      finalWindowMax=dict(zip(keys,np.max(residual[-20:],axis=0).tolist())),
                      finite=finite,maxWallSpeedOverU=wall_max,minAbsolutePressurePa=min_pressure,
                      pressureAccepted=bool(finite and wall_max<1e-9 and min_pressure>0 and
                                            len(rows)>=20 and np.max(residual[-20:]) < -8))
        if record["pressureAccepted"]:
            radius=radii/radii[-1]
            sample=np.array([np.interp([0,.286,.572,.858],radius,row) for row in base])
            # 双线性压力在平底三角扇片内采用径向面积权重，周向相邻两侧各占一半。
            # 底心是所有方位共用的同一压力自由度。
            a,b=radius[:-1],radius[1:]
            ray_mean=np.sum((b-a)*((2*a+b)*base[:,:-1]+(a+2*b)*base[:,1:]),axis=1)/3
            record.update(tapMilliPa=sample.tolist(),baseWedgeMilliPa=(.5*(ray_mean+np.roll(ray_mean,-1))).tolist())
    (out/"step.json").write_text(json.dumps(record,indent=2,allow_nan=False),encoding="utf8")
    print(json.dumps({k:record[k] for k in ("branch","processSeconds","exitCode","pressureAccepted")}))
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("branch",choices=("control","disturbed"))
    parser.add_argument("--scheme",choices=("muscl","firstOrder"),default="muscl")
    parser.add_argument("--control",type=Path,help="扰动支的无扰动对照目录")
    args=parser.parse_args()
    result=run(args.prepared.resolve(),args.output.resolve(),args.branch,args.scheme,args.control)
    sys.exit(0 if result["pressureAccepted"] else 2)
