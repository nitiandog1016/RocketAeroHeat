"""用官方 35×25 亚声速平板实测本机 SU2 的 SST-Vm 壁面剪力与残差。"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
BASE=ROOT/"outputs/v359_validation"
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe

EXE=ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
EXE_SHA="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
RESIDUALS=("rms[Rho]","rms[RhoU]","rms[RhoV]","rms[RhoE]","rms[k]","rms[w]")


def read(p):
    return json.loads(p.read_bytes())


def save(p,r):
    p.write_bytes((json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode())


def numbers(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        keys=[k.strip() for k in next(csv.reader(f,skipinitialspace=True))]
        data=np.loadtxt(f,delimiter=",",ndmin=2)
    return {k:v for k,v in zip(keys,data.T)}


def config(p):
    result={}
    for line in p.read_text("ascii").splitlines():
        key,sep,value=line.split("%",1)[0].partition("=")
        if sep:result[key.strip()]=value.strip()
    return result


def execute(folder,c,locks,limit,mib,threads):
    assert sha(EXE)==EXE_SHA
    cfg=folder/"case.cfg"
    cfg.write_bytes(("".join(f"{k}= {v}\n" for k,v in c.items())).encode("ascii"))
    inputs=[cfg,folder/"mesh.su2",EXE,DOC/"PLAN.md",DOC/"sources.json",*locks]
    receipt=observe([str(EXE),cfg.name,"--threads",str(threads)],folder,{str(p):sha(p) for p in inputs},limit,mib*1024**2)
    assert receipt["status"]=="terminal" and all(receipt["inputsStable"].values())
    if receipt["stopReason"] or receipt["debuggerExitCode"] or receipt["survivingProcessIDs"] or receipt["supervisorError"]:
        raise RuntimeError(f"原生进程中止；见 {folder/'debugger_process.json'}")
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")


def plate(folder):
    sources=BASE/"sources";ref=read(DOC/"sources.json")
    assert all(sha(sources/name)==item["sha256"] for name,item in ref.items())
    folder.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(sources/"plate.su2",folder/"mesh.su2")
    c=config(sources/"plate.cfg")
    c.update(MESH_FILENAME="mesh.su2",RESTART_SOL="NO",SST_OPTIONS="(V1994m, VORTICITY)",
             CFL_NUMBER="20",ITER="6000",CONV_RESIDUAL_MINVAL="-10",CONV_STARTITER="60",
             CONV_FIELD="(RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_ENERGY, RMS_TKE, RMS_DISSIPATION)",
             SCREEN_WRT_FREQ_INNER="100",OUTPUT_FILES="(RESTART_ASCII, SURFACE_CSV)",
             OUTPUT_WRT_FREQ="(500, 500)",WRT_RESTART_OVERWRITE="NO",WRT_SURFACE_OVERWRITE="NO",
             VOLUME_OUTPUT="(COORDINATES, SOLUTION, PRIMITIVE, SKIN_FRICTION-X, Y_PLUS, HEAT_FLUX)")
    save(folder/"initial.json",dict(kind="native uniform freestream initialization",restart=False,sourceOptions=c))
    execute(folder,c,[sources/"plate.cfg",sources/"cf.dat",Path(__file__).resolve(),folder/"initial.json"],180,60,1)


def warm(folder,previous):
    proof=read(previous/"checked.json");seed=previous/f"restart_flow_{proof['lastIteration']:06d}.csv"
    assert not proof["accepted"] and sha(seed)==proof["frameSHA256"] and sha(previous/"mesh.su2")==proof["meshSHA256"]
    folder.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(previous/"mesh.su2",folder/"mesh.su2");shutil.copyfile(seed,folder/"seed.csv")
    c=config(previous/"case.cfg")
    c.update(CFL_NUMBER=config(BASE/"sources/plate.cfg")["CFL_NUMBER"],RESTART_SOL="YES",READ_BINARY_RESTART="NO",
             SOLUTION_FILENAME="seed",ITER="10000",OUTPUT_WRT_FREQ="(1000,1000)",WRT_RESTART_COMPACT="NO",
             OUTPUT_FILES="(RESTART_ASCII, SURFACE_PARAVIEW_ASCII)",VOLUME_OUTPUT="(COORDINATES, SOLUTION, PRIMITIVE)")
    execute(folder,c,[folder/"seed.csv",previous/"checked.json",DOC/"RECOVERY.md",Path(__file__).resolve()],180,60,1)


def audit(folder):
    hist=numbers(folder/"history.csv")
    iterations=hist["Inner_Iter"].astype(int)
    np.testing.assert_array_equal(iterations,np.arange(len(iterations)))
    assert len(iterations)>=20 and set(k for k in hist if k.startswith("rms["))==set(RESIDUALS)
    end=int(iterations[-1]);frame=folder/f"restart_flow_{end:06d}.csv"
    data=numbers(frame)
    # 原生结构网在恒定 x 处用壁面与第一排速度复核剪力。零壁速平板在此处 d(v)/dx=0。
    ids=np.arange(6,35);first=ids+35
    np.testing.assert_allclose(data["x"][ids],data["x"][first],atol=1e-14,rtol=0)
    np.testing.assert_array_equal(data["y"][ids],np.zeros(len(ids)))
    rho=data["Density"];u=data["Momentum_x"]/rho;v=data["Momentum_y"]/rho
    R=287.058;gamma=1.4
    temp=(data["Energy"]/rho-(u*u+v*v)/2-data["Turb_Kin_Energy"])*(gamma-1)/R
    def mu(t):
        return 1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
    speed=.2*np.sqrt(gamma*R*300);rho_inf=5e6*mu(300)/speed
    shear=mu(temp[ids])*(u[first]-u[ids])/(data["y"][first]-data["y"][ids])
    cf_nodes=shear/(.5*rho_inf*speed**2)
    x=data["x"][ids];station=.970084071
    gradient_cf=float(np.interp(station,x,cf_nodes))
    yplus=(data["y"][first]-data["y"][ids])*np.sqrt(abs(shear)*rho[ids])/mu(temp[ids])
    wall=folder/f"surface_flow_{end:06d}.vtk"
    native=wall.is_file()
    if native:
        from auditSU2SSTState import vtk_tokens
        xyz,fields=vtk_tokens(wall);order=np.argsort(xyz[:,0])
        np.testing.assert_allclose(xyz[order,0],x,atol=6e-6,rtol=0)
        cf=float(np.interp(station,xyz[order,0],fields["Skin_Friction_Coefficient"][order,0]))
        yplus=fields["Y_Plus"][:,0]
    else:
        wall=folder/f"surface_flow_{end:06d}.csv";surface=numbers(wall)
        np.testing.assert_array_equal(surface["PointID"],data["PointID"][ids])
        cf=gradient_cf
    assert len(x)==29 and np.all(np.diff(x)>0) and x[0]<=.970084071<=x[-1]
    # 只读独立参考文件中两个代码的同一 816 单元行。
    refs={};name=None
    for line in (BASE/"sources/cf.dat").read_text().splitlines():
        if line.startswith("zone"):name=line.split('"')[1]
        elif line.split() and line.split()[0]=="816.":refs[name]=float(line.split()[3])
    assert set(refs)=={"FUN3D","CFL3D"}
    maxima={k:float(max(hist[k][-20:])) for k in RESIDUALS}
    errors={k:abs(cf/v-1) for k,v in refs.items()}
    rc=read(folder/"debugger_process.json")
    gates=dict(finite=all(np.isfinite(data[k]).all() for k in ("Density","Momentum_x","Momentum_y","Energy","Omega")),
               physicalState=bool(rho.min()>0 and temp.min()>0 and data["Omega"].min()>0),
               nativeFrictionAvailable=native,sixResiduals=max(maxima.values())< -8,
               friction=max(errors.values())<.05,nativeSuccess=rc["debuggerExitCode"]==0 and rc["stopReason"] is None)
    proof=dict(accepted=all(gates.values()),gates=gates,requestedMach=.2,gridCells=816,wallNodes=len(x),
               stationX=.970084071,stationNodes=x[np.argsort(abs(x-.970084071))[:2]].tolist(),
               skinFriction=cf,stationGradientEstimate=gradient_cf,references=refs,relativeErrors=errors,last20Maximum=maxima,lastIteration=end,
               minimumWallYPlus=float(min(yplus)),maximumWallYPlus=float(max(yplus)),nativeWallRead=native,
               thermalBoundary="adiabatic",heatTransferValidation=False,seconds=rc["elapsedSeconds"],
               peakRSSBytes=rc["peakFamilyRSSBytes"],sourceSHA256=sha(folder/"case.cfg"),meshSHA256=sha(folder/"mesh.su2"),
               historySHA256=sha(folder/"history.csv"),frameSHA256=sha(frame),wallSHA256=sha(wall))
    save(folder/"checked.json",proof)
    print(json.dumps(proof,indent=2))
    return proof


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path);p.add_argument("stage",choices=("solve","read","warm"));args=p.parse_args()
    folder=args.folder.resolve();folder.relative_to(BASE)
    if args.stage=="solve":plate(folder)
    elif args.stage=="warm":warm(folder,BASE/"plate01")
    else:audit(folder)
