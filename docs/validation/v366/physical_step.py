"""从冻结细网末帧按物理时间推进，逐步审核原生内迭代、终帧和可传递的初场。"""
import argparse
import hashlib
from pathlib import Path
import shutil
import sys
import zipfile

import numpy as np

DOC = Path(__file__).resolve().parent
ROOT = DOC.parents[2]
BASE = ROOT/"outputs/v366_validation"
OLD = ROOT/"outputs/v365_validation/warm2"
sys.path.insert(0, str(DOC.parent/"v365"))
from read_wall import RES, numbers, read, save, sha, terminal, vtk_tokens
from run_wall import EXE, EXE_SHA, observe
from plate_check import config

ORIGIN = OLD/"restart_009999.csv"
ARCHIVE = DOC.parent/"v365/research/runs.json"
STEPS = {"full": (50e-6,1), "half1": (25e-6,1), "half2": (25e-6,2)}
SCHEMES = ("adaptive","fixed","newton")


def step_config(scheme, kind):
    """保留几何、物性和二阶空间离散，只为单个后向欧拉物理步配置伪时间收敛。"""
    dt,index=STEPS[kind]
    cfg=config(OLD/"case.cfg")
    for key in ("ITER","WRT_SURFACE_OVERWRITE","WRT_RESTART_OVERWRITE"):
        del cfg[key]
    cfg.update(TIME_DOMAIN="YES",TIME_MARCHING="DUAL_TIME_STEPPING-1ST_ORDER",
               TIME_STEP=str(dt),RESTART_ITER=str(index),TIME_ITER=str(index+1),INNER_ITER="2500",
               MAX_TIME="0.1",CONV_STARTITER="50",CONV_RESIDUAL_MINVAL="-10",
               HISTORY_WRT_FREQ_INNER="1",HISTORY_OUTPUT="(ITER,TIME_DOMAIN,RMS_RES,HEAT)",
               OUTPUT_WRT_FREQ="(1,1)",SCREEN_WRT_FREQ_INNER="50",RESTART_FILENAME="restart",SURFACE_FILENAME="surface")
    if scheme=="fixed":
        cfg.update(CFL_ADAPT="NO",HISTORY_OUTPUT="(ITER,TIME_DOMAIN,RMS_RES,HEAT,CFL_NUMBER,LINSOL)")
    elif scheme=="newton":
        cfg.update(NEWTON_KRYLOV="YES",NEWTON_KRYLOV_IPARAM="(0,3,1)",NEWTON_KRYLOV_DPARAM="(0,0.1,-6,1e-6,1)",
                   CFL_ADAPT="NO",CFL_NUMBER="100",LINEAR_SOLVER_ERROR="0.05",INNER_ITER="250",
                   HISTORY_OUTPUT="(ITER,TIME_DOMAIN,RMS_RES,MAX_RES,HEAT,CFL_NUMBER,LINSOL)",
                   VOLUME_OUTPUT="(COORDINATES,SOLUTION,PRIMITIVE,RESIDUAL)")
    else:
        assert scheme=="adaptive"
    return cfg


def check_source(path,digest):
    """本机输入按原字节核对；已合并的三份执行脚本依其运行时快照核对。"""
    path=path if path.is_absolute() else ROOT/path
    if path.is_file() and sha(path)==digest:
        return
    archive=DOC/"research"
    fingerprints=read(archive/"executed.json")
    assert path.parent==DOC and fingerprints.get(path.name)==digest,(path,digest)
    with zipfile.ZipFile(archive/"executed.zip") as z:
        assert hashlib.sha256(z.read(path.name)).hexdigest()==digest


def inspect(folder, kind, scheme):
    dt,index=STEPS[kind]
    receipt=terminal(folder)
    for path,h in receipt["inputSHA256"].items():
        check_source(Path(path),h)
    cfg=config(folder/"case.cfg")
    expected=step_config(scheme,kind)
    assert cfg==expected
    # 各次原运行的键次序及 Windows 换行也必须能由交付入口重现。
    text="".join(f"{k}= {v}\n" for k,v in expected.items())
    assert (folder/"case.cfg").read_text(encoding="ascii")==text
    hist_path=folder/f"history_{index:05d}.csv"
    h=numbers(hist_path)
    n=len(h["Inner_Iter"])
    np.testing.assert_array_equal(h["Inner_Iter"],np.arange(n))
    assert (h["Time_Iter"]==index).all() and set(k for k in h if k.startswith("rms["))==set(RES)
    np.testing.assert_allclose(h["Time_Step"],dt,rtol=1e-7,atol=0)
    assert len(np.unique(h["Cur_Time"]))==1
    norms=np.column_stack([h[k] for k in RES])
    assert np.isfinite(norms).all()
    frame=folder/f"restart_{index:05d}.csv"
    surf=folder/f"surface_{index:05d}.vtk"
    seed=folder/f"seed_{index-1:05d}.csv"
    a,b=numbers(frame),numbers(seed)
    for key in ("PointID","x","y"):
        np.testing.assert_array_equal(a[key],b[key])
    assert len(a["x"])==9701
    rho=a["Density"]
    p=.4*(a["Energy"]-.5*(a["Momentum_x"]**2+a["Momentum_y"]**2)/rho)
    np.testing.assert_allclose(p,a["Pressure"],rtol=1e-10)
    np.testing.assert_allclose(p/rho,a["Temperature"],rtol=1e-10)
    assert np.isfinite(np.stack(list(a.values()))).all() and min(rho)>0 and min(p)>0
    assert not np.array_equal(a["Energy"],b["Energy"])
    assert max(abs(a["Momentum_x"][:109]))==max(abs(a["Momentum_y"][:109]))==0
    xyz,fields=vtk_tokens(surf)
    assert len(xyz)==109 and "Heat_Flux" in fields and "Skin_Friction_Coefficient" in fields
    last=norms[-20:].max(axis=0)
    return dict(branch=kind,step=index,dtSeconds=dt,elapsedFromSeedSeconds=dt*index,
                nativeCurTimeSeconds=float(h["Cur_Time"][-1]),innerIterations=n,
                last20Max=dict(zip(RES,last.tolist())),innerPassed=bool(n>=20 and max(last)<-8),
                seedSHA256=sha(seed),frameSHA256=sha(frame),wallSHA256=sha(surf),historySHA256=sha(hist_path),
                cfgSHA256=sha(folder/"case.cfg"),meshSHA256=sha(folder/"mesh.su2"),
                seconds=receipt["elapsedSeconds"],peakRSSBytes=receipt["peakFamilyRSSBytes"])


def execute(out,kind,scheme):
    out=out.resolve()
    out.relative_to(BASE)
    dt,index=STEPS[kind]
    old=read(ARCHIVE)["warm2"];terminal(OLD)
    assert sha(EXE)==EXE_SHA and sha(ORIGIN)==old["windows"][-1]["sourceRestartSHA256"]
    assert sha(OLD/"case.cfg")==old["configSHA256"] and sha(OLD/"mesh.su2")==old["meshSHA256"]
    assert sha(OLD/"condition.json")=="712d77252229a50d27592266654e0e819b18b3e57d34c23adcf2ad00a0384cf1"
    source=ORIGIN
    if kind=="half2":
        previous=inspect(out/"half1","half1",scheme)
        assert previous["innerPassed"] and previous["seedSHA256"]==sha(ORIGIN)
        source=out/"half1/restart_00001.csv"
    folder=out/kind
    folder.mkdir(parents=True,exist_ok=False)
    cfg=folder/"case.cfg"
    cfg.write_text("".join(f"{k}= {v}\n" for k,v in step_config(scheme,kind).items()),encoding="ascii")
    mesh=folder/"mesh.su2";seed=folder/f"seed_{index-1:05d}.csv"
    shutil.copyfile(OLD/"mesh.su2",mesh);shutil.copyfile(source,seed)
    assert sha(seed)==sha(source)
    locks=[EXE,mesh,cfg,seed,source,OLD/"condition.json",OLD/"case.cfg",ARCHIVE,DOC/"PLAN.md",
           Path(__file__),ROOT/"scripts/runSU2NativeDebugger.py"]
    if scheme!="adaptive":
        locks.append(DOC/("RECOVERY.md" if scheme=="fixed" else "NONLINEAR.md"))
    result=observe([str(EXE),cfg.name,"--threads","2"],folder,
                   {str(p):sha(p) for p in locks},150,64*1024**2)
    if result["debuggerExitCode"] or result["stopReason"]:
        raise RuntimeError(f"原生运行中止，见 {folder}")
    record=inspect(folder,kind,scheme)
    save(folder/"step.json",record)
    print(kind,record["innerIterations"],record["last20Max"],record["innerPassed"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path);parser.add_argument("kind",choices=STEPS)
    parser.add_argument("--scheme",choices=SCHEMES,default="adaptive")
    args=parser.parse_args();execute(args.out,args.kind,args.scheme)
