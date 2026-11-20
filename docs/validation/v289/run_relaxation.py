"""B 型珠体 6 m/s 全周初场的扰动对照、定常松弛及有条件的原生首步。"""
import argparse
import gzip
import hashlib
import os
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/"v288"
sys.path.insert(0,str(OLD))
from matched_step import ENGINE, ENGINE_SHA, FULL, FULL_CONFIG_SHA, ROOT, frame, observe, read, save, sha
from compare_residuals import one_step

RUN=ROOT/"outputs/v289_validation"
LAST=ROOT/"outputs/v288_validation"
ASSETS=read(OLD/"evidence/local_assets.json")["files"]
FULL_ASSETS=read(OLD.parent/"v287/evidence/local_assets.json")["files"]
PLAN=HERE/"PLAN.md"
PROBE=LAST/"probe_mirror01"
EXECUTED_SHA="14bc79e8b14d359a5d6f6bab06aba27056c9e0bac3bb86ef36a1dcc3e47e3eff"


def terminal(folder):
    """逐项验证终态输入；历史运行精确绑定其归档源码，保留当时的回执。"""
    rec=read(folder/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and rec["stopReason"] is None
    assert rec["supervisorError"] is None and not rec["survivingProcessIDs"] and all(rec["inputsStable"].values())
    for name,digest in rec["inputSHA256"].items():
        if Path(name)==Path(__file__).resolve() and digest==EXECUTED_SHA:
            original=gzip.decompress((HERE/"evidence/executed_driver.py.gz").read_bytes())
            assert hashlib.sha256(original).hexdigest()==digest
        else:
            frozen(Path(name),digest)
    return rec


def frozen(path,digest):
    if sha(path)!=digest:
        raise ValueError(f"输入散列不符 {path}")
    return path


def mode(values,ids):
    rings=values[ids[:-1]]
    return float(2*np.abs(np.fft.rfft(rings,axis=0)[1]).max()/len(rings))


def write_seed(path,a):
    header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"'
    with path.open("x",encoding="ascii") as f:
        np.savetxt(f,a,fmt=["%d"]+["%.17g"]*9,header=header,comments="",delimiter=",")


def calm(out):
    """直接取扇区对应的轴向速度，独立核对原扰动公式以及逐节点还原。"""
    seed=frozen(LAST/"mirror01/seed_00000.csv",ASSETS["mirror01/seed_00000.csv"]["sha256"])
    base=frozen(LAST/"base01/restart_001505.csv",ASSETS["base01/restart_001505.csv"]["sha256"])
    with np.load(frozen(FULL/"nodes.npz",FULL_ASSETS["first02/nodes.npz"]["sha256"])) as grid:
        xyz,ids=grid["xyz"],grid["ids"]
    with np.load(frozen(LAST/"base01/nodes.npz","8d3911634283e2df8d4fe61d3ee50a89cba4503b2eb16721cee365c97de01629")) as grid:
        parent=frame(base,grid["xyz"]);sector=grid["ids"]
    data=frame(seed,xyz)
    nodes,first=np.unique(ids[:-1].ravel(),return_index=True)
    np.testing.assert_array_equal(nodes,np.arange(len(xyz)))
    ray,station=np.divmod(first,ids.shape[1])
    folded=ray%8
    ux=parent[sector[np.minimum(folded,8-folded),station],5]
    delta=data[:,5]-ux
    x,r=xyz[ids[0,station],:2].T
    ellipse=((x+1.7)/.8)**2+((r-.55)/.45)**2
    expected=1e-3*np.maximum(1-ellipse,0.)**2*np.cos(ray*np.pi/24)
    np.testing.assert_allclose(delta,expected,atol=2e-16,rtol=0)
    np.testing.assert_array_equal(ux+delta,data[:,5])
    assert mode(ux,ids)<1e-14 and np.all(delta[ids[:,:65]]==0)
    data[:,5]=ux
    write_seed(out/"seed_00000.csv",data)
    np.save(out/"impulse.npy",delta)
    save(out/"calm.json",dict(seedSHA256=sha(out/"seed_00000.csv"),sourceSeedSHA256=sha(seed),
         baseSHA256=sha(base),impulseSHA256=sha(out/"impulse.npy"),calmM1=mode(ux,ids),
         impulseM1=mode(delta,ids),changedNodes=int(np.count_nonzero(delta))))


def control(folder):
    """用原生节点和历史确认控制量，并取下游、轴线两个条件的交集。"""
    rec=terminal(folder)
    if rec["inputSHA256"].get(str(Path(__file__).resolve()))==EXECUTED_SHA:
        # 已交付的原始观察按输入与输出指纹回读。新运行按原生场重算。
        result=read(HERE/"evidence/control.json")
        frozen(folder/"seed_00000.csv",result["seedSHA256"])
        frozen(folder/"restart_00001.csv",result["frameSHA256"])
        frozen(folder/"history_00001.csv",read(HERE/"evidence/local_assets.json")["files"]["probe01/history_00001.csv"]["sha256"])
    else:
        result=one_step(folder)
    with np.load(folder/"nodes.npz") as grid:
        xyz=grid["xyz"]
    import csv
    with (folder/"restart_00001.csv").open(encoding="ascii") as f:
        keys=next(csv.reader(f))
        h=np.loadtxt(f,delimiter=",",usecols=[keys.index("Residual_Enthalpy")])
    mask=(np.hypot(xyz[:,1],xyz[:,2])==0)&(xyz[:,0]<-.5*1.179248562795375)
    share=float(np.sum(h[mask]**2)/np.sum(h*h))
    result["downstreamAxisSquaredShare"]=share
    native=result["nativeLogRMS"]
    result["needsRelaxation"]=bool(result["pressureEnthalpyLogMatch"] and result["initialMode"]<1e-14 and
                                  max(native["rms[P]"],native["rms[h]"])>=-8 and share>.10)
    return result


def steady_config():
    source=frozen(LAST/"base01/case.cfg",ASSETS["base01/case.cfg"]["sha256"])
    cfg=source.read_text("ascii")
    for old,new in [("MARKER_SYM= (sym0, sym1)\n",""),("ITER= 7000\n","ITER= 2000\n"),
                    ("OUTPUT_WRT_FREQ= 800","OUTPUT_WRT_FREQ= 400"),
                    ("VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)",
                     "VOLUME_OUTPUT= (COORDINATES, SOLUTION, RESIDUAL)\nWRT_RESTART_COMPACT= NO")]:
        assert cfg.count(old)==1,old
        cfg=cfg.replace(old,new)
    return cfg


def source_geometry(calm_meta,seed):
    """把无扰动种子的实际父场与体网自带的旧来源分开写入工况记录。"""
    geom=read(FULL/"geometry.json")
    geom.update(initialGuessFrameSHA256=calm_meta["baseSHA256"],initialSeedSHA256=sha(seed),
                peakAxialPerturbationOverU=0.,initialGuessTransfer="unperturbed five-ray mirror")
    assert geom["initialSeedSHA256"]==calm_meta["seedSHA256"] and calm_meta["calmM1"]<1e-14
    return geom


def run(stage,source,out):
    """所有数值运行在独立目录中进行。共享体网及种子只作为只读输入。"""
    files=[PLAN,Path(__file__).resolve(),OLD/"matched_step.py",OLD/"compare_residuals.py",
           OLD.parent/"v287/audit_step.py",ROOT/"scripts/runSU2NativeDebugger.py"]
    out.mkdir(parents=True,exist_ok=False)
    if stage in ("seed","control"):
        command=[sys.executable,str(Path(__file__).resolve()),"_"+stage,str(source),str(out)]
        if stage=="seed":
            files.extend([LAST/"mirror01/seed_00000.csv",LAST/"base01/restart_001505.csv",
                          LAST/"base01/nodes.npz",FULL/"nodes.npz"])
        else:
            files.extend(source/k for k in ("restart_00001.csv","history_00001.csv","debugger_process.json"))
        limits=(90,128)
    else:
        frozen(ENGINE,ENGINE_SHA)
        terminal(source)
        for name in ("mesh.su2","body.json","nodes.npz"):
            os.link(frozen(FULL/name,FULL_ASSETS["first02/"+name]["sha256"]),out/name)
        geom=read(FULL/"geometry.json")
        if stage=="steady":
            check=read(source/"control.json")
            assert check["needsRelaxation"]
            initial=frozen(RUN/"calm01/seed_00000.csv",check["seedSHA256"])
            geom=source_geometry(read(RUN/"calm01/calm.json"),initial)
            input_name="initial_flow.csv";cfg=steady_config();limits=(1200,768)
            files.extend([source/"control.json",LAST/"base01/case.cfg",RUN/"calm01/calm.json"])
            threads=2
        elif stage=="probe":
            info=read(source/"calm.json")
            initial=frozen(source/"seed_00000.csv",info["seedSHA256"])
            cfg=frozen(PROBE/"case.cfg",ASSETS["probe_mirror01/case.cfg"]["sha256"]).read_text("ascii")
            input_name="seed_00000.csv";limits=(60,180);threads=1
            geom=source_geometry(info,initial)
            files.extend([source/"calm.json",PROBE/"case.cfg"])
        else:
            info=read(source/"relaxed_seed.json")
            assert info["acceptedForInitialization"]
            initial=frozen(source/"seed_00000.csv",info["seedSHA256"])
            cfg=frozen(FULL/"case.cfg",FULL_CONFIG_SHA).read_text("ascii")
            before="VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)"
            assert cfg.count(before)==1
            cfg=cfg.replace(before,before[:-1]+", RESIDUAL)")
            geom["initialGuessFrameSHA256"]=info["relaxationFrameSHA256"]
            geom["initialSeedSHA256"]=info["seedSHA256"]
            geom["baseAxialM1"]=info["baseM1"];geom["addedAxialM1"]=info["addedM1"]
            input_name="seed_00000.csv";limits=(360,180);threads=1
            files.extend([source/"relaxed_seed.json",FULL/"case.cfg"])
        os.link(initial,out/input_name)
        (out/"case.cfg").write_text(cfg,encoding="ascii")
        save(out/"geometry.json",geom)
        command=[str(ENGINE),"case.cfg","--threads",str(threads)]
        files.extend([ENGINE,out/input_name]+[out/k for k in ("mesh.su2","nodes.npz","body.json","case.cfg","geometry.json")])
    for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[key]="1"
    rec=observe(command,out,{str(p):sha(p) for p in files},limits[0],limits[1]*1024**2)
    return rec["debuggerExitCode"]==0 and rec["stopReason"] is None


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("seed","probe","control","steady","step","_seed","_control"))
    cli.add_argument("source",type=Path)
    cli.add_argument("output",type=Path)
    args=cli.parse_args();src=args.source.resolve();dest=args.output.resolve()
    dest.relative_to(RUN)
    if args.stage=="_seed":
        calm(dest)
    elif args.stage=="_control":
        save(dest/"control.json",control(src))
    else:
        raise SystemExit(0 if run(args.stage,src,dest) else 1)
