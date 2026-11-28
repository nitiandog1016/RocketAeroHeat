"""仅为第 249 次组装取得官方更新前的状态，然后检查一次原生残差。"""
import argparse
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from build_observer import BASE,BIN,DOC,ROOT,sha
from replay import PARENT
sys.path.insert(0,str(DOC.parent/"v366"))
from physical_step import EXE,EXE_SHA,RES,inspect,numbers,observe,read,save,terminal


def run(out,mode,build=None,prefix=None):
    out=out.resolve();out.relative_to(BASE)
    prior=inspect(PARENT,"full","newton")
    assert prior["innerIterations"]==250 and sha(EXE)==EXE_SHA
    out.mkdir(parents=True,exist_ok=False)
    locks={}
    for name in ("mesh.su2","seed_00000.csv"):
        shutil.copyfile(PARENT/name,out/name)
        locks[str(out/name)]=sha(PARENT/name)
    cfg=out/"case.cfg"
    text=(PARENT/"case.cfg").read_text(encoding="ascii")
    assert text.count("INNER_ITER= 250\n")==1
    cfg.write_text(text.replace("INNER_ITER= 250\n",f"INNER_ITER= {249 if mode=='prefix' else 1}\n"),encoding="ascii")
    paths=[cfg,Path(__file__).resolve(),DOC/"STATE_MATCH.md",DOC/"PLAN.md",PARENT/"case.cfg",
           PARENT/"history_00001.csv",ROOT/"scripts/runSU2NativeDebugger.py"]
    exe=EXE
    if mode=="sample":
        build=build.resolve();prefix=prefix.resolve()
        exe=build/"tail_observe.exe"
        assert sha(exe)==read(build/"build.json")["executableSHA256"]
        marker=read(prefix/"prefix.json")
        checkpoint=prefix/"restart_00001.csv"
        assert sha(checkpoint)==marker["frameSHA256"] and marker["maxPrefixLogDifference"]<1e-9
        state=numbers(checkpoint);n=len(state["PointID"])
        np.testing.assert_array_equal(state["PointID"],np.arange(n))
        assert n==9701
        # 文件只含固定网格的四个无量纲守恒变量，点序已与原生原面对应关系核对。
        buf=np.column_stack([state[k] for k in ("Density","Momentum_x","Momentum_y","Energy")])
        buf.astype("<f8").tofile(out/"assembly_state.bin")
        paths.extend([checkpoint,out/"assembly_state.bin",build/"build.json",DOC/"match.hpp"])
    paths.append(exe)
    locks.update({str(p):sha(p) for p in paths})
    old=os.environ.copy()
    try:
        for key in ("RAH_M030_TAIL_BUDGET","RAH_M030_ASSEMBLY_STATE"):
            os.environ.pop(key,None)
        if mode=="sample":
            os.environ.update(PATH=str(BIN)+os.pathsep+old["PATH"],RAH_M030_TAIL_BUDGET="1",RAH_M030_ASSEMBLY_STATE="1")
        observe([str(exe),"case.cfg","--threads","2"],out,locks,150,64*1024**2)
    finally:
        for key in ("PATH","RAH_M030_TAIL_BUDGET","RAH_M030_ASSEMBLY_STATE"):
            if key in old:os.environ[key]=old[key]
            else:os.environ.pop(key,None)
    rc=terminal(out)
    h=numbers(out/"history_00001.csv");reference=numbers(PARENT/"history_00001.csv")
    assert (h["Time_Step"]==50e-6).all() and (h["Time_Iter"]==1).all()
    count=249 if mode=="prefix" else 1
    np.testing.assert_array_equal(h["Inner_Iter"],np.arange(count))
    idx=slice(0,249) if mode=="prefix" else slice(249,250)
    errors={k:float(max(abs(h[k]-reference[k][idx]))) for k in RES}
    r=dict(mode=mode,innerIterations=count,logDifferences=errors,maxPrefixLogDifference=max(errors.values()) if mode=="prefix" else None,
           referenceIteration=248 if mode=="prefix" else 249,frameSHA256=sha(out/"restart_00001.csv"),
           seedSHA256=sha(out/"seed_00000.csv"),seconds=rc["elapsedSeconds"],peakRSSBytes=rc["peakFamilyRSSBytes"])
    save(out/f"{mode}.json",r);print(r,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("out",type=Path);p.add_argument("mode",choices=("prefix","sample"))
    p.add_argument("--build",type=Path);p.add_argument("--prefix",type=Path)
    a=p.parse_args();run(a.out,a.mode,a.build,a.prefix)
