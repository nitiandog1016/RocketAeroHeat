"""从同一已绑定初场核实连续 20 行物理步验收，不跨物理时刻续接旧试探态。"""
import argparse
import os
from pathlib import Path
import sys
import numpy as np
from check_fields import columns,rhs_at_endpoint
from run_relaxation import ENGINE, ENGINE_SHA, FULL, HERE, PLAN, RUN, frozen, observe, read, save, sha, terminal
from matched_step import RMS, audit


def run(source,out):
    terminal(source)
    prior=read(RUN/"finish01/first_step.json")
    assert not prior["acceptedFirstStep"] and prior["innerRows"]==180
    hist=frozen(source/"history_00001.csv","20b9a134b492fc2d81be1b179439d76358e73e0838d4ee0e21ce87aca32a4b2c")
    assert np.all(columns(hist,RMS)[-12:]<-8)
    frozen(source/"seed_00000.csv","ea94eb8a3b668a47e1d3e66e59065e704ab05c57ae752fcd6ab75fa112f42445")
    frozen(ENGINE,ENGINE_SHA)
    cfg=(source/"case.cfg").read_text("ascii")
    assert cfg.count("INNER_ITER= 180\n")==1
    out.mkdir(parents=True,exist_ok=False)
    files=[PLAN,HERE/"WINDOW_CHECK.md",Path(__file__).resolve(),ROOT_SUPERVISOR,ENGINE,hist,source/"case.cfg",
           source/"debugger_process.json",RUN/"finish01/first_step.json"]
    for key in ("seed_00000.csv","body.json","geometry.json","mesh.su2","nodes.npz"):
        os.link(source/key,out/key)
        files.append(out/key)
    (out/"case.cfg").write_text(cfg.replace("INNER_ITER= 180\n","INNER_ITER= 200\n"),encoding="ascii")
    files.append(out/"case.cfg")
    rec=observe([str(ENGINE),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in files},360,180*1024**2)
    return rec["debuggerExitCode"]==0 and rec["stopReason"] is None


def verify(source,out):
    terminal(source)
    old=RUN/"step01"
    assert (source/"case.cfg").read_text("ascii")==(old/"case.cfg").read_text("ascii").replace("INNER_ITER= 180\n","INNER_ITER= 200\n")
    for name in ("mesh.su2","body.json","seed_00000.csv","nodes.npz"):
        assert sha(source/name)==sha(old/name)
    hist=source/"history_00001.csv"
    original=frozen(old/hist.name,"20b9a134b492fc2d81be1b179439d76358e73e0838d4ee0e21ce87aca32a4b2c")
    prefix=hist.read_bytes().startswith(original.read_bytes())
    assert prefix,"旧物理步的 180 行前缀不一致"
    state=audit(source)
    with np.load(source/"nodes.npz") as z:
        xyz=z["xyz"]
    h=columns(hist,["Inner_Iter"]+RMS)
    rhs=rhs_at_endpoint(source/"restart_00001.csv",xyz,h[-1,1:])
    good=np.all(h[:,1:]<-8,axis=1)
    windows=np.flatnonzero(np.convolve(good.astype(int),np.ones(20,int),mode="valid")==20)
    end=int(h[windows[0]+19,0]) if len(windows) else None
    state.update(acceptedFirstStep=bool(state["acceptedFirstStep"] and rhs["maxLogDifferenceFromHistory"]<1e-7),
                 firstCompleteWindowEndIndex=end,referencePrefixByteIdentical=prefix,referenceRows=180,nodeRHS=rhs)
    save(out/"window.json",state)


ROOT_SUPERVISOR=FULL.parents[2]/"scripts/runSU2NativeDebugger.py"
if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("run","audit","_audit"))
    cli.add_argument("source",type=Path)
    cli.add_argument("output",type=Path)
    args=cli.parse_args()
    folder,out=args.source.resolve(),args.output.resolve()
    folder.relative_to(RUN);out.relative_to(RUN)
    if args.stage=="run":
        raise SystemExit(0 if run(folder,out) else 1)
    elif args.stage=="_audit":
        verify(folder,out)
    else:
        out.mkdir(parents=True,exist_ok=False)
        files=[PLAN,HERE/"WINDOW_CHECK.md",Path(__file__).resolve(),HERE/"check_fields.py",
               folder/"restart_00001.csv",folder/"history_00001.csv",folder/"debugger_process.json"]
        for key in ("MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
            os.environ[key]="1"
        rec=observe([sys.executable,str(Path(__file__).resolve()),"_audit",str(folder),str(out)],out,
                    {str(p):sha(p) for p in files},90,128*1024**2)
        raise SystemExit(0 if rec["debuggerExitCode"]==0 and rec["stopReason"] is None else 1)
