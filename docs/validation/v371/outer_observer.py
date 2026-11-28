"""在已收敛的 91×109 细网原生终态上取得与原网同定义的质量和能量各项。"""
import argparse
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from outer_step import BASE,DOC,SHAPE,CFG_HASH,RES,observe,numbers,read,save,sha,terminal,CONS
sys.path.insert(0,str(DOC.parent/"v370"))
from build_ledger import BIN,compile_observer


def header():
    source=DOC.parent/"v370/ledger.hpp"
    assert sha(source)==read(DOC.parent/"v370/evidence/checks.json")["successfulBuild"]["headerSHA256"]
    original=source.read_bytes()
    assert original.count(b"n!=9701")==1
    return original.replace(b"n!=9701",f"n!={np.prod(SHAPE)}".encode("ascii"))


def sample(build,step,out):
    record=read(step/"step.json")
    terminal(step)
    assert record["passed"] and record["elapsedUs"]==6.25 and (build/"ledger.hpp").read_bytes()==header()
    exe=build/"ledger.exe";manifest=read(build/"build.json")
    assert sha(exe)==manifest["executableSHA256"] and sha(build/"ledger.hpp")==manifest["headerSHA256"]
    files={"frame":step/"restart_00001.csv","mesh":step/"mesh.su2","seed":step/"seed_00000.csv"}
    assert all(sha(p)==record[k+"SHA256"] for k,p in files.items())
    assert sha(step/"case.cfg")==CFG_HASH
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(files["frame"],out/"seed_00000.csv")
    shutil.copyfile(files["mesh"],out/"mesh.su2")
    states=[]
    for name in ("frame","seed"):
        native=numbers(files[name])
        np.testing.assert_array_equal(native["PointID"],np.arange(np.prod(SHAPE)))
        states.append(np.column_stack([native[k] for k in CONS]))
    states[1].astype("<f8").tofile(out/"previous.bin")
    text=(step/"case.cfg").read_text(encoding="ascii")
    assert text.count("INNER_ITER= 2500\n")==1
    (out/"case.cfg").write_text(text.replace("INNER_ITER= 2500\n","INNER_ITER= 1\n"),encoding="ascii")
    locked=[exe,build/"build.json",build/"ledger.hpp",DOC/"PLAN.md",Path(__file__).resolve(),step/"step.json",step/"case.cfg",
            *files.values(),*(out/k for k in ("case.cfg","seed_00000.csv","mesh.su2","previous.bin"))]
    os.environ.update(PATH=str(BIN)+os.pathsep+os.environ["PATH"],RAH_M030_STEP_LEDGER="1")
    observe([str(exe),"case.cfg","--threads","2"],out,{str(p):sha(p) for p in locked},120,80*1024**2)
    proc=terminal(out)
    raw=np.fromfile(out/"state.bin",dtype="<f8").reshape(-1,11)
    for offset,state in zip((3,7),states):
        np.testing.assert_array_equal(raw[:,offset:offset+4],state)
    h=numbers(out/"history_00001.csv")
    assert h["Inner_Iter"].tolist()==[0] and h["Time_Iter"].tolist()==[1] and h["Time_Step"].tolist()==[6.25e-6]
    rhs=np.fromfile(out/"rhs.bin",dtype="<f8").reshape(np.prod(SHAPE),4)
    log_rms=np.log10(np.sqrt(np.mean(rhs**2,axis=0)))
    np.testing.assert_allclose(log_rms,[h[k][0] for k in RES],atol=1e-7,rtol=0)
    save(out/"replay.json",dict(frameSHA256=record["frameSHA256"],previousSHA256=record["seedSHA256"],meshSHA256=record["meshSHA256"],
         historySHA256=sha(out/"history_00001.csv"),postUpdateLogRMS=dict(zip(RES,log_rms.tolist())),seconds=proc["elapsedSeconds"],
         peakRSSBytes=proc["peakFamilyRSSBytes"],stageFiles={p.name:sha(p) for p in sorted(out.glob("*.bin"))}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("build",type=Path)
    p.add_argument("--sample",nargs=2,type=Path,metavar=("STEP","OUTPUT"))
    a=p.parse_args();build=a.build.resolve();build.relative_to(BASE)
    if a.sample:
        step,out=[v.resolve() for v in a.sample];step.relative_to(BASE);out.relative_to(BASE)
        sample(build,step,out)
    else:
        compile_observer(build,header())
