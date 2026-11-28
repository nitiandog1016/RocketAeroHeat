"""用官方步末的当前状态和原物理上一步，回放一次相同方程的原生组装。"""
import argparse
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from build_ledger import BASE,BIN,DOC,ROOT,sha
sys.path.insert(0,str(DOC.parent/"v369"))
from short_step import GRID_HASH,SEED_HASH,RES,numbers,observe,read,save,terminal

PARENT=ROOT/"outputs/v369_validation/pair01/half1"
FRAME_HASH="7c783fbdee774c07729c86f6cc39daa778f84955a3111ef5de0fd10bea88f970"
KEYS=("Density","Momentum_x","Momentum_y","Energy")


def conserved(path):
    a=numbers(path)
    np.testing.assert_array_equal(a["PointID"],np.arange(9701))
    return np.column_stack([a[k] for k in KEYS])


def execute(build,out):
    build=build.resolve();out=out.resolve()
    out.relative_to(BASE);build.relative_to(BASE)
    exe=build/"ledger.exe"
    assert sha(exe)==read(build/"build.json")["executableSHA256"]
    assert sha(DOC/"ledger.hpp")==read(build/"build.json")["headerSHA256"]
    terminal(PARENT)
    frame=PARENT/"restart_00001.csv";prior=PARENT/"seed_00000.csv";mesh=PARENT/"mesh.su2";cfg=PARENT/"case.cfg"
    record=read(DOC.parent/"v369/evidence/runs.json")["half1"]
    assert record["passed"] and record["iterations"]==740 and sha(cfg)==record["cfgSHA256"]
    assert sha(prior)==SEED_HASH and sha(frame)==FRAME_HASH and sha(mesh)==GRID_HASH
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(frame,out/"seed_00000.csv")
    shutil.copyfile(mesh,out/mesh.name)
    conserved(prior).astype("<f8").tofile(out/"previous.bin")
    text=cfg.read_text(encoding="ascii")
    assert text.count("INNER_ITER= 2500\n")==1
    (out/"case.cfg").write_text(text.replace("INNER_ITER= 2500\n","INNER_ITER= 1\n"),encoding="ascii")
    inputs=[exe,build/"build.json",DOC/"PLAN.md",DOC/"ledger.hpp",Path(__file__).resolve(),
            frame,prior,mesh,cfg,*(out/n for n in ("case.cfg","seed_00000.csv","mesh.su2","previous.bin"))]
    os.environ.update(PATH=str(BIN)+os.pathsep+os.environ["PATH"],RAH_M030_STEP_LEDGER="1")
    observe([str(exe),"case.cfg","--threads","2"],out,{str(p):sha(p) for p in inputs},120,80*1024**2)
    proc=terminal(out)
    raw=np.fromfile(out/"state.bin",dtype="<f8").reshape(9701,11)
    np.testing.assert_array_equal(raw[:,3:7],conserved(frame))
    np.testing.assert_array_equal(raw[:,7:],conserved(prior))
    h=numbers(out/"history_00001.csv")
    assert h["Inner_Iter"].tolist()==[0] and h["Time_Iter"].tolist()==[1] and h["Time_Step"].tolist()==[6.25e-6]
    rhs=np.fromfile(out/"rhs.bin",dtype="<f8").reshape(9701,4)
    log_rms=np.log10(np.sqrt(np.mean(rhs**2,axis=0)))
    np.testing.assert_allclose(log_rms,[h[k][0] for k in RES],atol=1e-7,rtol=0)
    save(out/"replay.json",dict(frameSHA256=FRAME_HASH,previousSHA256=SEED_HASH,meshSHA256=GRID_HASH,
         historySHA256=sha(out/"history_00001.csv"),postUpdateLogRMS=dict(zip(RES,log_rms.tolist())),
         seconds=proc["elapsedSeconds"],peakRSSBytes=proc["peakFamilyRSSBytes"],
         stageFiles={p.name:sha(p) for p in sorted(out.glob("*.bin"))}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("build",type=Path);p.add_argument("out",type=Path)
    args=p.parse_args()
    execute(args.build,args.out)
