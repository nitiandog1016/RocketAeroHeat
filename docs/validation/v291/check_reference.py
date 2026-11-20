"""从 v290 固定体场核对可指定步号的热量账和面元映射是否逐字节复现原发布包。"""
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from run_half import HERE, MANIFEST, ROOT, RUN, SOURCE, input_file, observe, sha
from thermistorStartupSurface import write


if __name__=="__main__":
    out=Path(sys.argv[1]).resolve()
    out.relative_to(RUN)
    out.mkdir(parents=True,exist_ok=False)
    script=ROOT/"scripts/thermistorStartupSurface.py"
    old=ROOT/"outputs/v290_validation"
    reference=ROOT/"cases/thermistor_b/startup_6ms.zip"
    files=[HERE/"PLAN.md",Path(__file__).resolve(),script,MANIFEST,reference,
           ROOT/"scripts/beadSurfaceResult.py",input_file("mesh.su2"),input_file("seed_00000.csv"),
           input_file("restart_00001.csv"),input_file("body.json"),input_file("nodes.npz")]
    files.extend(old/name/"debugger_process.json" for name in ("observe01","release01"))
    command=[sys.executable,str(script),str(SOURCE),str(old/"observe01"),str(old/"release01"),str(out/"result")]
    rec=observe(command,out,{str(p):sha(p) for p in files},120,100*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
    rebuilt=out/"result/startup_6ms.zip"
    assert rebuilt.read_bytes()==reference.read_bytes()
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<100*1024**2
    write(out/"reference.json",dict(byteIdentical=True,archiveSHA256=sha(rebuilt),bytes=rebuilt.stat().st_size,
                                   elapsedSeconds=rec["elapsedSeconds"],nativeSolverRerun=False))
