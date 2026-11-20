"""监督民用珠体嵌套网格、物理初场传递和几何回读，绑定每个真实输入。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorNearWall import prepare
from thermistorBead import read,sha
from runSU2NativeDebugger import observe

SOURCE=ROOT/"outputs/v289_validation/window01"
RUN=ROOT/"outputs/v293_validation"
NATIVE=ROOT/"outputs/v292_validation/endpoint_observe/before_wall.csv"
COMPARE=HERE.parent/"v292/evidence/comparison.json"


def bound_source(name):
    path=SOURCE/name
    manifest=read(HERE.parent/"v290/evidence/local_assets.json")
    assert sha(path)==manifest[path.relative_to(ROOT).as_posix()]["sha256"]
    return path


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output",type=Path);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();out=a.output.resolve();out.relative_to(RUN)
    if a.worker:
        prepare(SOURCE,NATIVE,COMPARE,out)
    else:
        out.mkdir(parents=True,exist_ok=False)
        files=[HERE/"PLAN.md",Path(__file__).resolve(),ROOT/"scripts/thermistorNearWall.py",
               ROOT/"scripts/cm011TunnelSector.py",ROOT/"scripts/thermistorStartupSurface.py",
               ROOT/"scripts/runSU2NativeDebugger.py",COMPARE,NATIVE]
        assert sha(NATIVE)==read(HERE.parent/"v292/evidence/first_energy.json")["energy"]["observedWallSHA256"]
        files.extend(bound_source(n) for n in ("mesh.su2","seed_00000.csv","nodes.npz","geometry.json","body.json"))
        rec=observe([sys.executable,str(Path(__file__).resolve()),str(out),"--worker"],out,
                    {str(p):sha(p) for p in files},180,200*1024**2)
        assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
