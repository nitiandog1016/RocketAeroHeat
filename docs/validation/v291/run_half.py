"""在民用热敏珠的固定全周网格上组织一个半步，并绑定实际旧态和原生终场。"""
import argparse
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT/"scripts"), str(HERE.parent/"v290")]
from observe_heat import BIN, ENGINE_SHA, OFFICIAL, SOURCE, replace
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe
from thermistorStartupSurface import read

RUN = ROOT/"outputs/v291_validation"
MANIFEST = HERE.parent/"v290/evidence/local_assets.json"
MAX_INNER = 240


def input_file(name):
    path = SOURCE/name
    assert sha(path) == read(MANIFEST)[path.relative_to(ROOT).as_posix()]["sha256"], path
    return path


def config(index):
    text = input_file("case.cfg").read_text("ascii")
    text = replace(text, "TIME_STEP= 0.05\n", "TIME_STEP= 0.025\n")
    text = replace(text, "INNER_ITER= 200\n", f"INNER_ITER= {MAX_INNER}\n")
    text = replace(text, "RESTART_ITER= 1\n", f"RESTART_ITER= {index}\n")
    return replace(text, "TIME_ITER= 2\n", f"TIME_ITER= {index+1}\n")


def run(index, out, previous):
    # 每次新建独立目录。第二个半步只从已经核对过的中间场和储热账出发。
    assert index in (1, 2) and sha(OFFICIAL) == ENGINE_SHA
    inputs = [HERE/"PLAN.md", Path(__file__).resolve(), MANIFEST,
              ROOT/"scripts/runSU2NativeDebugger.py", input_file("case.cfg"), OFFICIAL]
    if index == 1:
        assert previous is None
        old = input_file("seed_00000.csv")
    else:
        assert previous is not None
        proof = read(previous)
        assert proof["accepted"] and proof["step"] == 1 and proof["energy"]["relativeImbalance"] < .005
        old = Path(proof["acceptedFrame"])
        assert sha(old) == proof["energy"]["acceptedFrameSHA256"]
        inputs.append(previous)
    out.relative_to(RUN)
    out.mkdir(parents=True, exist_ok=False)
    os.link(old, out/f"seed_{index-1:05}.csv")
    os.link(input_file("mesh.su2"), out/"mesh.su2")
    (out/"case.cfg").write_text(config(index),encoding="ascii")
    inputs.extend([out/f"seed_{index-1:05}.csv", out/"mesh.su2", out/"case.cfg"])
    os.environ["PATH"] = str(BIN)+os.pathsep+os.environ["PATH"]
    receipt = observe([str(OFFICIAL),"case.cfg","--threads","1"],out,
                      {str(p):sha(p) for p in inputs},420,180*1024**2)
    assert not receipt["stopReason"] and receipt["debuggerExitCode"] == 0, out/"debugger.log"
    assert "Exit Success (SU2_CFD)" in (out/"debugger.log").read_text(errors="replace")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("step", choices=[1,2], type=int)
    cli.add_argument("output", type=Path)
    cli.add_argument("--previous", type=Path)
    args = cli.parse_args()
    run(args.step, args.output.resolve(), args.previous.resolve() if args.previous else None)
