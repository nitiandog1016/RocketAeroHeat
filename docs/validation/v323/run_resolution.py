"""执行登记的银球只读诊断，绑定原始输入并保存资源回执。"""
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from buildThermistorWireDriver import invoke
from shortSphere65 import sha, write_json


if __name__ == "__main__":
    out = ROOT/"outputs/v323_validation"/sys.argv[1]
    assert sys.argv[1].isidentifier()
    out.mkdir(parents=True, exist_ok=False)
    work = ROOT/"outputs/v322_validation/run01"
    script = ROOT/"scripts/shortSphereWakeResolution.py"
    paths = [script, HERE/"PLAN.md", Path(__file__).resolve(), ROOT/"cases/short_sphere/test65.json"]
    for branch, n in [("coarse",32), ("fine",64)]:
        d = work/branch/f"step{n:03d}"
        paths.extend(d/name for name in ("mesh.su2","case.cfg","debugger_process.json","native.json","debugger.log",
                     "readout/accepted.json",f"restart_{n:05d}.csv",f"history_{n:05d}.csv",f"seed_{n-1:05d}.csv"))
        paths.extend(ROOT/"cases/short_sphere"/("wake_"+branch+ext) for ext in (".json",".zip"))
    locks = {p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    target = out/"wakeResolution.json"
    invoke("resolution", [sys.executable,script,work,target], ROOT, out, 120, 2)
    assert all(sha(ROOT/p)==h for p,h in locks.items())
    assert sum(p.stat().st_size for p in out.iterdir()) < 100*1024**2
    write_json(out/"inputs.json",dict(inputs=locks,output=sha(target)))
