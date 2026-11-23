"""在真实任务被占用时验证并发接续和旧温差任务均被拒绝，保留已有进程与初场。"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from thermistorBead import read, sha, write_json


def check(run, prepared, result):
    progress = read(run / "progress.json")
    stage = run / progress["stage"]
    assert progress["status"] == "running"
    receipt = read(stage / "engine_process.json")
    assert receipt["status"] == "running"
    before = sorted(p.name for p in run.glob("steady*"))
    original = sha(run / "cold/initial_flow.csv")
    def reject(folder, expected):
        command = [sys.executable, str(ROOT / "scripts/runThermistorMeasuredCold.py"), "resume", str(folder),
                   "--prepared", str(prepared)]
        r = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=30,
                           env=dict(os.environ,PYTHONIOENCODING="utf-8"))
        error = r.stderr.decode("utf8")
        assert r.returncode != 0 and expected in error, error
        return error.strip().splitlines()[-1]
    occupied = reject(run, "该目录仍有本任务运行")
    wrong = reject(ROOT / "outputs/v315_validation/native01", "所选目录没有本工况的冷启动任务记录")
    after = read(stage / "engine_process.json")
    assert receipt["processID"] == after["processID"] and after["status"] == "running"
    assert sha(run / "cold/initial_flow.csv") == original and before == sorted(p.name for p in run.glob("steady*"))
    write_json(result, dict(occupied=occupied, wrongJob=wrong, processID=after["processID"], stage=stage.name,
               initialSHA256=original, sourceJob=progress["job"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("run","prepared","result"):p.add_argument(name,type=Path)
    args = p.parse_args()
    check(args.run.resolve(),args.prepared.resolve(),args.result.resolve())
