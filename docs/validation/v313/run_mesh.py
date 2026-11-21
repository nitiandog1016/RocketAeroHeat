"""监督热敏珠细网的建网或独立读回，保存命令、耗时、资源与实际终态。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def watch(action, folder, base, seed):
    if action == "prepare":
        script = ROOT / "scripts/thermistorLeadTipSpacing.py"
        command = [sys.executable, str(script), str(base), str(seed), str(folder)]
    else:
        command = [sys.executable, str(HERE / "audit_mesh.py"), str(folder)]
    record = folder.with_name(folder.name + "_" + action + ".json")
    record.parent.mkdir(parents=True, exist_ok=True)
    assert not record.exists()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    start = time.monotonic()
    peak = size = 0
    stop = None
    with record.with_suffix(".log").open("xb") as log:
        p = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log, stderr=log)
        handle = psutil.Process(p.pid)
        while p.poll() is None:
            try:
                peak = max(peak, handle.memory_info().rss)
            except psutil.NoSuchProcess:
                pass
            size = max(size, sum(v.stat().st_size for v in folder.glob("*") if v.is_file()))
            if time.monotonic()-start > 180 or peak > 3*1024**3 or size > 300*1024**2:
                stop = "resource limit"
                p.kill()
                break
            time.sleep(.15)
        code = p.wait()
    duration = time.monotonic()-start
    size = max(size, sum(v.stat().st_size for v in folder.glob("*") if v.is_file()))
    result = dict(command=command, pid=p.pid, status="terminal", exitCode=code, stopReason=stop,
                  seconds=duration, peakSampledRSSBytes=peak, outputBytes=size)
    record.write_text(json.dumps(result, indent=2), encoding="utf8")
    assert code == 0 and not stop and size <= 300*1024**2


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "audit"])
    p.add_argument("folder", type=Path)
    p.add_argument("--base", type=Path, default=ROOT / "outputs/v311_validation/cold01")
    p.add_argument("--seed", type=Path, default=ROOT / "outputs/v311_validation/steady02")
    args = p.parse_args()
    folder = args.folder.resolve()
    folder.relative_to(ROOT / "outputs")
    watch(args.action, folder, args.base.resolve(), args.seed.resolve())
