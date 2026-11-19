"""在 66 mm 加热圆柱相同初态上运行一次原生自适应对照，记录完整行出现的时间。"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "v279"))
from tunnel_step import EXE, EXE_SHA, save, sha

OLD = ROOT / "outputs/v279_validation/flow04"
RUN = ROOT / "outputs/v280_validation/adapt01"
INPUT_MANIFEST = "4398894cdef2cfcbd060d1090b1dcdd0cf3742bc64a11ddb82cb91279e5ed551"


def prepare():
    manifest = OLD / "inputs.json"
    assert sha(manifest) == INPUT_MANIFEST and sha(EXE) == EXE_SHA
    known = {Path(k).name: v for k, v in json.loads(manifest.read_bytes()).items()}
    RUN.mkdir(parents=True, exist_ok=False)
    bound = [manifest, EXE, HERE / "PLAN.md", Path(__file__), HERE.parent / "v279/tunnel_step.py",
             ROOT / "scripts/runEmbeddedSU2.py", OLD / "case.cfg"]
    for name in ("mesh.su2", "seed_00000.csv", "condition.json", "wall.npz"):
        assert sha(OLD / name) == known[name]
        target = RUN / name
        os.link(OLD / name, target)
        bound.append(target)
    assert sha(OLD / "case.cfg") == known["case.cfg"]
    cfg = (OLD / "case.cfg").read_text("ascii")
    probes = {"Pfront": "0,.033,0", "Pside": "0,0,.033", "Prear": "0,-.033,0", "Pref": "0,.33,0"}
    expressions = "; ".join(f"{key} : Probe{{PRESSURE}}[{xyz}]" for key, xyz in probes.items())
    changes = {"INNER_ITER= 120": "INNER_ITER= 720",
               "CFL_NUMBER= 3": "CFL_NUMBER= 3\nCFL_ADAPT= YES\nCFL_ADAPT_PARAM= (0.5, 1.15, 1, 48, 0.005, 8)",
               "HEAT_FLUX, Y_PLUS)": "HEAT_FLUX, Y_PLUS, RESIDUAL)",
               "HISTORY_OUTPUT= (ITER, RMS_RES, HEAT)":
                   f"CUSTOM_OUTPUTS= '{expressions}'\nHISTORY_OUTPUT= (ITER, RMS_RES, HEAT, MAX_RES, LINSOL, CFL_NUMBER, CUSTOM)"}
    for before, after in changes.items():
        assert cfg.count(before) == 1
        cfg = cfg.replace(before, after)
    target = RUN / "case.cfg"
    target.write_text(cfg, encoding="ascii")
    bound.append(target)
    save(RUN / "inputs.json", {str(p): sha(p) for p in bound})


def last_visible():
    """只认落盘且有换行的完整数值行，返回从 1 开始的已完成次数。"""
    try:
        with (RUN / "history_00001.csv").open("rb") as f:
            f.seek(max(0, f.seek(0, 2) - 8192))
            lines = f.read().split(b"\n")
        if len(lines) >= 2:
            cols = lines[-2].strip().split(b",")
            if len(cols) >= 3 and cols[2].strip().isdigit():
                return int(cols[2]) + 1
    except (FileNotFoundError, PermissionError):
        pass
    return 0


def solve():
    bindings = json.loads((RUN / "inputs.json").read_bytes())
    assert all(sha(Path(p)) == h for p, h in bindings.items())
    command = [str(EXE), "case.cfg", "--threads", "2"]
    start = time.monotonic()
    peak = disk = 0
    reason = None
    with (RUN / "native.log").open("xb") as log, (RUN / "visible.csv").open("x") as visible:
        visible.write("elapsedS,completeRows\n")
        proc = subprocess.Popen(command, cwd=RUN, stdout=log, stderr=log,
                                creationflags=subprocess.CREATE_NO_WINDOW,
                                env=dict(os.environ, OMP_NUM_THREADS="2"))
        try:
            handle = psutil.Process(proc.pid)
            save(RUN / "started.json", dict(pid=proc.pid, created=handle.create_time(), command=command,
                                            inputsSHA256=sha(RUN / "inputs.json")))
            while proc.poll() is None:
                try:
                    peak = max(peak, handle.memory_info().rss)
                except psutil.NoSuchProcess:
                    break
                disk = max(disk, sum(p.stat().st_size for p in RUN.iterdir() if p.is_file()))
                count = last_visible()
                elapsed = time.monotonic() - start
                visible.write(f"{elapsed:.6f},{count}\n")
                visible.flush()
                if elapsed >= 420 or peak > 2 * 1024**3 or disk > 240 * 1024**2:
                    reason = "resource budget"
                    break
                time.sleep(.2)
        finally:
            if proc.poll() is None:
                proc.kill()
            proc.wait()
        elapsed = time.monotonic() - start
        visible.write(f"{elapsed:.6f},{last_visible()}\n")
    disk = max(disk, sum(p.stat().st_size for p in RUN.iterdir() if p.is_file()))
    if disk > 240 * 1024**2:
        reason = "final output exceeds disk budget"
    stable = all(sha(Path(p)) == h for p, h in bindings.items())
    save(RUN / "terminal.json", dict(pid=proc.pid, command=command, exitCode=proc.returncode,
         elapsedSeconds=elapsed, stopReason=reason, inputsStable=stable,
         peakSampledRssBytes=peak, peakSampledFolderBytes=disk,
         outputs={p.name: sha(p) for p in RUN.iterdir() if p.is_file() and str(p) not in bindings}))
    if reason or proc.returncode or not stable:
        raise RuntimeError("请核对终态回执和原生日志")


if __name__ == "__main__":
    prepare()
    solve()
