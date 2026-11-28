"""按登记的 M0.30 时间步配对运行固定版 SU2，并逐步绑定实际的起止状态。"""
import argparse
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy.spatial import cKDTree

DOC = Path(__file__).resolve().parent
ROOT = DOC.parents[2]
BASE = ROOT / "outputs/v369_validation"
sys.path.insert(0, str(DOC.parent / "v366"))
import physical_step as old
from energy_budget import budget
from physical_step import EXE, EXE_SHA, OLD, ORIGIN, RES, numbers, observe, read, save, sha, terminal, vtk_tokens

TEMPLATE = ROOT / "outputs/v366_validation/time02/full/case.cfg"
GRID_HASH = "f5656636b80e5a47a7d0cdb04b52f7a97e24fefb274f298ae1c31f3a42e4da68"
SEED_HASH = "998600a5f54616298783e53f78370ecd51b203f64e639835f598d4d8a7871346"
STEPS = {"full": (12.5e-6, 1), "half1": (6.25e-6, 1), "half2": (6.25e-6, 2)}


def configuration(kind):
    dt, index = STEPS[kind]
    assert sha(TEMPLATE) == "2988b783e0eb9583d12f8f598b662ee994f1b60da3a954e6193a1a68e911563b"
    cfg = old.config(TEMPLATE)
    assert cfg == old.step_config("fixed", "full")
    cfg.update(TIME_STEP=str(dt), RESTART_ITER=str(index), TIME_ITER=str(index+1))
    return cfg


def residual_window(values):
    """在密度、两动量和能量上共同检查连续 20 行，保留每个窗口最差的分量。"""
    assert values.ndim == 2 and values.shape[1] == 4 and len(values) >= 20 and np.isfinite(values).all()
    worst = values.max(axis=1)
    windows = np.lib.stride_tricks.sliding_window_view(worst, 20).max(axis=1)
    best = int(np.argmin(windows))
    last = values[-20:].max(axis=0)
    return dict(iterations=len(values), last20Max=dict(zip(RES,last.tolist())), passed=bool(max(last)<-8),
                bestWindowFirstIter=best, bestWindowWorst=float(windows[best]),
                passingWindows=int(np.count_nonzero(windows < -8)))


def inspect(folder, kind):
    dt, index = STEPS[kind]
    proc = terminal(folder)
    assert all(sha(Path(path)) == h for path,h in proc["inputSHA256"].items())
    assert old.config(folder/"case.cfg") == configuration(kind)
    c = read(OLD/"condition.json")
    seed = folder / f"seed_{index-1:05d}.csv"
    frame = folder / f"restart_{index:05d}.csv"
    wall = folder / f"surface_{index:05d}.vtk"
    history = folder / f"history_{index:05d}.csv"
    h = numbers(history)
    np.testing.assert_array_equal(h["Inner_Iter"], np.arange(len(h["Inner_Iter"])))
    assert len(h["Inner_Iter"]) <= 2500 and (h["Time_Iter"] == index).all()
    assert (h["Time_Step"] == dt).all() and np.unique(h["Cur_Time"]).size == 1
    assert min(abs(h["Cur_Time"][0] - t) for t in (dt, dt*index)) < 1e-12
    assert {k for k in h if k.startswith("rms[")} == set(RES)
    values = np.column_stack([h[k] for k in RES])
    result = residual_window(values)
    a, b = numbers(frame), numbers(seed)
    assert len(a["x"]) == 9701 and np.isfinite(np.stack(list(a.values()))).all()
    for key in ("x", "y", "PointID"):
        np.testing.assert_array_equal(a[key], b[key])
    rho = a["Density"]
    pressure = .4 * (a["Energy"] - .5*(a["Momentum_x"]**2+a["Momentum_y"]**2)/rho)
    np.testing.assert_allclose(pressure, a["Pressure"], rtol=1e-10)
    np.testing.assert_allclose(pressure/rho, a["Temperature"], rtol=1e-10)
    assert min(rho)>0 and min(pressure)>0 and np.any(a["Energy"] != b["Energy"])
    assert not np.any(a["Momentum_x"][:109]) and not np.any(a["Momentum_y"][:109])
    xyz, fields = vtk_tokens(wall)
    distance, order = cKDTree(xyz[:,:2]).query(np.column_stack((a["x"][:109], a["y"][:109])))
    assert len(xyz) == len(np.unique(order)) == 109 and max(distance)<1e-6
    for name, key in (("Pressure","Pressure"), ("Heat_Flux","Heat_Flux"), ("Temperature","Temperature"),
                      ("Skin_Friction_Coefficient","Skin_Friction_Coefficient_x")):
        np.testing.assert_allclose(fields[name][order,0], a[key][:109], rtol=6e-6, atol=1e-9)
    source = ORIGIN if index==1 else folder.parent/"half1/restart_00001.csv"
    assert sha(seed) == sha(source) and sha(folder/"mesh.su2") == GRID_HASH
    result.update(kind=kind, deltaUs=dt*1e6, nominalEndUs=index*dt*1e6, nativeTimeS=float(h["Cur_Time"][0]),
                  meshSHA256=GRID_HASH, sourceSHA256=sha(seed), frameSHA256=sha(frame), surfaceSHA256=sha(wall),
                  cfgSHA256=sha(folder/"case.cfg"), historySHA256=sha(history), seconds=proc["elapsedSeconds"],
                  peakRSSBytes=proc["peakFamilyRSSBytes"], fieldBalance=budget(b,a,dt,c))
    return result


def execute(out, kind):
    out = out.resolve()
    out.relative_to(BASE)
    terminal(OLD)
    assert sha(EXE)==EXE_SHA and sha(ORIGIN)==SEED_HASH and sha(OLD/"mesh.su2")==GRID_HASH
    assert sha(OLD/"condition.json")=="712d77252229a50d27592266654e0e819b18b3e57d34c23adcf2ad00a0384cf1"
    cfg = configuration(kind)
    source = ORIGIN
    if kind == "half2":
        for prior in ("full","half1"):
            r = inspect(out/prior, prior)
            if not r["passed"] or r["sourceSHA256"]!=SEED_HASH:
                raise ValueError(f"{prior} 未通过，无法接续第二半步")
        source = out/"half1/restart_00001.csv"
    folder = out/kind
    folder.mkdir(parents=True, exist_ok=False)
    index = STEPS[kind][1]
    config_file = folder/"case.cfg"
    config_file.write_text("".join(f"{k}= {v}\n" for k,v in cfg.items()), encoding="ascii")
    seed, mesh = folder/f"seed_{index-1:05d}.csv", folder/"mesh.su2"
    shutil.copyfile(source, seed)
    shutil.copyfile(OLD/"mesh.su2", mesh)
    assert sha(seed)==sha(source) and sha(mesh)==GRID_HASH
    locked = (EXE, seed, mesh, source, config_file, TEMPLATE, OLD/"condition.json", DOC/"PLAN.md",
              Path(__file__), Path(old.__file__), ROOT/"scripts/runSU2NativeDebugger.py")
    observe([str(EXE), config_file.name, "--threads", "2"], folder, {str(p):sha(p) for p in locked}, 150,64*1024**2)
    r = inspect(folder, kind)
    save(folder/"step.json",r)
    print(kind,r["last20Max"],r["passed"],flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("kind",choices=STEPS)
    a = p.parse_args()
    execute(a.directory,a.kind)
