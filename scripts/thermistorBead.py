"""坪内、佐藤 B 型热敏珠的名义等温珠体，低速空气单工况。"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

import numpy as np
import psutil

from cm011Tunnel import laminar_config, surface, volume, write_json, write_volume_mesh
from runEmbeddedSU2 import sha

ROOT = Path(__file__).resolve().parents[1]
AREA, ECC = 3.094e-6, .53
RATIO = math.sqrt(1-ECC**2)
A = math.sqrt(AREA / (2*math.pi*RATIO**2*(1+math.asin(ECC)/(RATIO*ECC))))
B = A*RATIO
D = 2*B
U, RHO, MU, CP, K = 3., 1.180, 1.842e-5, 1006., .0262
P0, T0, DT = 101325., 293.15, 12.
RE, PE = RHO*U*D/MU, RHO*CP*U*D/K
GRIDS = {"coarse": (40, 30, 12., 1), "fine": (64, 42, 12., 2), "outer": (64, 48, 20., 2)}
ENGINE = ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
ENGINE_SHA = "01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"


def wall(stations):
    angle = np.linspace(math.pi, 0, stations+1)
    w = np.column_stack((A/D*np.cos(angle), B/D*np.sin(angle)))
    w[[0,-1],1] = 0
    return w


def body():
    # 三套体网格共用同一显示面。保留从三角顶点到子午极角的对应。
    w = wall(64)
    return surface(w, np.zeros(64,int), 1, D, ["热敏珠完整外表面"])


def grid_points(grid):
    m, n, radius, level = GRIDS[grid]
    return volume(wall(m), n, radius, level)


def read(path):
    return json.loads(path.read_bytes())


def prepare(folder, grid):
    assert sha(ENGINE) == ENGINE_SHA
    folder.mkdir(parents=True, exist_ok=False)
    xy, cells, markers = grid_points(grid)
    write_volume_mesh(folder/"mesh.su2",xy,cells,markers)
    write_json(folder/"body.json",body())
    cfg = laminar_config(RE, PE, 36000).replace("OUTPUT_WRT_FREQ= (1500, 500, 18000)",
                                               "OUTPUT_WRT_FREQ= (2000, 2000, 36000)")
    cfg = cfg.replace("MUSCL_FLOW= YES", "MUSCL_FLOW= NO")
    m, n, extent, _ = GRIDS[grid]
    spec = dict(grid=grid, areaSourceM2=AREA, eccentricity=ECC, semiMajorM=A, semiMinorM=B,
                diameterM=D, measuredLengthM=.3430e-3, impliedLengthM=8*math.pi*A*B*B/(3*AREA),
                freestreamMps=U, pressurePa=P0, ambientK=T0, wallK=T0+DT, rho=RHO, mu=MU, cp=CP, k=K,
                reynolds=RE, peclet=PE, richardson=9.80665*DT*D/((T0+DT/2)*U**2),
                meridianIntervals=m, radialIntervals=n, outerDiameters=extent)
    np.savez_compressed(folder/"initial.npz",xy=xy,pressure=np.zeros(len(xy)),velocity=np.tile([-1.,0.],(len(xy),1)),
                        temperature=np.ones(len(xy)),wallNodes=np.unique(markers["wall"]),wallTemperature=2.)
    (folder/"case.cfg").write_text(cfg,encoding="ascii")
    write_json(folder/"case.json",spec)
    inputs = [Path(__file__).resolve(),ROOT/"scripts/cm011Tunnel.py",folder/"case.json",folder/"body.json",
              folder/"initial.npz"]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               threads=1,maximumSeconds=360,jobToken=str(uuid.uuid4()),
               sourceFiles=[p.relative_to(ROOT).as_posix() for p in inputs]))


def solve(folder):
    request = folder/"engine_request.json"
    # 工作台关闭时停止其专属运行。单独启动阶段脚本时仍沿用原有监督方式。
    flag = os.environ.get("RAH_THERMISTOR_CANCEL_FILE")
    if flag and Path(flag).exists():
        raise RuntimeError("用户已关闭本次热敏珠计算窗口")
    cmd = [sys.executable,str(ROOT/"scripts/runEmbeddedSU2.py"),"--request",str(request)]
    with (folder/"supervisor.log").open("xb") as f:
        proc = subprocess.Popen(cmd, stdout=f, stderr=f, creationflags=subprocess.CREATE_NO_WINDOW)
        handle = psutil.Process(proc.pid)
        peak = size = 0
        exceed = False
        try:
            while proc.poll() is None:
                try:
                    family = [handle]+handle.children(recursive=True)
                    peak = max(peak,sum(p.memory_info().rss for p in family if p.is_running()))
                except psutil.NoSuchProcess:
                    pass
                size = max(size,sum(p.stat().st_size for p in folder.iterdir() if p.is_file()))
                if not exceed and (peak > 2*1024**3 or size > 180*1024**2):
                    exceed = True
                    write_json(folder/"engine_cancel.json",dict(jobToken=read(request)["jobToken"]))
                if flag and Path(flag).exists() and not (folder/"engine_cancel.json").exists():
                    write_json(folder/"engine_cancel.json",dict(jobToken=read(request)["jobToken"]))
                time.sleep(.25)
        finally:
            if proc.poll() is None and not (folder/"engine_cancel.json").exists():
                write_json(folder/"engine_cancel.json",dict(jobToken=read(request)["jobToken"]))
            proc.wait()
    size = max(size,sum(p.stat().st_size for p in folder.iterdir() if p.is_file()))
    write_json(folder/"resources.json",dict(peakSampledFamilyRSSBytes=peak,directoryBytes=size,
                                            resourceExceeded=bool(exceed or size>180*1024**2)))
    if proc.returncode or exceed or size>180*1024**2 or read(folder/"engine_process.json")["cancelled"]:
        raise RuntimeError("热敏珠原生计算结束，请核对原生回执")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("grid",choices=GRIDS)
    args=parser.parse_args()
    folder=args.folder.resolve()
    folder.relative_to(ROOT/"outputs")
    prepare(folder,args.grid)
    solve(folder)
