"""简化运载火箭米制模型的固定低速风洞算例。"""
import argparse
import json
from pathlib import Path

import numpy as np

from cm011Tunnel import laminar_config, volume, write_json, write_volume_mesh
from runEmbeddedSU2 import ROOT, sha

CASE = ROOT / "cases/simple_launcher"
SPEED, RHO, MU, CP, K = 0.05, 1.177, 1.846e-5, 1007.0, 0.02624
T_AIR, T_WALL, P_REF = 300.0, 300.001, 101325.0
GRIDS = {"coarse": (1, 36, 14), "fine": (2, 52, 14), "outer": (2, 56, 22)}


def geometry():
    """按原 STL 面序读取，抽取圆环与底面径向母线。"""
    definition = json.loads((CASE / "geometry.json").read_bytes())
    path = CASE / definition["file"]
    if sha(path) != definition["fileSHA256"]:
        raise ValueError("当前 STL 与已登记的几何不一致")
    dtype = np.dtype([("normal", "<f4", 3), ("corners", "<f4", (3, 3)), ("attribute", "<u2")])
    triangles = np.fromfile(path, dtype=dtype, offset=84)["corners"].astype(float)
    if len(triangles) != definition["faces"]:
        raise ValueError("STL 面数与几何记录不一致")
    vertices = np.unique(triangles.reshape(-1, 3), axis=0)
    xs = np.unique(vertices[:, 0])
    radii = np.hypot(vertices[:, 1], vertices[:, 2])
    radius = np.array([radii[vertices[:, 0] == x].max() for x in xs])
    diameter = 2 * radius.max()
    base = np.column_stack((np.zeros(13), np.linspace(0., radius[0], 13)))
    wall = np.vstack((base, np.column_stack((xs[1:], radius[1:])))) / diameter
    return definition, triangles, diameter, wall


def prepare(folder, grid):
    folder = folder.resolve()
    folder.relative_to(ROOT / "outputs")
    level, layers, extent = GRIDS[grid]
    definition, _, diameter, wall = geometry()
    refined = np.concatenate([a + (b-a) * np.arange(level)[:, None] / level
                              for a, b in zip(wall[:-1], wall[1:])] + [wall[-1:]])
    xy, cells, markers = volume(refined, layers, extent, level)
    folder.mkdir(parents=True, exist_ok=False)
    write_volume_mesh(folder / "mesh.su2", xy, cells, markers)
    re = RHO*SPEED*diameter/MU
    pe = RHO*CP*SPEED*diameter/K
    config = laminar_config(re, pe, 16000)
    config = config.replace("OUTPUT_WRT_FREQ= (1500, 500, 18000)", "OUTPUT_WRT_FREQ= (500, 500, 16000)")
    (folder / "case.cfg").write_text(config, encoding="ascii")
    np.savez_compressed(folder / "initial.npz", xy=xy, pressure=np.zeros(len(xy)),
                        velocity=np.tile([-1., 0.], (len(xy), 1)), theta=np.ones(len(xy)),
                        wallNodes=np.unique(markers["wall"]), wallTheta=2.)
    write_json(folder / "condition.json", dict(grid=grid, layers=layers, extentD=extent,
               surfaceSHA256=definition["fileSHA256"], diameterM=diameter, reynoldsD=re, pecletD=pe,
               speedMps=SPEED, airK=T_AIR, wallK=T_WALL, pRefPa=P_REF, wall=refined.tolist()))
    sources=[CASE / "geometry.json", CASE / definition["file"], Path(__file__), Path(__file__).with_name("cm011Tunnel.py"),
             folder / "condition.json", folder / "initial.npz"]
    write_json(folder / "engine_request.json", dict(workDirectory=str(folder), configFile="case.cfg",
               meshFile="mesh.su2", maximumSeconds=720, threads=1,
               sourceFiles=[p.resolve().relative_to(ROOT).as_posix() for p in sources]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("grid", choices=GRIDS)
    args = parser.parse_args()
    prepare(args.folder, args.grid)
