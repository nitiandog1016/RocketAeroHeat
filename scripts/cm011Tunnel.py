"""CM-011 缩比等温空气试件的母线、体网格及固定 SU2 工况。"""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq

from runEmbeddedSU2 import sha

ROOT = Path(__file__).resolve().parents[1]
D = 155.481 * 0.0254 / 1000
RHO, U, MU, CP, K = 1.177, 0.1, 1.846e-5, 1007., 0.02624
P0, T0, DT = 101325., 300., 0.1
RE, PE = RHO * U * D / MU, RHO * CP * U * D / K
GRIDS = {"coarse": (1, 24, 12.), "fine": (2, 42, 12.), "outer": (2, 48, 20.)}
REGIONS = ["小端球帽", "锥面", "肩部", "迎风盾面"]


def write_json(path, value):
    with path.open("x", encoding="utf8") as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def meridian(level):
    # 所有长度先除以试件直径。由球面与圆角内切、锥面与圆弧相切确定连接点。
    front, tail, shoulder = np.array([180.890, 9.0118, (7.46+7.17)/2]) / 155.481
    alpha = math.radians(32+31/60+23/3600)
    join = math.asin((.5-shoulder)/(front-shoulder))
    cone = math.pi/2 + alpha
    tail_end = np.array([tail*(1+math.cos(cone)), tail*math.sin(cone)])
    rc = .5-shoulder
    xc = tail_end[0] + (rc+shoulder*math.sin(cone)-tail_end[1])/math.tan(alpha) - shoulder*math.cos(cone)
    fc = xc-(front-shoulder)*math.cos(join)
    points, strips = [], []

    def append(coords, tag):
        if points and not np.allclose(points[-1], coords[0], atol=1e-13, rtol=0):
            raise ValueError("母线接缝不闭合")
        points.extend(coords if not points else coords[1:])
        strips.extend([tag]*(len(coords)-1))

    def arc(cx, cr, radius, a, b, count, tag):
        t = np.linspace(a, b, count*level+1)
        append(np.column_stack((cx+radius*np.cos(t), cr+radius*np.sin(t))), tag)

    arc(tail, 0, tail, math.pi, cone, 5, 0)
    append(np.linspace(tail_end, [xc+shoulder*math.cos(cone), rc+shoulder*math.sin(cone)], 18*level+1), 1)
    arc(xc, rc, shoulder, cone, math.pi/2, 4, 2)
    arc(xc, rc, shoulder, math.pi/2, join, 7, 2)
    arc(fc, 0, front, join, 0, 14, 3)
    w = np.array(points)
    w[[0, -1], 1] = 0
    return w, np.array(strips)


def surface(w, bands, level, diameter=D, region_names=REGIONS):
    # 不同半径上用不同环点数；面编号与母线源节点一同返回给 MATLAB。
    vertices, rings, source, faces, tags = [], [], [], [], []
    for i, (x, r) in enumerate(w):
        count = 1 if r == 0 else max(6, round(64*level*r/.5))
        rings.append(np.arange(len(vertices), len(vertices)+count))
        angles = np.arange(count)*2*np.pi/count
        vertices.extend(np.column_stack((np.full(count, x), r*np.cos(angles), r*np.sin(angles))))
        source.extend([i]*count)
    for k, (a, b) in enumerate(zip(rings, rings[1:])):
        na, nb = len(a), len(b)
        if na == 1:
            strip = [(a[0], b[(j+1) % nb], b[j]) for j in range(nb)]
        elif nb == 1:
            strip = [(a[i], a[(i+1) % na], b[0]) for i in range(na)]
        else:
            strip, i, j = [], 0, 0
            while i < na or j < nb:
                if i < na and (j == nb or (i+1)*nb <= (j+1)*na):
                    strip.append((a[i % na], a[(i+1) % na], b[j % nb]))
                    i += 1
                else:
                    strip.append((a[i % na], b[(j+1) % nb], b[j % nb]))
                    j += 1
        faces.extend(strip)
        tags.extend([int(bands[k])]*len(strip))
    v, f = np.array(vertices)*diameter, np.array(faces)
    signed = np.einsum("ij,ij->i", v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]]))
    if signed.sum() <= 0:
        raise ValueError("三维试件表面朝向错误")
    return dict(V=v.tolist(), F=(f+1).tolist(), nodeRow=source, regionCode=tags, regionNames=region_names)


def volume(w, layers, extent, level):
    center = np.array([(w[0, 0]+w[-1, 0])*.5, 0])
    ray = w-center
    length = np.linalg.norm(ray, axis=1)
    direction = ray/length[:, None]
    # 逐射线保持首层 0.004D/level，向外按等比数列扩张。
    h = .004/level
    coords = []
    for wall, unit, far in zip(w, direction, extent-length):
        ratios = np.arange(layers)
        growth = brentq(lambda a: h*np.sum(a**ratios)-far, 1, 2)
        dist = np.r_[0., np.cumsum(h*growth**ratios)]
        coords.append(wall+dist[:, None]*unit)
    pts = np.array(coords).transpose(1, 0, 2)
    m = len(w)
    ids = np.arange(pts.shape[0]*m).reshape(-1, m)
    cells = np.stack((ids[:-1, :-1], ids[:-1, 1:], ids[1:, 1:], ids[1:, :-1]), axis=-1).reshape(-1, 4)
    uv = pts.reshape(-1, 2)
    corners = uv[cells]
    edge = np.roll(corners, -1, axis=1)-corners
    if not np.all(np.cross(edge, np.roll(edge, -1, axis=1)) > 0):
        raise ValueError("流体单元翻转")
    axis = np.vstack((np.column_stack((ids[1:, 0], ids[:-1, 0])), np.column_stack((ids[:-1, -1], ids[1:, -1]))))
    boundaries = {"wall": np.column_stack((ids[0, :-1], ids[0, 1:])),
                  "far": np.column_stack((ids[-1, 1:], ids[-1, :-1])), "axis": axis}
    return uv, cells, boundaries


def extend_initial_field(previous, folder, w, xy, layers):
    # 扩域只复用细网稳态作初猜。每条母线射线独立插值，远场外补初始均匀流。
    from readCM011Tunnel import evaluate, match_nodes
    summary, _ = evaluate(previous)
    if not summary["acceptedNumericalWall"] or summary["grid"] != "fine":
        raise ValueError("扩大外域需要已核对的细网数值终场")
    data = np.loadtxt(previous/"restart.csv", delimiter=",", skiprows=1)
    old_grid, _, _ = volume(w, GRIDS["fine"][1], GRIDS["fine"][2], 2)
    order = match_nodes(data[:, 1:3], old_grid)
    old = data[order, 3:].reshape(-1, len(w), 4)
    old_r = np.linalg.norm(old_grid.reshape(-1, len(w), 2)-w, axis=2)
    new_r = np.linalg.norm(xy.reshape(-1, len(w), 2)-w, axis=2)
    # 固定版本 CONSTANT_DENSITY 的 h=cp*(T-298.15)。冷启动远场归一化 T=1。
    far = [0., -1., 0., 1.-298.15]
    new = np.empty((layers+1, len(w), 4))
    for i in range(len(w)):
        for k, initial in enumerate(far):
            new[:, i, k] = np.interp(new_r[:, i], old_r[:, i], old[:, i, k], right=initial)
    np.savetxt(folder/"initial_flow.csv", np.column_stack((np.arange(len(xy)), xy, new.reshape(-1, 4))),
               fmt=["%d"]+["%.17e"]*6, delimiter=",", comments="",
               header='"PointID","x","y","Pressure","Velocity_x","Velocity_y","Enthalpy"')


def write_volume_mesh(path, xy, cells, markers):
    with path.open("x", encoding="ascii") as stream:
        stream.write(f"NDIME= 2\nNELEM= {len(cells)}\n")
        np.savetxt(stream, cells, fmt="9 %d %d %d %d")
        stream.write(f"NPOIN= {len(xy)}\n")
        np.savetxt(stream, xy, fmt="%.17g")
        stream.write(f"NMARK= {len(markers)}\n")
        for name, edges in markers.items():
            stream.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(edges)}\n")
            np.savetxt(stream, edges, fmt="3 %d %d")


def laminar_config(reynolds, peclet, iterations):
    return f"""SOLVER= INC_NAVIER_STOKES
AXISYMMETRIC= YES
MATH_PROBLEM= DIRECT
INC_DENSITY_MODEL= CONSTANT
INC_ENERGY_EQUATION= YES
INC_NONDIM= DIMENSIONAL
FLUID_MODEL= CONSTANT_DENSITY
INC_DENSITY_INIT= 1
INC_VELOCITY_INIT= (-1, 0, 0)
INC_TEMPERATURE_INIT= 1
SPECIFIC_HEAT_CP= 1
VISCOSITY_MODEL= CONSTANT_VISCOSITY
MU_CONSTANT= {1/reynolds:.17g}
CONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY
THERMAL_CONDUCTIVITY_CONSTANT= {1/peclet:.17g}
MARKER_ISOTHERMAL= (wall, 2)
MARKER_SYM= (axis)
MARKER_FAR= (far)
MARKER_MONITORING= (wall)
MARKER_PLOTTING= (wall)
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CONV_NUM_METHOD_FLOW= FDS
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= NONE
TIME_DISCRE_FLOW= EULER_IMPLICIT
CFL_NUMBER= 10
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 0.001
LINEAR_SOLVER_ITER= 20
ITER= {iterations}
CONV_FIELD= (RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_ENTHALPY)
CONV_RESIDUAL_MINVAL= -11
CONV_STARTITER= 50
MESH_FILENAME= mesh.su2
MESH_FORMAT= SU2
SCREEN_OUTPUT= (INNER_ITER, RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_ENTHALPY, TOTAL_HEATFLUX)
SCREEN_WRT_FREQ_INNER= 100
HISTORY_OUTPUT= (ITER, RMS_RES, HEAT)
OUTPUT_FILES= (RESTART_ASCII, SURFACE_PARAVIEW_ASCII, PARAVIEW_ASCII)
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
OUTPUT_WRT_FREQ= (1500, 500, 18000)
WRT_RESTART_OVERWRITE= NO
WRT_SURFACE_OVERWRITE= NO
TABULAR_FORMAT= CSV
"""


def prepare(folder, grid, restart=None, seed=None):
    folder.mkdir(parents=True, exist_ok=False)
    level, layers, extent = GRIDS[grid]
    w, bands = meridian(level)
    surf = surface(w, bands, level)
    surf.update(meridian=w.tolist(), diameterM=D, grid=grid, layers=layers, extentD=extent)
    write_json(folder/"body.json", surf)
    xy, cells, markers = volume(w, layers, extent, level)
    write_volume_mesh(folder/"mesh.su2", xy, cells, markers)
    config = laminar_config(RE, PE, 10000 if grid == "coarse" else 18000)
    if restart is not None:
        previous = restart.resolve()
        receipt = json.loads((previous/"engine_process.json").read_bytes())
        if receipt["status"] != "terminal" or receipt["exitCode"] != 0 or not all(receipt["inputsStable"].values()):
            raise ValueError("上一次求解终态身份未确认")
        if sha(previous/"mesh.su2") != sha(folder/"mesh.su2") or sha(previous/"body.json") != sha(folder/"body.json"):
            raise ValueError("重新启动要求完全相同的体网格与试件")
        old_config = (previous/"case.cfg").read_text("ascii").partition("CFL_NUMBER=")[0]
        if old_config != config.partition("CFL_NUMBER=")[0]:
            raise ValueError("重新启动的方程或物理边界不同")
        if sha(previous/"restart.csv") != receipt["outputSHA256"]["restart.csv"]:
            raise ValueError("重启动场与原终态回执不符")
        shutil.copyfile(previous/"restart.csv", folder/"initial_flow.csv")
    if seed is not None:
        if grid != "outer" or restart is not None:
            raise ValueError("扩域初猜只用于同工况的外域对照")
        extend_initial_field(seed.resolve(), folder, w, xy, layers)
    if restart is not None or seed is not None:
        config += "RESTART_SOL= YES\nREAD_BINARY_RESTART= NO\nSOLUTION_FILENAME= initial_flow\nRESTART_FILENAME= restart\n"
    (folder/"case.cfg").write_text(config, encoding="ascii")
    # 保存指定的完整均匀初猜与壁面约束，原生首个已推进帧由 SU2 自行写出。
    if restart is None and seed is None:
        np.savez_compressed(folder/"initial.npz", xy=xy, pressure=np.zeros(len(xy)),
                            velocity=np.tile([-1., 0.], (len(xy), 1)), theta=np.ones(len(xy)),
                            wallNodes=markers["wall"].ravel(), wallTheta=2.)
    request = dict(workDirectory=str(folder.resolve()), configFile="case.cfg", meshFile="mesh.su2",
                   maximumSeconds=180, threads=1, sourceFiles=[Path(__file__).resolve().relative_to(ROOT).as_posix(),
                   (folder/"body.json").resolve().relative_to(ROOT).as_posix()])
    write_json(folder/"engine_request.json", request)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("grid", choices=GRIDS)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--seed", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    folder = args.folder.resolve()
    folder.relative_to(ROOT/"outputs")
    prepare(folder, args.grid, args.restart, args.seed)
    if not args.prepare_only:
        subprocess.run([sys.executable, str(ROOT/"scripts/runEmbeddedSU2.py"), "--request",
                        str(folder/"engine_request.json")], check=True)
