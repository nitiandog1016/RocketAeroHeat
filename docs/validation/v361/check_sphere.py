"""在标准滑移球上直接核对不可压无黏引擎的压力、速度和残差。"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / "outputs/v361_validation"
sys.path.insert(0, str(ROOT / "scripts"))
from cm011Tunnel import volume, write_volume_mesh
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe


def grid(level):
    theta = np.linspace(np.pi, 0, 64 * level + 1)
    wall = 0.5 * np.column_stack((np.cos(theta), np.sin(theta)))
    wall[[0, -1], 1] = 0
    return volume(wall, 36 * level, 14, level)


def run(folder, level):
    """只启动固定球体对照；既有运行目录不能覆盖。"""
    exe = ROOT / "tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
    if sha(exe) != "01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7":
        raise ValueError("SU2 版本指纹改变")
    xy, cells, markers = grid(level)
    folder.mkdir(exist_ok=False, parents=True)
    write_volume_mesh(folder / "mesh.su2", xy, cells, markers)
    shutil.copyfile(HERE / "sphere.cfg", folder / "sphere.cfg")
    np.savez_compressed(folder / "initial.npz", xy=xy, pressure=np.zeros(len(xy)),
                        velocity=np.tile([-1., 0.], (len(xy), 1)))
    inputs = [exe, Path(__file__).resolve(), HERE / "PLAN.md", HERE / "sphere.cfg",
              ROOT / "scripts/cm011Tunnel.py", folder / "mesh.su2", folder / "initial.npz", folder / "sphere.cfg"]
    if level == 2:
        inputs.append(HERE / "GRID_CHECK.md")
    receipt = observe([str(exe), "sphere.cfg", "--threads", "1"], folder,
                      {str(p): sha(p) for p in inputs}, 360 if level == 1 else 180, 80 * 1024**2)
    check_process(folder, receipt)


def check_process(folder, rc):
    if (rc["status"] != "terminal" or rc["debuggerExitCode"] != 0 or rc["stopReason"]
            or rc["supervisorError"] or rc["survivingProcessIDs"] or not all(rc["inputsStable"].values())
            or "Exit Success (SU2_CFD)" not in (folder / "debugger.log").read_text(errors="replace")):
        raise ValueError("运行或保存未正常完成，请查原生回执")


def columns(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        labels = [k.strip() for k in next(csv.reader(stream, skipinitialspace=True))]
        values = np.loadtxt(stream, delimiter=",", ndmin=2)
    if len(labels) != values.shape[1]:
        raise ValueError("数据列与字段不符")
    return dict(zip(labels, values.T))


def review(folder):
    """用真实末帧与无界域球体解析解比较，定压和无滑移均不作为球壁观测值。"""
    rc = json.loads((folder / "debugger_process.json").read_bytes())
    check_process(folder, rc)
    hist = columns(folder / "history.csv")
    iterations = hist["Inner_Iter"]
    np.testing.assert_array_equal(iterations, np.arange(len(iterations)))
    assert len(iterations) >= 20
    last = int(iterations[-1])
    frame = folder / f"restart_{last:06}.csv"
    surface = folder / f"surface_{last:06}.csv"
    data, wall_data = columns(frame), columns(surface)
    n_wall = len(wall_data["PointID"])
    if n_wall not in (65, 129):
        raise ValueError("原生球壁点数不属于本次两档网格")
    xy, _, _ = grid((n_wall-1)//64)
    order = np.argsort(data["PointID"])
    np.testing.assert_array_equal(data["PointID"][order], np.arange(len(xy)))
    np.testing.assert_allclose(np.column_stack((data["x"], data["y"]))[order], xy, rtol=0, atol=1e-9)
    wall_ids = np.argsort(wall_data["PointID"])
    np.testing.assert_array_equal(wall_data["PointID"][wall_ids], np.arange(n_wall))
    velocity = np.column_stack((data["Velocity_x"], data["Velocity_y"]))[order]
    press = data["Pressure"][order]
    for key in ("Pressure", "Velocity_x", "Velocity_y"):
        np.testing.assert_allclose(data[key][order[:n_wall]], wall_data[key][wall_ids], rtol=1e-9, atol=1e-10)
    assert np.all(np.isfinite(velocity)) and np.all(np.isfinite(press))
    normal = xy[:n_wall] / 0.5
    tangent = np.column_stack((-normal[:, 1], normal[:, 0]))
    speed_t = np.sum(velocity[:n_wall] * tangent, axis=1)
    reference_t = 1.5 * normal[:, 1]
    reference_cp = 1 - reference_t**2
    cp = 2 * press[:n_wall]
    cp_error = cp - reference_cp
    # 球环的面积正比于相邻极角的 cos 之差；每个节点占两侧各半个球环。
    band = np.diff(normal[:, 0])
    weight = np.r_[band, 0] / 2 + np.r_[0, band] / 2
    rms = lambda error: float(np.sqrt(np.average(error**2, weights=weight)))
    residuals = {field: float(np.max(hist[field][-20:])) for field in ("rms[P]", "rms[U]", "rms[V]")}
    # 唯一外界的节点速度面积求积，以来流量除球体大圆面积归一化。
    v_far_normal = np.sum(velocity[-n_wall:] * normal, axis=1)
    far_net = float(2 * 14**2 * np.dot(weight, v_far_normal) / 0.5**2)
    tests = {"pressureRMS": rms(cp_error) < 0.04, "pressureMaximum": max(abs(cp_error)) < 0.10,
             "tangentialRMS": rms(speed_t-reference_t) < 0.03, "residuals": max(residuals.values()) < -8}
    return dict(passed=bool(all(tests.values())), tests={k: bool(v) for k, v in tests.items()},
                finalIteration=last, iterations=len(iterations), seconds=rc["elapsedSeconds"],
                peakRSSBytes=rc["peakFamilyRSSBytes"], pressureCoefficientRMS=rms(cp_error),
                pressureCoefficientMax=float(max(abs(cp_error))), tangentialVelocityRMS=rms(speed_t-reference_t),
                worstPressureAngleDeg=float(np.rad2deg(np.arccos(normal[np.argmax(abs(cp_error)), 0]))),
                frontPressureCoefficient=float(cp[-1]), rearPressureCoefficient=float(cp[0]),
                equatorPressureCoefficient=float(cp[n_wall//2]),
                maxWallNormalSpeed=float(max(abs(np.sum(velocity[:n_wall]*normal, axis=1)))),
                farNetFluxOverFreestreamDiskFlux=far_net, last20LogRMSMax=residuals,
                files={p.name: sha(p) for p in (frame, surface, folder/"history.csv", folder/"mesh.su2",
                                                 folder/"sphere.cfg", folder/"debugger.log", folder/"debugger_process.json")})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "read"))
    parser.add_argument("folder", type=Path)
    parser.add_argument("--level", type=int, choices=(1, 2), default=1, help="仅在新算时选网格")
    args = parser.parse_args()
    folder = args.folder.resolve()
    folder.relative_to(BASE)
    if args.action == "run":
        run(folder, args.level)
    else:
        print(json.dumps(review(folder), indent=2, allow_nan=False))
