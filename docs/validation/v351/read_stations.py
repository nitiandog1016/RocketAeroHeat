"""从固定低速算例的真实保存帧分别读取点值、鼻端面平均与底面平均。"""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/simple_launcher"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect(job):
    geometry = json.loads((CASE / "geometry.json").read_bytes())
    condition = json.loads((CASE / "tunnel.json").read_bytes())
    archive = CASE / "wind_tunnel.zip"
    assert digest(archive) == json.loads((CASE / "receipt.json").read_bytes())["sha256"]
    with ZipFile(archive) as bundle:
        ref = json.loads(bundle.read("summary.json"))
        assert ref["accepted"] and ref["condition"] == condition
        assert hashlib.sha256(bundle.read("wall.csv")).hexdigest() == ref["wallSHA256"]
    stl = CASE / geometry["file"]
    assert digest(stl) == geometry["fileSHA256"] == ref["surfaceSHA256"]
    dtype = np.dtype([("normal", "<f4", 3), ("tri", "<f4", (3, 3)), ("attr", "<u2")])
    tri = np.fromfile(stl, dtype=dtype, offset=84)["tri"].astype(float)
    n = len(tri)
    assert n == geometry["faces"]
    area = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)/2
    x = tri[:, :, 0]
    tip = x.max(axis=1) == x.max()
    base = (x == 0).all(axis=1)
    assert tip.sum() == base.sum() == 48
    delta = condition["wallTemperatureK"] - condition["ambientTemperatureK"]
    speed = np.linalg.norm(condition["velocityMps"])
    pressure_scale = condition["densityKgM3"] * speed**2
    # 图 2 的 R* 以底面直径归一化。中心和三圈在本试件只作虚拟测点。
    radial_fractions = np.array([0, .286, .572, .858])
    result = dict(surfaceSHA256=digest(stl), archiveSHA256=digest(archive),
                  tapRadiusOverBaseRadius=radial_fractions.tolist(),
                  noseTipFaceIDs=[int(np.flatnonzero(tip)[0]+1), int(np.flatnonzero(tip)[-1]+1)],
                  noseTipAxialSpanMm=float(np.ptp(x[tip])*1000), noseTipAreaM2=float(area[tip].sum()),
                  baseAreaM2=float(area[base].sum()), grids={})
    for grid in ("coarse", "fine", "outer"):
        flow = job / grid
        receipt = json.loads((flow / "engine_process.json").read_bytes())
        info = json.loads((flow / "condition.json").read_bytes())
        assert receipt["status"] == "terminal" and receipt["exitCode"] == 0
        assert not (receipt["cancelled"] or receipt["timedOut"]) and all(receipt["inputsStable"].values())
        assert digest(flow / "condition.json") == receipt["inputSHA256"][str((flow / "condition.json").resolve())]
        assert info["surfaceSHA256"] == digest(stl) and info["grid"] == grid
        assert (info["speedMps"], info["pRefPa"], info["airK"], info["wallK"]) == (
            speed, condition["ambientPressurePa"], condition["ambientTemperatureK"], condition["wallTemperatureK"])
        seal = ref["runs"][grid]
        frame = flow / "restart.csv"
        assert digest(frame) == receipt["outputSHA256"]["restart.csv"] == seal["restartSHA256"]
        table = job / ("heat_" + grid) / "wall.csv"
        assert digest(table) == seal["wallSHA256"]
        if grid == "outer":
            assert digest(table) == ref["wallSHA256"]
        faces = np.loadtxt(table, delimiter=",", skiprows=1)
        np.testing.assert_array_equal(faces[:, 0], np.arange(1, n+1))
        np.testing.assert_allclose(faces[:, 1], area, rtol=0, atol=1e-15)
        wall = np.asarray(info["wall"])
        count = np.count_nonzero(wall[:, 0] == 0)
        nodes = np.loadtxt(frame, delimiter=",", skiprows=1, max_rows=count)
        np.testing.assert_allclose(nodes[:, 1:3], wall[:count], rtol=0, atol=1e-13)
        np.testing.assert_array_equal(nodes[:, 0], np.arange(count))
        radius = nodes[:, 2] * info["diameterM"]
        assert radius[0] == 0 and np.all(np.diff(radius) > 0)
        np.testing.assert_allclose(radius[-1], np.hypot(tri[base, :, 1], tri[base, :, 2]).max(), atol=1e-12, rtol=0)
        gauge = nodes[:, 3] * pressure_scale * 1000
        # 分段线性压力乘以环带面积权重，精确积分后独立核对每个底面扇片。
        a, b = radius[:-1], radius[1:]
        mean = np.sum((b-a)*((2*a+b)*gauge[:-1] + (a+2*b)*gauge[1:])) / (3*radius[-1]**2)
        displayed = (faces[base, 2] - info["pRefPa"])*1000
        error = np.max(np.abs(displayed-mean))
        assert error < 1e-6
        q = -np.average(faces[tip, 3], weights=area[tip])
        result["grids"][grid] = dict(restartSHA256=seal["restartSHA256"], wallSHA256=digest(table),
            baseRadiusMm=1000*radius[-1], baseRe=condition["densityKgM3"]*speed*2*radius[-1]/condition["dynamicViscosityPaS"],
            tapGaugeMilliPa=np.interp(radial_fractions*radius[-1],radius,gauge).tolist(),
            baseMeanGaugeMilliPa=mean, maxFaceDifferenceMilliPa=error,
            noseTipOutwardMilliWm2=q*1000, noseTipH_Wm2K=q/delta,
            peakOutwardMilliWm2=-1000*faces[:, 3].min(), peakInTip=bool(tip[np.argmin(faces[:, 3])]))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job", type=Path)
    print(json.dumps(inspect(parser.parse_args().job.resolve()), allow_nan=False, indent=2))
