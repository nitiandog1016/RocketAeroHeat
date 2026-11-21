"""逐个检查带线热敏珠实际保存的冷态、周期帧和原生终场。"""
import argparse
import csv
import numpy as np

from audit_third_mesh import FINE, read, sha, store, terminal
from thermistorBead import P0, RHO, U
from beadSurfaceResult import thermal_frame


def inspect(path, xyz, wall):
    raw = np.loadtxt(path, delimiter=",", skiprows=1)
    assert raw.shape == (len(xyz), 10) and np.isfinite(raw).all()
    np.testing.assert_array_equal(raw[:, 0], np.arange(len(xyz)))
    np.testing.assert_allclose(raw[:, 1:4], xyz, rtol=0, atol=1e-13)
    theta = raw[:, 9]
    item = dict(file=path.name, sha256=sha(path), thetaRange=[float(theta.min()), float(theta.max())],
                outside=int(np.count_nonzero((theta < -1e-12) | (theta > 1+1e-12))),
                minPressurePa=float(P0 + RHO*U**2*raw[:, 4].min()),
                maxWallSpeedMps=float(U*np.linalg.norm(raw[wall, 5:8], axis=1).max()),
                maxWallThetaError=float(abs(theta[wall]-1).max()))
    return item


def check(live):
    reg = read(FINE / "inputs.json")
    assert reg["allTargetsPassed"] and all(sha(FINE/k) == v for k,v in reg["inputs"].items())
    with (FINE / "history.csv").open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    idx = rows[0].index("Inner_Iter")
    last = int([r for r in rows[1:] if len(r) == len(rows[0])][-1][idx])
    rec = None if live else terminal(FINE)
    with np.load(FINE / "nodes.npz") as data:
        xyz, wall = data["xyz"], np.unique(data["wallMap"])
    frames = [dict(step=None, **inspect(FINE/"initial_flow.csv", xyz, wall))]
    for path in sorted(FINE.glob("restart_*.csv")):
        step = int(path.stem.split("_")[-1])
        if live and step >= last:
            continue
        surf = FINE / f"surface_{step:06d}.vtk"
        item = dict(step=step, **inspect(path, xyz, wall), surfaceSHA256=sha(surf))
        if rec:
            assert item["sha256"] == rec["outputSHA256"][path.name]
            assert item["surfaceSHA256"] == rec["outputSHA256"][surf.name]
        frames.append(item)
    if rec:
        assert [item["step"] for item in frames] == [None, *range(500, last, 500), last]
        assert frames[-1]["sha256"] == sha(FINE/"restart.csv") and frames[-1]["surfaceSHA256"] == sha(FINE/"surface.vtk")
    good = all(x["outside"] == 0 and x["minPressurePa"] > 0 and x["maxWallSpeedMps"] < 1e-10*U and x["maxWallThetaError"] == 0 for x in frames)
    result = dict(frames=frames, lastHistoryIteration=last, fieldChecksPassed=good)
    if rec:
        result["residualWindow"] = thermal_frame(FINE)[1]
        store(FINE.parent/"checkpoints.json", result)
    print(last, [(x["step"], x["thetaRange"], x["outside"]) for x in frames])
    assert good and (not rec or max(result["residualWindow"].values()) < -8)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    check(parser.parse_args().live)
