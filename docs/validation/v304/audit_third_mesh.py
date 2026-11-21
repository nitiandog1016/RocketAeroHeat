"""从成品 SU2 网格和原面映射判断第三档中段网能否进入热敏珠冷态求解。"""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v300"))
from mesh_audit import ROOT, WIRED, NEW as MEDIUM, cells, crossings, extent, meridian, read, sha, terminal
from thermistorOuterDomain import load_mesh

FINE = ROOT / "outputs/v304_validation/wired01"


def store(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf8") as stream:
        json.dump(data, stream, indent=2, ensure_ascii=False, allow_nan=False)


def check():
    previous = read(HERE.parent / "v300/evidence/inputs.json")
    for folder, reference in ((MEDIUM, previous), (WIRED, read(HERE.parent / "v297/evidence/fine_inputs.json"))):
        assert all(sha(folder / name) == value for name, value in reference["inputs"].items())
        terminal(folder)
    prep = read(FINE / "debugger_process.json")
    assert prep["status"] == "terminal" and prep["debuggerExitCode"] == 0 and all(prep["inputsStable"].values())
    assert not (prep["survivingProcessIDs"] or prep["stopReason"] or prep["supervisorError"])
    locked = {name: sha(FINE / name) for name in previous["inputs"]}
    for name in ("body.json", "case.cfg", "case.json"):
        assert locked[name] == previous["inputs"][name]
    # 壁面、首层和各近壁层接口按真实坐标逐点比较；外边界按有序顶点坐标核对。
    geom = read(FINE / "geometry.json")
    with np.load(MEDIUM / "nodes.npz") as old, np.load(FINE / "nodes.npz") as new:
        for field in ("wallMap", "firstMap"):
            np.testing.assert_array_equal(old["xyz"][old[field]], new["xyz"][new[field]])
        np.testing.assert_array_equal(old["area"], new["area"])
        assert old["ids"].shape[0] == new["ids"].shape[0] == 5
        layers = 8 * (geom["meridianIntervals"] + 1)
        np.testing.assert_array_equal(old["xyz"][old["ids"][:, :layers]], new["xyz"][new["ids"][:, :layers]])
        cold = np.loadtxt(FINE / "initial_flow.csv", delimiter=",", skiprows=1)
        np.testing.assert_array_equal(cold[:, 0], np.arange(len(new["xyz"])))
        np.testing.assert_array_equal(cold[:, 1:4], new["xyz"])
        theta = np.zeros(len(cold))
        theta[np.unique(new["wallMap"])] = 1
        np.testing.assert_array_equal(cold[:, 9], theta)
        np.testing.assert_array_equal(cold[:, 8], 1 - 298.15 + theta)
        np.testing.assert_array_equal(cold[:, 5], theta - 1)
        assert not np.any(cold[:, [4, 6, 7]])
    def far_faces(folder):
        _, xyz, marks = load_mesh(folder / "mesh.su2")
        return {tuple(map(tuple, xyz[list(face)])) for face in marks["far"]}
    far = far_faces(MEDIUM)
    assert len(far) == 128 and far == far_faces(FINE)
    mid_cells, fine_cells = (set.union(*cells(folder)) for folder in (MEDIUM, FINE))
    wall = meridian()[0]
    removed = extent(mid_cells - fine_cells, wall)
    added = extent(fine_cells - mid_cells, wall)
    assert max(map(abs, removed["xOverA"] + added["xOverA"])) < previous["jointInnerAbsXOverA"]
    for change in (added, removed):
        assert change["vertexWallDistanceD"][0] > .014 and change["vertexWallDistanceD"][1] < .16
    readings = [crossings(folder) for folder in (WIRED, MEDIUM, FINE)]
    assert readings[0] == previous["coarse"] and readings[1] == previous["refined"]
    target = [m["degFromUpstream"] for m, f in zip(*readings[1:]) if
              m["degFromUpstream"] in (60, 85, 120) and f["maxCrossingIntervalOverD"] < m["maxCrossingIntervalOverD"]]
    request = read(FINE / "engine_request.json")
    assert request["threads"] == 1 and request["maximumSeconds"] == 720
    result = dict(inputs=locked, geometry=geom, sourceSHA256={p:sha(ROOT/p) for p in request["sourceFiles"]},
                  readingsByGrid=readings, reducedAtDegrees=target, allTargetsPassed=target == [60,85,120],
                  identicalCells=len(mid_cells & fine_cells), removed=removed, added=added,
                  preparationSeconds=prep["elapsedSeconds"], unchangedFarFaces=len(far))
    store(FINE / "inputs.json", result)
    print("Three-grid maximum intervals at 60, 85, 120 degrees:", [[v["maxCrossingIntervalOverD"] for v in row[1:4]] for row in readings])
    assert result["allTargetsPassed"], "目标法向仍有间距未缩小，停止在几何核对"


if __name__ == "__main__":
    check()
