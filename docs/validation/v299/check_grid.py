"""启动前从实际 SU2 文件核对原面、首层、变化的单元及三条目标法向。"""
from pathlib import Path
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/"v298"))
from check_mesh import cells
from mesh_readings import crossings
from thermistorBead import A, D, read, sha
from thermistorLeadMesh import meridian
from thermistorMidbodyMesh import wall_distance
from thermistorReaction import terminal

OLD = ROOT/"outputs/v298_validation/smooth01"
NEW = ROOT/"outputs/v299_validation/smooth01"


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf8") as out:
        json.dump(data, out, indent=2, ensure_ascii=False, allow_nan=False)


def extent(changed, wall):
    coords = np.asarray(sorted({v for _, cell in changed for v in cell}))
    meridian_points = np.column_stack((coords[:, 0], np.hypot(coords[:, 1], coords[:, 2])))
    distance = wall_distance(meridian_points, wall)
    axial = coords[:, 0]*D/A
    return dict(cells=len(changed), xOverA=[float(axial.min()), float(axial.max())],
                vertexWallDistanceD=[float(distance.min()), float(distance.max())])


def run():
    base = read(HERE.parent/"v298/evidence/inputs.json")
    assert all(sha(OLD/name) == value for name, value in base["inputs"].items())
    terminal(OLD)
    prepare = read(NEW/"debugger_process.json")
    assert prepare["status"] == "terminal" and prepare["debuggerExitCode"] == 0
    assert not prepare["stopReason"] and not prepare["supervisorError"] and all(prepare["inputsStable"].values())
    inputs = {name: sha(NEW/name) for name in base["inputs"]}
    assert all(inputs[name] == base["inputs"][name] for name in ("case.cfg", "case.json", "body.json"))
    request = read(NEW/"engine_request.json")
    assert request["threads"] == 1 and request["maximumSeconds"] == 720
    with np.load(OLD/"nodes.npz") as a, np.load(NEW/"nodes.npz") as b:
        for k in ("wallMap", "firstMap"):
            np.testing.assert_array_equal(a["xyz"][a[k]], b["xyz"][b[k]])
        np.testing.assert_array_equal(a["area"], b["area"])
    old = set.union(*cells(OLD))
    new = set.union(*cells(NEW))
    wall = meridian(True)[0]
    removed, added = extent(old-new, wall), extent(new-old, wall)
    for changed in (removed, added):
        assert changed["vertexWallDistanceD"][0] > .014 and changed["vertexWallDistanceD"][1] < .25
        assert max(map(abs, changed["xOverA"])) < 1
    coarse, refined = crossings(OLD), crossings(NEW)
    assert coarse == read(HERE.parent/"v298/evidence/spatial.json")["newMixed"]
    for a, b in zip(coarse, refined):
        if a["degFromUpstream"] in (60, 85, 120):
            assert b["maxCrossingIntervalOverD"] < a["maxCrossingIntervalOverD"]
    save(HERE/"evidence/inputs.json", dict(baseSHA=base["inputs"]["mesh.su2"], inputs=inputs,
         geometry=read(NEW/"geometry.json"), commonCells=len(old & new), oldCells=len(old), newCells=len(new),
         removed=removed, added=added, coarse=coarse, refined=refined,
         sources={p:sha(ROOT/p) for p in request["sourceFiles"]}, prepareSeconds=prepare["elapsedSeconds"]))


if __name__ == "__main__":
    run()
