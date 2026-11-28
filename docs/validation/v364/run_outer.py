"""复用同形试件的子午网，在固定低亚声速来流下计算稳态轴对称滑移外流。"""
import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(HERE.parent / "v362"), str(ROOT / "scripts")]
from sphere_case import BIN, SOURCE, TOOL_SHA, foam, initial_field
from runEilmerProcess import sha
from runOpenFoamProcess import run

BASE = Path.home() / ".cache/RocketAeroHeat/v364_validation"
FLIGHT = HERE.parent / "v358/evidence/cases.json"
MODEL = ROOT / "cases/simple_launcher"
GRIDS = {
    "fine": "418bdd26fb24d7dea23f8babc81b14eb062d7b428d89ca7da9d3b0cf2a17c3a9",
    "outer": "9149e7e74d828e3a0b143c44e001a51cc229faa318e16a91740080b6e36ef1ea",
}


def flight():
    archive = json.loads(FLIGHT.read_bytes())
    model = json.loads((MODEL / "geometry.json").read_bytes())
    assert sha(MODEL / model["file"]) == model["fileSHA256"] == archive["surfaceSHA256"]
    case, = [c for c in archive["cases"] if c["mach"] == 0.3]
    assert case["altitudeM"] == 5000 and case["wallGivenK"] == 300
    return case, model


def read_meridian(path):
    """读取已锁定的 SU2 四边形、坐标及边界段，保留原壁段次序。"""
    lines = iter(path.read_text().splitlines())
    def value(key):
        label, value = next(lines).split("=")
        assert label == key
        return value.strip()
    assert value("NDIME") == "2"
    quads = np.array([next(lines).split() for _ in range(int(value("NELEM")))], int)
    assert np.all(quads[:, 0] == 9)
    points = np.array([next(lines).split() for _ in range(int(value("NPOIN")))], float)
    edges = {}
    for _ in range(int(value("NMARK"))):
        key = value("MARKER_TAG")
        rows = np.array([next(lines).split() for _ in range(int(value("MARKER_ELEMS")))], int)
        assert np.all(rows[:, 0] == 3)
        edges[key] = rows[:, 1:]
    assert next(lines, None) is None and set(edges) == {"wall", "axis", "far"}
    return points, quads[:, 1:], edges


def write_wedge(mesh, coords, cells, edges):
    """子午四边形两侧各旋转一度。轴点由两侧共享，轴线不生成零面积面。"""
    xy = coords
    axis = xy[:, 1] == 0
    np.testing.assert_array_equal(np.flatnonzero(axis), np.unique(edges["axis"]))
    assert np.all(xy[:, 1] >= 0)
    angle = np.deg2rad(1)
    lower = np.column_stack((xy[:, 0], xy[:, 1]*np.cos(angle), -xy[:, 1]*np.sin(angle)))
    upper = lower[~axis].copy()
    upper[:, 2] *= -1
    points = np.vstack((lower, upper))
    hi = np.arange(len(xy))
    hi[~axis] = len(xy) + np.arange(len(upper))
    patches = {name: [] for name in ("body", "far", "minus", "plus")}
    names = {tuple(sorted(pair)): ("body" if key == "wall" else key, i)
             for key, pairs in edges.items() for i, pair in enumerate(pairs)}
    opened, internal = {}, []
    for cell, quad in enumerate(cells):
        patches["minus"].append((cell, list(quad[::-1])))
        patches["plus"].append((cell, list(hi[quad])))
        for a, b in zip(quad, np.roll(quad, -1)):
            key = tuple(sorted((a, b)))
            face = list(dict.fromkeys((a, b, hi[b], hi[a])))
            if key in names:
                name, index = names[key]
                if name == "axis":
                    assert len(face) == 2
                else:
                    patches[name].append((index, cell, face))
            elif key in opened:
                owner, first = opened.pop(key)
                internal.append((owner, cell, first))
            else:
                opened[key] = (cell, face)
    assert not opened
    internal.sort(key=lambda v: v[:2])
    faces = [v[2] for v in internal]
    owners = [v[0] for v in internal]
    boundary = []
    for name, rows in patches.items():
        if name in ("body", "far"):
            rows = [(cell, face) for _, cell, face in sorted(rows)]
        kind = {"body": "wall", "far": "patch"}.get(name, "wedge")
        boundary.append(f"{name} {{ type {kind}; nFaces {len(rows)}; startFace {len(faces)}; }}")
        owners.extend(r[0] for r in rows)
        faces.extend(r[1] for r in rows)
    mesh.mkdir(parents=True)
    def write(name, kind, lines):
        rows = list(lines)
        (mesh/name).write_text(foam(name, f"{len(rows)}\n(\n" + "\n".join(rows) + "\n)", kind), encoding="ascii")
    write("points", "vectorField", ("(" + " ".join(f"{v:.17g}" for v in p) + ")" for p in points))
    write("faces", "faceList", (f"{len(f)}(" + " ".join(map(str, f)) + ")" for f in faces))
    write("owner", "labelList", map(str, owners))
    write("neighbour", "labelList", (str(v[1]) for v in internal))
    write("boundary", "polyBoundaryMesh", boundary)


def compute(grid, job):
    """从当前来流冷启动，按相同资源预算监督网格检查和一次稳态投影。"""
    case = job / "case"
    conditions, model = flight()
    path = ROOT / f"outputs/v350_validation/{grid}01/mesh.su2"
    assert sha(path) == GRIDS[grid]
    xy, cells, edges = read_meridian(path)
    xy *= model["maxDiameterM"]
    np.testing.assert_allclose(xy[edges["wall"][-1, 1]], [model["lengthM"], 0], atol=1e-12)
    write_wedge(case / "constant/polyMesh", xy, cells, edges)
    system = case / "system"
    system.mkdir()
    dictionaries = {
        "controlDict": "application potentialFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 1; "
                       "deltaT 1; writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 16;",
        "fvSchemes": "ddtSchemes { default steadyState; } gradSchemes { default Gauss linear; } "
                     "divSchemes { default none; } laplacianSchemes { default Gauss linear corrected; } "
                     "interpolationSchemes { default linear; } snGradSchemes { default corrected; }",
        "fvSolution": "solvers { Phi { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; maxIter 10000; } } "
                      "potentialFlow { nNonOrthogonalCorrectors 4; PhiRefCell 0; PhiRefValue 0; }",
    }
    for key, text in dictionaries.items():
        (system/key).write_text(foam(key, text), encoding="ascii")
    u = f"({-conditions['speedMps']:.17g} 0 0)"
    initial = job / "initial"
    initial.mkdir()
    for field, args in {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", u, "type slip;", f"type fixedValue; value uniform {u};"),
        "Phi": ("volScalarField", "[0 2 -1 0 0 0 0]", "0", "type zeroGradient;", "type zeroGradient;"),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "0", "type zeroGradient;", "type fixedValue; value uniform 0;"),
    }.items():
        (initial/field).write_text(initial_field(field, *args), encoding="ascii")
    shutil.copytree(initial, case / "0")
    assert not (case/"0/phi").exists()
    locks = [HERE/"PLAN.md", Path(__file__), HERE.parent/"v362/sphere_case.py", FLIGHT, path,
             MODEL/"geometry.json", MODEL/model["file"], SOURCE/"potentialFoam.C", SOURCE/"createFields.H",
             *system.iterdir(), *initial.iterdir(), *(case/"constant/polyMesh").iterdir()]
    for tool, args in (("checkMesh", []), ("potentialFoam", ["-initialiseUBCs", "-writePhi"])):
        assert sha(BIN/tool) == TOOL_SHA[tool]
        receipt = run([BIN/tool, *args], case, job/tool, 180, [*locks, BIN/tool])
        log = (job/tool/"stdout.txt").read_text()
        print(f"{grid} {tool} {receipt['elapsedSeconds']:.2f} s", flush=True)
        if not receipt["passed"] or not log.rstrip().endswith("End") or (tool == "checkMesh" and "Mesh OK." not in log):
            raise RuntimeError(f"需要检查 {job/tool}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid", choices=GRIDS)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    job = args.folder.resolve()
    job.relative_to(BASE)
    job.mkdir(parents=True, exist_ok=False)
    compute(args.grid, job)
