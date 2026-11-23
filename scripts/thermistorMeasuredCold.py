"""B 型热敏珠风洞测点：仅由 96 档几何和已确认的 12 K 边界构造初场、推进及审核保存帧。"""
import csv
import os
from pathlib import Path
import uuid

import numpy as np

from runThermistorLeadCase import configuration
from thermistorBead import ENGINE, ENGINE_SHA, ROOT, read, sha, solve, write_json
from thermistorReaction import terminal


GRID_KEY = ROOT / "cases/thermistor_b/sidewindAngular96.json"
MESH_FILES = {"mesh.su2": "meshSHA256", "nodes.npz": "nodeMapSHA256", "body.json": "bodySHA256"}
RESIDUALS = ["rms[" + k + "]" for k in ("P", "U", "V", "W", "h", "rho*Y_0")]
HEADER = '"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"'
STEPS = 320


def grid(prepared):
    """允许工作台重建同一网，身份由实际内容决定。"""
    key = read(GRID_KEY)
    for name, field in dict(MESH_FILES, **{"case.json": "caseSHA256", "mesh_audit.json": "mappingAuditSHA256"}).items():
        assert sha(prepared / name) == key[field], name
    spec = read(prepared / "case.json")
    assert spec["caseKind"] == "thermistor-B-sidewind-two-leads" and spec["velocityMps"] == [0, -3, 0]
    assert spec["airPa"] == 101325 and spec["airK"] == 293.15 and spec["wallK"] == 299.15
    spec.update(wallK=305.15, scalarReferenceDeltaK=24)
    return spec


def config(steps):
    """从项目配置生成侧向来流、强标量及弱焓边界，不访问旧运行目录。"""
    cfg = configuration()
    changes = {"MARKER_SYM= (sym0, sym1)\n": "", "INC_VELOCITY_INIT= (-1, 0, 0)": "INC_VELOCITY_INIT= (0, -1, 0)",
               "(wall, VALUE, 1.0)": "(wall, VALUE, 0.5)", "MARKER_ISOTHERMAL= (wall, 2)": "MARKER_ISOTHERMAL= (wall, 1.5)",
               "OUTPUT_WRT_FREQ= (500, 500)": "OUTPUT_WRT_FREQ= (100, 100)", "ITER= 7000\n": f"ITER= {steps}\n"}
    for old, new in changes.items():
        assert cfg.count(old) == 1, old
        cfg = cfg.replace(old, new)
    return cfg


def frame(path, xyz, wall, cold=False):
    """核对完整节点顺序、强温界、无滑移与绝压，冷态另核对每个自由度。"""
    with path.open(encoding="ascii") as stream:
        assert stream.readline().strip() == HEADER
        state = np.loadtxt(stream, delimiter=",")
    assert state.shape == (len(xyz), 10) and np.isfinite(state).all()
    np.testing.assert_array_equal(state[:, 0], np.arange(len(xyz)))
    np.testing.assert_allclose(state[:, 1:4], xyz, atol=1e-12, rtol=0)
    theta = state[:, 9]
    pressure = 101325 + 1.18 * 9 * state[:, 4]
    assert pressure.min() > 0 and theta.min() >= -1e-12 and theta.max() <= .5 + 1e-12
    assert np.all(theta[wall] == .5) and np.max(abs(state[wall, 5:8])) < 1e-10
    if cold:
        interior = np.ones(len(xyz), dtype=bool)
        interior[wall] = False
        assert np.all(state[:, [4, 5, 7]] == 0) and np.all(state[wall, 6] == 0)
        assert np.all(state[interior, 6] == -1) and np.all(theta[interior] == 0)
        assert np.all(state[interior, 8] == 1 - 298.15) and np.all(state[wall, 8] == 1.5 - 298.15)
    return dict(sha256=sha(path), pressurePa=[float(pressure.min()), float(pressure.max())],
                thetaRange=[float(theta.min()), float(theta.max())])


def seed(prepared, folder, spec):
    """保存生成冷态所需的实际网与初值，供独立接续。"""
    folder.mkdir()
    for name in MESH_FILES:
        os.link(prepared / name, folder / name)
    write_json(folder / "case.json", spec)
    with np.load(folder / "nodes.npz") as mapping:
        xyz, wall = mapping["xyz"], np.unique(mapping["wallMap"])
    state = np.zeros((len(xyz), 10))
    state[:, 0] = np.arange(len(xyz))
    state[:, 1:4] = xyz
    state[:, 6], state[wall, 6] = -1, 0
    state[:, 8], state[wall, 8], state[wall, 9] = 1 - 298.15, 1.5 - 298.15, .5
    path = folder / "initial_flow.csv"
    np.savetxt(path, state, fmt=["%d"] + ["%.17g"] * 9, delimiter=",", comments="", header=HEADER)
    checked = frame(path, xyz, wall, cold=True)
    write_json(folder / "cold.json", dict(initial=checked, geometry=read(GRID_KEY)["meshSHA256"], nodes=len(xyz),
               wallNodes=len(wall), caseSHA256=sha(folder / "case.json"), initialSource="geometry-and-boundary"))


def advance(run, dest, initial, steps):
    """每段使用独立目录，保存源绑定、原生监督回执及周期场。"""
    assert sha(ENGINE) == ENGINE_SHA and 1 <= steps <= STEPS
    dest.mkdir()
    for name in (*MESH_FILES, "case.json"):
        os.link(run / "cold" / name, dest / name)
    os.link(initial, dest / "initial_flow.csv")
    (dest / "case.cfg").write_text(config(steps), encoding="ascii")
    files = [Path(__file__).resolve(), ROOT / "scripts/runThermistorLeadCase.py", ROOT / "scripts/thermistorBead.py",
             GRID_KEY, run / "job.json", run / "cold/cold.json", *(dest / n for n in (*MESH_FILES, "case.json"))]
    write_json(dest / "engine_request.json", dict(workDirectory=str(dest), configFile="case.cfg", meshFile="mesh.su2",
               threads=4, maximumSeconds=720, jobToken=str(uuid.uuid4()), sourceFiles=[p.relative_to(ROOT).as_posix() for p in files]))
    solve(dest, memory_bytes=4 * 1024**3, output_bytes=1536 * 1024**2)


def audit(run, dest, initial, previous, steps):
    """从独立原生记录重建可接续步数，验证所有计划保存帧与末窗。"""
    receipt = terminal(dest)
    assert receipt["request"]["workDirectory"] == str(dest) and receipt["inputSHA256"][str(run / "job.json")] == sha(run / "job.json")
    assert receipt["inputSHA256"][str(dest / "initial_flow.csv")] == sha(initial)
    assert (dest / "case.cfg").read_text("ascii") == config(steps)
    for name in (*MESH_FILES, "case.json"):
        assert sha(dest / name) == receipt["inputSHA256"][str(dest / name)] == sha(run / "cold" / name)
    with (dest / "history.csv").open(encoding="ascii") as stream:
        header = [k.strip() for k in next(csv.reader(stream, skipinitialspace=True))]
        history = np.loadtxt(stream, delimiter=",", ndmin=2)
    assert 1 <= len(history) <= steps and np.isfinite(history).all()
    np.testing.assert_array_equal(history[:, header.index("Inner_Iter")], np.arange(len(history)))
    last = len(history) - 1
    assert last == receipt["progress"]["Inner_Iter"]
    window = {k: float(history[-20:, header.index(k)].max()) for k in RESIDUALS}
    assert len(history) == steps or (len(history) >= 20 and max(window.values()) < -8)
    with np.load(dest / "nodes.npz") as mapping:
        xyz, wall = mapping["xyz"], np.unique(mapping["wallMap"])
    frames = {}
    for i in sorted(set(range(100, last + 1, 100)) | {last}):
        name = f"restart_{i:06d}.csv"
        info = frame(dest / name, xyz, wall)
        surface = f"surface_{i:06d}.vtk"
        info["surfaceSHA256"] = sha(dest / surface)
        assert info["sha256"] == receipt["outputSHA256"][name] and info["surfaceSHA256"] == receipt["outputSHA256"][surface]
        frames[name] = info
    assert frames[name]["sha256"] == receipt["outputSHA256"]["restart.csv"]
    result = dict(steps=len(history), cumulativeSteps=previous + len(history), initialSHA256=sha(initial),
                  lastFrame=frames[name]["sha256"], frames=frames, last20Residual=window,
                  allResidualsBelowMinus8=len(history) >= 20 and max(window.values()) < -8, seconds=receipt["elapsedSeconds"])
    path = dest / "field_audit.json"
    if path.exists():
        assert read(path) == result
    else:
        write_json(path, result)
    return result
