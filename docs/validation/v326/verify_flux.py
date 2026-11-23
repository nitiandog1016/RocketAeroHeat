"""从原生历史、全节点更新和面通量重建 Short Test 65 的拒绝证据。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "cases/short_sphere"
BASE = ROOT / "outputs/v326_validation"
OLD = ROOT / "outputs/v325_validation"
sys.path.insert(0, str(ROOT / "scripts"))
from shortSphere65 import CASE, configuration, read, sha


def table(folder, name, header=False):
    data = np.loadtxt(folder / name, delimiter=",", skiprows=int(header), ndmin=2)
    assert np.isfinite(data).all(), name
    return data


def terminal(folder):
    record = read(folder / "debugger_process.json")
    assert record["status"] == "terminal" and record["debuggerExitCode"] == 0
    assert all(record["inputsStable"].values()) and record["observedProcessIDs"]
    assert not (record["survivingProcessIDs"] or record["stopReason"] or record["supervisorError"])
    assert "Exit Success (SU2_CFD)" in (folder / "debugger.log").read_text(errors="replace")
    return record


def sorted_by_id(rows):
    rows = rows[np.argsort(rows[:, 0])]
    ids = rows[:, 0].astype(int)
    np.testing.assert_array_equal(ids, rows[:, 0])
    assert len(np.unique(ids)) == len(ids)
    return rows


def reconstruct(detail, xyz, wet, ceiling, pe):
    trace = table(detail, "extrema.csv")
    np.testing.assert_array_equal(trace[:, 0], np.arange(200))
    step = int(np.flatnonzero((trace[:, 4:6] > 0).any(axis=1))[0])
    states = sorted_by_id(table(detail, "first_update.csv", header=True))
    np.testing.assert_array_equal(states[:, 0], np.arange(len(xyz)))
    before, after, spatial, delta, relax = states[:, 1:].T
    below = np.flatnonzero(after < -1e-12)
    assert before.min() >= -1e-12 and before.max() <= ceiling + 1e-12
    assert len(below) == trace[step, 4] and (after > ceiling + 1e-12).sum() == trace[step, 5] == 0
    point = int(np.argmin(after))
    np.testing.assert_array_equal(trace[step, [1, 2, 6, 7, 8, 9]],
                                  [point, after[point], before[point], spatial[point], delta[point], relax[point]])
    np.testing.assert_allclose(before[below] + relax[below] * delta[below], after[below], atol=1e-26, rtol=0)
    assert not np.intersect1d(below, wet).size
    first = dict(iteration=step, update=step + 1, belowIDs=below.tolist(), worstID=point,
                 xyzD=xyz[point].tolist(), nearestWallNodeD=float(np.linalg.norm(xyz[wet]-xyz[point], axis=1).min()),
                 before=before[point], after=after[point], delta=delta[point], relaxation=relax[point])

    # 按流动与标量各自记录的边号、节点号对齐，正号为目标节点的空间残差。
    edges = sorted_by_id(table(detail, "scalar_edges.csv"))
    flow = sorted_by_id(table(detail, "flow_edges.csv"))
    np.testing.assert_array_equal(flow, edges[:, [0, 1, 2, 6]])
    i, j = edges[:, 1:3].astype(int).T
    assert np.all((i == point) ^ (j == point))
    sign = np.where(i == point, 1., -1.)
    mass, left, right, adv_i, adv_j, diff_i, diff_j, total_i, total_j = edges[:, 6:].T
    np.testing.assert_array_equal(left, before[i])
    np.testing.assert_array_equal(right, before[j])
    np.testing.assert_array_equal(adv_i, mass * np.where(mass >= 0, left, right))
    np.testing.assert_array_equal(edges[:, [9, 11, 13]], -edges[:, [10, 12, 14]])
    np.testing.assert_allclose(adv_i + diff_i, total_i, atol=1e-28, rtol=0)
    adv, diff = sign * adv_i, sign * diff_i
    divergence = -sign * mass * before[point]
    dx = xyz[j] - xyz[i]
    projection = (edges[:, 3:6] * dx).sum(axis=1) / (dx * dx).sum(axis=1)
    two_point = -sign * (right-left) * projection / pe
    correction = diff - two_point
    # 冻结当前场，逐边比较修正项与两点项的幅值。这个估算不更新未知量。
    limited = np.sign(correction) * np.minimum(abs(correction), abs(two_point))
    assembled = (adv + divergence + diff).sum()
    row, = table(detail, "row.csv")
    iteration, original, rhs, volume, pseudo_dt = row
    assert iteration == step and original == spatial[point] == -rhs
    np.testing.assert_allclose(assembled, original, atol=1e-27, rtol=0)
    matrix = sorted_by_id(table(detail, "matrix.csv"))
    ids = matrix[:, 0].astype(int)
    np.testing.assert_array_equal(ids, np.unique(np.r_[i, j]))
    np.testing.assert_array_equal(matrix[:, 2], delta[ids])
    products = matrix[:, 1] * matrix[:, 2]
    row_result = dict(edges=len(edges), advection=adv.sum(), massDivergence=divergence.sum(),
                      diffusion=diff.sum(), twoPoint=two_point.sum(), gradientCorrection=correction.sum(),
                      spatial=original, rhs=rhs, diagonal=matrix[ids == point, 1].item(),
                      offDiagonalProduct=products[ids != point].sum(), matrixTimesDelta=products.sum(),
                      linearResidual=rhs-products.sum(), volumeD3=volume, pseudoDt=pseudo_dt,
                      limitedEdges=int((abs(correction) > abs(two_point)).sum()),
                      frozenLimitedSpatial=(adv+divergence+two_point+limited).sum())

    # 对两种远场数值通量分别积分；量纲与迭代阶段同原生采样。
    far = sorted_by_id(table(detail, "scalar_far.csv"))
    flow_far = sorted_by_id(table(detail, "flow_far.csv"))
    np.testing.assert_array_equal(far[:, 0], flow_far[:, 0])
    mass = far[:, 1]
    np.testing.assert_allclose(mass, (far[:, 2:5]*far[:, 5:8]).sum(axis=1)*far[:, 8], atol=2e-15, rtol=1e-14)
    np.testing.assert_allclose(far[:, 11], mass*np.where(mass >= 0, far[:, 9], far[:, 10]), atol=1e-25, rtol=1e-14)
    boundary = dict(nodes=len(far), flowSignedSum=flow_far[:, 1].sum(), scalarSignedSum=mass.sum(),
                    maxAbsoluteDifference=float(abs(mass-flow_far[:, 1]).max()))
    return first, row_result, boundary


def audit():
    prep = OLD / "prep02/prepared"
    baseline, new = OLD / "run01/steady01", BASE / "run01/steady01"
    scout, detail = BASE / "scout01", BASE / "detail01"
    grid = read(DATA / "supportedGrid65.json")
    spec = read(prep / "prepared.json")
    audit_new = read(new.parent / "qualified.json")
    ref = read(baseline.parent / "qualified.json")
    assert audit_new == {**ref, "speciesConvection": "BOUNDED_SCALAR"}
    assert sha(prep / "prepared.json") == ref["prepared"]
    assert sha(CASE) == ref["caseSHA256"] == grid["caseSHA256"]
    assert grid["hybrid"]["meshSHA256"] == spec["files"]["mesh.su2"]
    assert sha(prep / "mapping.npz") == spec["files"]["mapping.npz"]
    with np.load(prep / "mapping.npz") as mapping:
        xyz, wet, theta = (mapping[k] for k in ("xyz", "wet", "theta"))
    ceiling = float(theta.max())
    assert sorted(p.name for p in new.glob("restart_*.csv")) == ["restart_000100.csv", "restart_000199.csv"]
    for folder in (baseline, new, scout, detail):
        receipt = terminal(folder)
        for local, key in (("mesh.su2", "mesh.su2"), ("initial_flow.csv", "initial.csv"), ("wire_boundary.dat", "profile.dat")):
            assert sha(folder/local) == receipt["inputSHA256"][str(folder/local)] == spec["files"][key]
        assert sha(folder/"restart.csv") == sha(folder/"restart_000199.csv")
        name = "SCALAR_UPWIND" if folder == baseline else "BOUNDED_SCALAR"
        cfg = (folder/"case.cfg").read_text("ascii")
        assert cfg == configuration(spec["condition"], 200, False, name)
        assert sha(folder/"case.cfg") == receipt["inputSHA256"][str(folder/"case.cfg")]
    old_cfg = (baseline/"case.cfg").read_bytes()
    assert old_cfg.count(b"SCALAR_UPWIND") == 1
    assert old_cfg.replace(b"SCALAR_UPWIND", b"BOUNDED_SCALAR") == (new/"case.cfg").read_bytes()
    smooth = ROOT/"outputs/v321_validation/run01"
    smooth_cfg = configuration(read(smooth/"prepared/prepared.json")["condition"], 400, False)
    assert (smooth/"steady01/case.cfg").read_text("ascii") == smooth_cfg
    history = table(new, "history.csv", True)
    np.testing.assert_array_equal(history[:, 2], np.arange(200))
    window = history[-20:, 3:9].max(axis=0)
    assert max(window) >= -8 and not (new.parent/"endpoints.json").exists()
    assert read(new.parent/"failure.json")["stage"] == "steady01"
    records = [read(folder/"replay.json") for folder in (scout, detail)]
    for folder, record in zip((scout, detail), records):
        assert all(sha(folder/key) == digest for key, digest in record["files"].items())
        assert record["build"] == sha(BASE/"build01/probe.exe")
        np.testing.assert_array_equal(table(folder, "history.csv", True)[:, 2:9], history[:, 2:9])
    for key in ("extrema.csv", "first_update.csv"):
        assert records[0]["files"][key] == records[1]["files"][key]
    first, row, far = reconstruct(detail, xyz, wet, ceiling, spec["condition"]["pe"])
    trace = table(detail, "extrema.csv")
    frames = []
    for ref_frame in grid["hybrid"]["stages"][-1]["frames"]:
        index = ref_frame["index"]
        file = f"restart_{index:06d}.csv"
        arr = table(new, file, True)
        old = table(baseline, file, True)
        assert sha(baseline/file) == ref_frame["sha256"]
        np.testing.assert_array_equal(arr[:, :9], old[:, :9])
        np.testing.assert_array_equal(arr[:, 0], np.arange(len(xyz)))
        np.testing.assert_allclose(arr[:, 1:4], xyz, atol=1e-12, rtol=0)
        values = arr[:, 9]
        np.testing.assert_allclose(values[wet], theta, atol=2e-14, rtol=0)
        below, above = int((values < -1e-12).sum()), int((values > ceiling+1e-12).sum())
        assert [below, above] == trace[index, 4:6].tolist() and below > 0
        for folder, record in zip((scout, detail), records):
            delta = float(abs(table(folder, file, True)-arr).max())
            assert delta == record["frameMaxDifference"][file] and delta < 1e-13
        frames.append(dict(index=index, below=below, above=above, minScalar=values.min(), maxScalar=values.max(),
                           sha256=sha(new/file), baselineSHA256=ref_frame["sha256"],
                           baselineBelow=ref_frame["below"], baselineMinimum=ref_frame["minScalar"]))
    assert records[0]["firstIteration"] == records[1]["firstIteration"] == first["iteration"]
    assert records[0]["worstGlobalID"] == records[1]["worstGlobalID"] == first["worstID"]
    return dict(caseSHA256=sha(CASE), meshSHA256=spec["files"]["mesh.su2"], coldSHA256=spec["files"]["initial.csv"],
                wallSHA256=spec["files"]["profile.dat"], preparedSHA256=ref["prepared"], configSHA256=sha(new/"case.cfg"),
                method="BOUNDED_SCALAR", status="rejected", iterations=len(history), lowerBound=-1e-12,
                scalarUnitK=spec["condition"]["scalarUnitK"], residualWindow=window.tolist(), frames=frames,
                first=first, row=row, far=far, originalFaces=ref["faces"], originalAreaM2=ref["physicalAreaM2"],
                acceptedFaces=0, acceptedAreaM2=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="首次生成工作台小数据，已有文件只复核")
    args = parser.parse_args()
    report = audit()
    output = DATA/"supportedFlux65.json"
    if args.write and not output.exists():
        with output.open("x", encoding="utf8", newline="\n") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
    assert read(output) == report
    print(sha(output))
