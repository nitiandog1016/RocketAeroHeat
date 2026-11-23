"""读取 Short 银球已验收的两个同刻场，量化回流与速度恢复带的原生网格分辨率。"""
import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from runEmbeddedSU2 import sha


ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "cases/short_sphere"
SECTIONS = (-.6, -.75, -1., -1.25, -1.5, -2.)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def meridian(mesh):
    """用射线分组给原节点排序，用相邻的真实棱柱及六面体核对每一个母面。"""
    with mesh.open(encoding="ascii") as f:
        assert next(f).strip() == "NDIME= 3"
        count = int(next(f).split("=")[1])
        cells = [list(map(int, next(f).split())) for _ in range(count)]
        points = int(next(f).split("=")[1])
        xyz = np.array([list(map(float, next(f).split())) for _ in range(points)])
    in_plane = (np.abs(xyz[:, 1]) < 1e-12) & (xyz[:, 2] >= -1e-13)
    ids = np.flatnonzero(in_plane)
    angle = np.round(np.arctan2(xyz[ids, 2], -xyz[ids, 0]), 10)
    distance = np.linalg.norm(xyz[ids], axis=1)
    rays = [ids[angle == a][np.argsort(distance[angle == a])] for a in np.unique(angle)]
    nodes = np.array(rays)
    assert nodes.shape == (49, 33)
    shared = Counter()
    for kind, *cell in cells:
        assert (kind, len(cell)) in ((12, 8), (13, 6))
        face = tuple(sorted(i for i in cell if in_plane[i]))
        if len(face) == 4:
            shared[face] += 1
    expected = {tuple(sorted(nodes[i:i+2, j:j+2].ravel()))
                for i in range(len(nodes)-1) for j in range(nodes.shape[1]-1)}
    assert set(shared) == expected and set(shared.values()) == {2}
    np.testing.assert_allclose(np.linalg.norm(xyz[nodes[:, 0]], axis=1), .5, atol=1e-13, rtol=0)
    return xyz, nodes


def axis_result(xyz, ray, velocity):
    distance = -xyz[ray, 0]-.5
    speed = -velocity[ray, 0]
    assert abs(speed[0]) < 1e-12 and np.all(np.diff(distance) > 0)
    result = dict(nodeIDs=ray.tolist(), distanceD=distance.tolist(), downstreamSpeedOverU=speed.tolist(),
                  lengthD=None, bracketD=[], poorAxialResolution=False)
    if speed[1] >= 0:
        result["status"] = "no-wall-connected-reversal"
        return result
    forward = np.flatnonzero(speed[1:] >= 0)
    if not len(forward):
        result["status"] = "extends-to-boundary"
        return result
    lo = int(forward[0])
    hi = lo+1
    span = distance[hi]-distance[lo]
    length = distance[lo] - span*speed[lo]/(speed[hi]-speed[lo])
    result.update(status="bracketed", lengthD=length, bracketD=distance[lo:hi+1].tolist(),
                  reverseNodes=lo, poorAxialResolution=bool(span > .1 and span > .1*length))
    return result


def recovery_band(radii, speed, axial_spans):
    """沿轴速亏损由谷底向外找首次恢复位置，记录跨带的独立供值线。"""
    core = np.flatnonzero(radii <= .6)
    i = core[np.argmin(speed[core])]
    low = speed[i]
    if low >= 1:
        return dict(status="no-deficit")
    bounds = []
    for fraction in (.1, .9):
        target = low + fraction*(1-low)
        up = np.flatnonzero((speed[i:-1] < target) & (speed[i+1:] >= target))
        if not len(up):
            return dict(status="unbracketed-recovery")
        j = i+up[0]
        bounds.append(radii[j]+(radii[j+1]-radii[j])*(target-speed[j])/(speed[j+1]-speed[j]))
    a, b = bounds
    segments = np.flatnonzero((radii[:-1] < b) & (radii[1:] > a))
    gaps = np.diff(radii)[segments]
    spans = axial_spans[segments[0]:segments[-1]+2]
    count = int(np.count_nonzero((a < radii) & (radii < b)))
    relative = (b-a)/max(gaps)
    return dict(status="bracketed", minSpeedOverU=low, minRadiusD=radii[i], r10D=a, r90D=b, widthD=b-a,
                innerRays=count, minRadialGapD=min(gaps), maxRadialGapD=max(gaps),
                minAxialSpanD=min(spans), maxAxialSpanD=max(spans), widthOverMaxGap=relative,
                sparse=bool(count < 4 or relative < 4))


def section(x, xyz, rays, velocity):
    samples = []
    for ray in rays:
        depth = -xyz[ray, 0]
        if not (depth[0] <= -x < depth[-1]):
            continue
        j = np.searchsorted(depth, -x, side="right")-1
        a, b = ray[j:j+2]
        span = depth[j+1]-depth[j]
        weight = (-x-depth[j])/span
        radius = (1-weight)*xyz[a, 2]+weight*xyz[b, 2]
        speed = -(1-weight)*velocity[a, 0]-weight*velocity[b, 0]
        samples.append((radius, speed, span, int(a), int(b)))
    radii, speed, spans, first, second = map(np.asarray, zip(*samples))
    assert radii[0] == 0 and np.all(np.diff(radii) > 0) and len(speed) > 2
    return dict(xD=x, radiusD=radii.tolist(), downstreamSpeedOverU=speed.tolist(), axialSpanD=spans.tolist(),
                sourceNode0=first.tolist(), sourceNode1=second.tolist(), band=recovery_band(radii, speed, spans))


def wall_stations(xyz, rays, velocity, table):
    # 极角从迎风点增大；在真实首层的各端点投影切向速度，再按实验角度插值。
    ids = rays[::-1, 1]
    phi = np.arctan2(xyz[ids, 2], xyz[ids, 0])
    tangent = -np.sin(phi)*velocity[ids, 0]+np.cos(phi)*velocity[ids, 2]
    stations = []
    for row in table[-2:]:
        angle = float(row["angleDeg"])
        j = np.searchsorted(np.degrees(phi), angle)-1
        assert 0 <= j < len(ids)-1
        value = np.interp(angle, np.degrees(phi), tangent)
        stations.append(dict(angleDeg=angle, firstTangentialSpeedOverU=value, adjacentFirstNodeIDs=ids[j:j+2].tolist(),
                             heightD=float(np.linalg.norm(xyz[ids[j]])-.5), heatInWm2=float(row["computedFluxInWm2"]),
                             computedNu=float(row["computedNu"]), measuredNu=float(row["reportedNu"])))
    return stations


def endpoint(work, branch, xyz, rays, mesh_hash):
    key = read(CASES / ("wake_"+branch+".json"))
    archive = CASES / ("wake_"+branch+".zip")
    assert sha(archive) == key["archiveSHA256"] and sha(CASES/"test65.json") == key["caseSHA256"]
    with ZipFile(archive) as z:
        summary = json.loads(z.read("summary.json"))
        stations = list(csv.DictReader(io.StringIO(z.read("stations.csv").decode("ascii"))))
    step = summary["physicalStep"]
    d = work/branch/f"step{step:03d}"
    frame = d/f"restart_{step:05d}.csv"
    proc = read(d/"debugger_process.json")
    native = read(d/"native.json")
    readout = read(d/"readout/accepted.json")
    assert proc["status"] == "terminal" and proc["debuggerExitCode"] == 0 and all(proc["inputsStable"].values())
    assert not (proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
    assert mesh_hash == summary["meshSHA256"] == sha(d/"mesh.su2") and summary["timeDU"] == 8
    assert sha(frame) == native["frame"] == readout["frame"] == summary["frameSHA256"]
    assert native["old"] == summary["previousFrameSHA256"] == sha(d/f"seed_{step-1:05d}.csv")
    assert native["history"] == sha(d/f"history_{step:05d}.csv") and max(native["window"].values()) < -8
    assert sha(d/"case.cfg") == proc["inputSHA256"][str(d/"case.cfg")]
    log = (d/"debugger.log").read_text("ascii", errors="replace")
    assert "First order integration in space." in log and "Exit Success (SU2_CFD)" in log
    config = dict(line.split("=", 1) for line in (d/"case.cfg").read_text("ascii").splitlines() if "=" in line)
    keys = ["CONV_NUM_METHOD_FLOW", "MUSCL_FLOW", "SLOPE_LIMITER_FLOW", "CONV_NUM_METHOD_SPECIES",
            "MUSCL_SPECIES", "NUM_METHOD_GRAD", "SPECIES_CLIPPING"]
    scheme = {k:config[k].strip() for k in keys}
    data = np.loadtxt(frame, skiprows=1, delimiter=",", usecols=range(8))
    np.testing.assert_array_equal(data[:, 0], np.arange(len(xyz)))
    np.testing.assert_allclose(data[:, 1:4], xyz, rtol=0, atol=1e-12)
    velocity = data[:, 5:8]
    assert np.isfinite(velocity).all()
    return dict(frameSHA256=native["frame"], caseConfigSHA256=sha(d/"case.cfg"), scheme=scheme,
                timeS=summary["elapsedPhysicalTimeS"], axis=axis_result(xyz, rays[0], velocity),
                sections=[section(x,xyz,rays,velocity) for x in SECTIONS],
                stations=wall_stations(xyz,rays,velocity,stations))


def analyze(work, output):
    mesh = work/"fine/step064/mesh.su2"
    mesh_hash = sha(mesh)
    xyz, rays = meridian(mesh)
    cases = {b:endpoint(work,b,xyz,rays,mesh_hash) for b in ("coarse","fine")}
    assert cases["coarse"]["timeS"] == cases["fine"]["timeS"]
    data = dict(meshSHA256=mesh_hash, caseSHA256=sha(CASES/"test65.json"), test=65,
                meridianQuads=(rays.shape[0]-1)*(rays.shape[1]-1), **cases)
    encoded = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",",":")).encode("utf8")
    assert len(encoded) < 100*1024**2
    with output.open("xb") as f:
        f.write(encoded)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    analyze(args.work.resolve(), args.output.resolve())
