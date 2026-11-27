"""固定已发布的低速动量帧，独立积分储热，逐原面比较两套热网并生成工作台复核包。"""
import argparse
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from refined import DOC, read, save, sha
import numpy as np
from finalize import audit_steps
from simpleLauncherTunnel import CASE, geometry


def stored_power(previous, current, mesh, dt, scale):
    """沿方位求均值；每个子午三角形用二次精确求积积分半径与温差之积。"""
    with np.load(previous) as src, np.load(current) as dst:
        before = src["theta"]
        if before.ndim == 3:
            before = before[0]
        change = (dst["theta"] - before).mean(axis=0)[mesh["triangles"]]
    points = mesh["xy"][mesh["triangles"]]
    area = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0]) / 2
    assert np.all(area > 0)
    integral = 0.
    for i in range(3):
        weights = np.full(3, 1/6)
        weights[i] = 2/3
        integral += area @ ((change @ weights) * (points[:, :, 1] @ weights)) / 3
    return float(integral * scale / dt)


def space_difference(coarse, fine, coarse_heat, fine_heat, regions):
    """所有相对差均以细网为分母，热流取原 STL 面积权重。"""
    area = fine[:, 1]
    np.testing.assert_array_equal(coarse[:, :3], fine[:, :3])
    q0, q1 = -coarse[:, 3], -fine[:, 3]
    region_errors = {}
    for name, s in regions.items():
        region_errors[name] = float(abs((area[s] @ q0[s]) / (area[s] @ q1[s]) - 1))
    nose = regions["sphericalNose"]
    coarse_peak, fine_peak = [float(q[nose].max()) for q in (q0, q1)]
    rings = np.array([h["baseAnnulusW"] for h in (coarse_heat, fine_heat)])
    np.testing.assert_allclose(rings.sum(axis=1),
                               [h["regionW"]["base"] for h in (coarse_heat, fine_heat)], atol=1e-14, rtol=0)
    ring_errors = abs(rings[0] / rings[1] - 1)
    total = float(abs((area @ q0) / (area @ q1) - 1))
    rms = float(np.sqrt((area @ (q0-q1)**2) / (area @ q1**2)))
    peak_error = abs(coarse_peak/fine_peak-1)
    return dict(outwardDifference=total, rms=rms, regionDifferences=region_errors,
                coarseAndFineNosePeakWm2=[coarse_peak, fine_peak], nosePeakDifference=peak_error,
                nosePeakFaceIDs=[int(fine[nose, 0][np.argmax(q[nose])]) for q in (q0, q1)],
                coarseAndFineBaseAnnulusW=rings.tolist(), baseAnnulusDifferences=ring_errors.tolist(),
                gates=dict(total=total < .05, faceRMS=rms < .05, sixRegions=max(region_errors.values()) < .1,
                           nosePeak=peak_error < .1, fourAnnuli=bool(max(ring_errors) < .1)))


def compare(work, flow, coarse, dest):
    definition, tri, diameter, _ = geometry()
    regions = {g["name"]: slice(g["firstFace"]-1, g["lastFace"]) for g in definition["regions"]}
    condition = read(CASE/"wake/condition.json")
    scale = (2*np.pi * diameter**3 * condition["densityKgM3"] * condition["specificHeatJkgK"] *
             (condition["wallTemperatureK"]-condition["ambientTemperatureK"]))
    reference = read(DOC.parent/"v354/evidence.json")["cases"]
    flow_index = read(CASE/"wake/flow_index.json")
    fine = work/"heat"
    fine_tables, records, time_errors = audit_steps(flow, fine)
    initial = read(work/"initial/results.json")
    mesh_file = work/"steady/temperature.npz"
    assert initial["accepted"] and all(initial["gates"].values())
    assert sha(mesh_file) == initial["sourceHeatMeshSHA256"]
    assert sha(work/"initial/fields.npz") == initial["fieldsSHA256"]
    with np.load(mesh_file) as source:
        mesh = {k: source[k] for k in ("xy", "triangles")}
    coarse_tables, spatial, steps = {}, {}, {}
    for name, record in records.items():
        f, h = record["flow"], record["heat"]
        assert all(f[k] == v for k, v in flow_index[name].items())
        old = read(coarse/name/"heat.json")
        frozen = reference[name]
        assert frozen["flowSHA256"] == old["sourceFrameSHA256"] == h["sourceFrameSHA256"]
        assert frozen["fieldSHA256"] == sha(coarse/name/"thermal.npz") == old["fieldSHA256"]
        assert frozen["wallSHA256"] == sha(coarse/name/"wall.csv") == old["wallSHA256"]
        assert old["elapsedPhysicalS"] == h["elapsedPhysicalS"] and h["heatMeshSHA256"] == initial["sourceHeatMeshSHA256"]
        assert record["heatPeakRSSBytes"] < 4*1024**3 and record["heatSeconds"] < 240
        previous = work/"initial/fields.npz" if f["step"] == 1 else fine/(f["branch"]+"_half_1")/"thermal.npz"
        assert sha(previous) == h["initialThetaSHA256"]
        power = stored_power(previous, fine/name/"thermal.npz", mesh, f["dtDU"]*diameter/.05, scale)
        assert abs(power-h["storageW"]) < 1e-14 and abs(h["outwardW"]-power-h["farW"]) < 1e-12
        coarse_tables[name] = np.loadtxt(coarse/name/"wall.csv", delimiter=",", skiprows=1)
        spatial[name] = space_difference(coarse_tables[name], fine_tables[name], old, h, regions)
        steps[name] = dict(heat=h, independentStorageW=power, elapsedSeconds=record["heatSeconds"],
                           peakBytes=record["heatPeakRSSBytes"], coarseWallSHA256=old["wallSHA256"],
                           coarseThermalSHA256=old["fieldSHA256"])

    base = regions["base"]
    mid = tri[base].mean(axis=1)
    angle = np.arctan2(mid[:, 2], mid[:, 1])
    assert len(angle) == 48
    phase = np.exp(-1j*angle)
    harmonics, coeff = {}, {}
    for grid, tables in (("coarse", coarse_tables), ("fine", fine_tables)):
        for suffix in ("full_1", "half_2"):
            response = tables["control_"+suffix][base, 3] - tables["disturbed_"+suffix][base, 3]
            harmonic = 2*np.mean(response*phase)
            key = grid+"_"+suffix
            coeff[key] = harmonic
            harmonics[key] = dict(realWm2=harmonic.real, imagWm2=harmonic.imag,
                                  magnitudeWm2=abs(harmonic), phaseDeg=np.angle(harmonic, deg=True))
    target = coeff["fine_half_2"]
    assert abs(target) > 0
    space_error = abs(coeff["coarse_half_2"]-target)/abs(target)
    time_error = abs(coeff["fine_full_1"]-target)/abs(target)
    harmonics.update(relativeSpace=space_error, relativeTime=time_error, identifiable=bool(max(space_error,time_error)<.2))
    name = "disturbed_half_2"
    h = records[name]["heat"]
    with ZipFile(CASE/"wake/wind_tunnel.zip") as z:
        primary = json.loads(z.read("summary.json"))
    assert primary["condition"] == condition and primary["wallSHA256"] == steps[name]["coarseWallSHA256"]
    assert primary["thermalFrameSHA256"] == steps[name]["coarseThermalSHA256"]
    assert primary["sourceFrameSHA256"] == h["sourceFrameSHA256"]
    assert primary["surfaceSHA256"] == definition["fileSHA256"] and primary["areaM2"] == definition["areaM2"]
    summary = dict(condition=condition, surfaceSHA256=primary["surfaceSHA256"], faces=len(tri), areaM2=primary["areaM2"],
                   coarseWallSHA256=primary["wallSHA256"], coarseThermalSHA256=primary["thermalFrameSHA256"],
                   flowSHA256=h["sourceFrameSHA256"], fineThermalSHA256=h["fieldSHA256"], fineWallSHA256=h["wallSHA256"],
                   heatMeshSHA256=h["heatMeshSHA256"], fineInitialSHA256=initial["fieldsSHA256"],
                   outwardW=h["outwardW"], spatial=spatial[name], baseHarmonic=harmonics)
    dest.mkdir(parents=True,exist_ok=False)
    bundle = dest/"wake_heat_grid.zip"
    with ZipFile(bundle, "x", ZIP_DEFLATED) as z:
        z.writestr("comparison.json", json.dumps(summary,ensure_ascii=False,allow_nan=False,indent=2))
        z.write(fine/name/"wall.csv","fine_wall.csv")
    save(dest/"evidence.json", dict(initial=initial, steps=steps, spatial=spatial, time=time_errors,
                                    baseHarmonic=harmonics, comparisonSHA256=sha(bundle)))
    print(json.dumps(dict(spatial=spatial[name], harmonic=harmonics, time=time_errors),indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work","flow","coarse","dest"):
        parser.add_argument(key,type=Path)
    a = parser.parse_args()
    compare(a.work,a.flow,a.coarse,a.dest)
