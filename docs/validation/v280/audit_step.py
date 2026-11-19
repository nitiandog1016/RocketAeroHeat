"""复核同一步的两条历史、原面功率及终态向量；按准确列名保存可读回初值。"""
import argparse
import gzip
from pathlib import Path
import re
import shutil
import time

import numpy as np
from scipy.spatial import cKDTree

from compare_step import HERE, OLD, ROOT, RUN, save, sha
from audit_tunnel import project, read, table, vtk_tokens

BASE = HERE.parent / "v279/evidence"
RMS = [f"rms[{k}]" for k in ("P", "U", "V", "W", "h", "k", "w", "LM_1", "LM_2")]


def compare_history():
    old_keys, fixed = table(OLD / "history_00001.csv")
    keys, adaptive = table(RUN / "history_00001.csv")
    def series(key):
        return adaptive[:, keys.index(key)]
    n = len(adaptive)
    assert np.isfinite(adaptive).all() and np.all(series("Time_Iter") == 1)
    assert np.array_equal(series("Inner_Iter"), np.arange(n)) and len(fixed) == 120
    assert np.all(series("Outer_Iter") == 0) and len(keys) == len(set(keys))
    np.testing.assert_array_equal(adaptive[0, [keys.index(k) for k in RMS]], fixed[0, [old_keys.index(k) for k in RMS]])
    visible = np.loadtxt(RUN / "visible.csv", delimiter=",", skiprows=1)
    elapsed = read(OLD / "terminal.json")["elapsedSeconds"]
    bracket = np.searchsorted(visible[:, 0], elapsed)
    assert 0 < bracket < len(visible) and visible[-1, 1] == n
    samples = visible[bracket-1:bracket+1]
    def record(i):
        return {k: float(series(k)[i]) for k in RMS + ["HF", "Min CFL", "Avg CFL", "Max CFL"]}
    lin = {}
    for name, suffix in (("flow", ""), ("SST", "Turb"), ("LM", "Trans")):
        residual, it = series("LinSolRes" + suffix), series("LinSolIter" + suffix)
        lin[name] = dict(finalLog10=residual[-1], finalIterations=int(it[-1]),
                         aboveTargetCount=int((residual > -3).sum()), atLimitCount=int((it >= 16).sum()),
                         maxLog10=float(residual.max()))
    spans = {k: float(np.ptp((series(k) - series("Pref"))[-30:])) for k in ("Pfront", "Pside", "Prear")}
    spans["HF_W"] = float(np.ptp(series("HF")[-30:]))
    result = dict(innerIterations=n, baselineElapsedS=elapsed,
                  fixed120={k: float(fixed[-1, old_keys.index(k)]) for k in RMS + ["HF"]},
                  adaptive120=record(119), sameElapsedVisible=[
                      dict(elapsedS=s, completeRows=int(i), history=record(int(i)-1)) for s, i in samples],
                  adaptiveFinal=record(-1), linear=lin, last30Range=spans,
                  passesNineResiduals=bool(max(series(k)[-1] for k in RMS) < -8),
                  terminalResiduals={k: float(series(k)[-1]) for k in RMS},
                  simultaneousPassCount=int(np.all(adaptive[:, [keys.index(k) for k in RMS]] < -8, axis=1).sum()))
    return result, dict(zip(keys, adaptive[-1])), fixed, adaptive, old_keys, keys


def wall_fields(volume, keys, condition, probes, monitor):
    wall = np.load(RUN / "wall.npz")
    ids, span = wall["nodeIDs"], wall["x"]
    xyz = volume[:, 1:4]
    center = xyz[ids]
    np.testing.assert_allclose(center, wall["xyz"], rtol=0, atol=1e-15)
    surface, fields = vtk_tokens(RUN / "cylinder_surface_00001.vtk")
    distance, idx = cKDTree(surface).query(center)
    assert len(np.unique(idx)) == len(surface) == 1800 and distance.max() < 1e-7
    def col(name):
        return volume[:, keys.index(name)]
    t = col("Temperature")
    np.testing.assert_allclose(t, col("Enthalpy") / condition["cpJkgK"] + 298.15, rtol=0, atol=2e-12)
    fluid_wall_error = np.max(abs(t[ids] - condition["wallTemperatureK"]))
    for name in ("Velocity_x", "Velocity_y", "Velocity_z"):
        assert abs(col(name)[ids]).max() < 1e-10
    normals = center.copy()
    normals[:, :, 0] = 0
    normals /= np.linalg.norm(normals, axis=2)[:, :, None]
    dn = np.sum((xyz[ids + 72] - center) * normals, axis=2)
    np.testing.assert_allclose(dn, 8e-6, rtol=0, atol=1e-14)
    flux = condition["conductivityWmK"] * (t[ids + 72] - condition["wallTemperatureK"]) / dn
    native = col("Heat_Flux")[ids]
    # 原生 VTK 标量只写小数点后六位科学记数法，逐项以相同精度比对。
    encoded = np.char.mod("%.6e", native).astype(float)
    np.testing.assert_array_equal(encoded, fields["Heat_Flux"][:, 0][idx])
    error = float(np.max(abs(native[1:-1] - flux[1:-1])))
    assert error < 1e-6 and np.all(native[[0, -1]] == 0)
    widths = np.diff(span)
    node_weights = (np.r_[0, widths] + np.r_[widths, 0]) / 2
    ring = node_weights @ flux.sum(axis=1)
    power = ring * .033 * np.sin(2*np.pi/72)
    assert abs(power - monitor["HF"]) < 1e-7
    for key, i in probes.items():
        np.testing.assert_allclose(col("Pressure")[i], monitor[key], rtol=1e-8, atol=1e-7)
    ref = col("Pressure")[probes["Pref"]]
    pressure = condition["referencePressurePa"] + project(span, col("Pressure")[ids] - ref)
    heat = project(span, flux)
    header, original = table(BASE / "trial_faces.csv")
    assert len(original) == len(heat) == 2160 and np.all(pressure > 0)
    mapped = original[:, 4] @ heat
    assert abs(mapped - ring * 360 * 2 * .033 * np.sin(np.pi/360) / 72) < 1e-10
    summary = dict(sharedNodes=144, nonSharedNativeHeatMaxDifferenceWm2=error, nativeWallPowerW=monitor["HF"],
                   wallNodeFluidTRangeK=[t[ids].min(), t[ids].max()], maxWallNodeFluidTMinusGivenAbsK=fluid_wall_error,
                   wallNodeAgreementAt1eMinus9K=bool(fluid_wall_error < 1e-9),
                   rebuiltDualWallPowerW=power, projectedTrialPowerW=mapped, trialFaceCount=len(heat),
                   trialPressureRangePa=[min(pressure), max(pressure)], trialHeatRangeWm2=[min(heat), max(heat)],
                   sharedCylinderHeatRangeWm2=[flux[[0, -1]].min(), flux[[0, -1]].max()])
    return header, np.column_stack((original[:, :5], pressure, heat)), summary


def locations(volume, keys, monitor):
    """逐区域列出终态线性系统右端量峰值。独立与原生最大值及 RMS 数值核对。"""
    columns = ["Pressure", "Velocity_x", "Velocity_y", "Velocity_z", "Enthalpy", "TKE",
               "Omega", "LM_intermittency", "LM_RE_THETA_T"]
    i = np.arange(len(volume))
    x, y, z = volume[:, 1:4].T
    regions = {"whole": np.ones(len(i), dtype=bool), "mountFirstTwoLayers": (i//2952 <= 2) | (i//2952 >= 22),
               "cylinderFirstTwoLayers": (i//72 % 41 <= 2), "downstreamWake": (y < -.033) & (abs(z) < .066)}
    rows = {}
    for key, column in zip(RMS, columns):
        v = abs(volume[:, keys.index("Residual_" + column)])
        rows[key] = {}
        for region, mask in regions.items():
            ids = i[mask]
            j = int(ids[v[mask].argmax()])
            rows[key][region] = dict(node=j, xyzM=volume[j, 1:4].tolist(), absRhs=v[j])
        rows[key]["historyLogMax"] = monitor[key.replace("rms", "max")]
        rows[key]["volumeLogMax"] = float(np.log10(np.max(v)))
        rows[key]["historyLogRms"] = monitor[key]
        rows[key]["volumeLogRms"] = float(np.log10(np.sqrt(np.mean(v*v))))
    return rows


def plot(path, fixed, adaptive, old_keys, keys):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), constrained_layout=True)
    for ax, key, label in zip(axes.flat, ["rms[h]", "rms[w]", "rms[LM_2]", "HF"],
                              ["Enthalpy", "Specific dissipation rate", "Transition Reynolds number", "Cylinder heat rate"]):
        ax.plot(np.arange(1, len(adaptive)+1), adaptive[:, keys.index(key)], label="Adaptive CFL")
        ax.plot(np.arange(1, len(fixed)+1), fixed[:, old_keys.index(key)], "--", label="Fixed CFL = 3")
        ax.set(title=label, xlabel="Inner iterations", ylabel="W, positive into wall" if key=="HF" else "log10 RMS")
        ax.grid(alpha=.2)
    axes[0, 0].legend()
    fig.savefig(path, dpi=135)
    plt.close(fig)


def audit(out):
    start = time.monotonic()
    old = read(BASE / "trial.json")
    assert sha(BASE / "trial_faces.csv") == old["faceSHA256"]
    assert sha(OLD / "history_00001.csv") == old["sourceHistorySHA256"]
    for folder in (OLD, RUN):
        receipt = read(folder / "terminal.json")
        assert receipt["exitCode"] == 0 and receipt["inputsStable"] and not receipt["stopReason"]
        assert all(sha(folder / p) == h for p, h in receipt["outputs"].items())
    assert all(sha(Path(p)) == h for p, h in read(RUN / "inputs.json").items())
    summary, monitor, fixed, adaptive, old_keys, history_keys = compare_history()
    keys, volume = table(RUN / "restart_00001.csv")
    seed_keys, seed = table(RUN / "seed_00000.csv")
    assert np.isfinite(volume).all() and len(volume) == 73800 and np.array_equal(volume[:, 0], np.arange(73800))
    np.testing.assert_allclose(volume[:, :4], seed[:, :4], atol=1e-15, rtol=0)
    ranges = {k: [volume[:, keys.index(k)].min(), volume[:, keys.index(k)].max()]
              for k in ("Temperature", "Turb_Kin_Energy", "Omega", "LM_gamma", "LM_Re_t")}
    assert min(lo for lo, hi in ranges.values()) > 0
    log = (RUN / "native.log").read_text("utf8", errors="replace")
    assert "Exit Success (SU2_CFD)" in log
    probes = {k: int(i) for k, i in re.findall(r"Probe (\w+) is using global point (\d+)", log)}
    assert probes == dict(Pfront=35424, Pside=35442, Prear=35460, Pref=38304)
    c = read(RUN / "condition.json")
    header, faces, wall = wall_fields(volume, keys, c, probes, monitor)
    spots = locations(volume, keys, monitor)
    out.mkdir(parents=True, exist_ok=False)
    ordered = volume[:, [keys.index(k) for k in seed_keys]]
    restart = out / "restart_with_lm_00001.csv"
    np.savetxt(restart, ordered, delimiter=",", fmt=["%d"]+["%.17g"]*(len(seed_keys)-1),
               header=",".join('"'+k+'"' for k in seed_keys), comments="")
    saved_keys, saved = table(restart)
    assert saved_keys == seed_keys and np.array_equal(saved, ordered)
    np.savetxt(out / "trial_faces.csv", faces, delimiter=",", comments="", fmt=["%d"]+["%.16g"]*6,
               header=",".join(header))
    save(out / "locations.json", spots)
    summary.update(wall, condition=c, physicalDtS=old["physicalDtS"], physicalStepsExecuted=1,
                   sourceFrameSHA256=sha(RUN / "restart_00001.csv"), sourceHistorySHA256=sha(RUN / "history_00001.csv"),
                   faceSHA256=sha(out / "trial_faces.csv"), seedSHA256=sha(RUN / "seed_00000.csv"),
                   extendedRestartSHA256=sha(restart), savedVolumeRanges=ranges, pressureProbeIDs=probes,
                   productionQualified=False, qualifiedAreaM2=0, process=receipt)
    for name in ("history_00001.csv", "visible.csv"):
        with gzip.GzipFile(filename=str(out / (name+".gz")), mode="xb", mtime=0) as archive:
            archive.write((RUN/name).read_bytes())
    for name in ("case.cfg", "inputs.json", "started.json"):
        shutil.copyfile(RUN / name, out / name)
    plot(out / "comparison.png", fixed, adaptive, old_keys, history_keys)
    summary["auditSeconds"] = time.monotonic() - start
    save(out / "trial.json", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("new_output", type=Path)
    audit(parser.parse_args().new_output)
