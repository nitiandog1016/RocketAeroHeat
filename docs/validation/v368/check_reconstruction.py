"""归档两种原生读数，或仅用随 Git 发布的小表重算对照结论。"""
import argparse
from pathlib import Path
import numpy as np
from reconstruction import DOC, MATCHED, OUT, PREVIOUS, read, save, sha, terminal
from read_assembly import IDS, RES, assembly, numbers

METHODS = {"weighted": "inverse-distance weighted Least-Squares", "gg": "Green-Gauss"}
BASELINE = DOC.parent / "v367"
PRIMITIVES = ["T", "u", "v", "p", "rho", "H"]
FACE_COLUMNS = (["edge", "left", "right", "nxArea", "nrArea", "firstOrderWeight", "musclRamp"] +
                [s + k for s in ("rawL", "rawR", "faceL", "faceR") for k in PRIMITIVES] +
                [s + k for s in ("limL", "limR") for k in PRIMITIVES[:4]] +
                ["centerMass", "centerEnergy", "roeMass", "roeEnergy"])


def write_array(path, a, names):
    with path.open("x", encoding="ascii", newline="\n") as f:
        np.savetxt(f, a, fmt="%.17g", delimiter=",", header=",".join(names), comments="")


def archive(folder):
    folder.mkdir(parents=True, exist_ok=False)
    reference = np.genfromtxt(BASELINE / "evidence/matched_nodes.csv", names=True, delimiter=",")
    baseline_edges = np.genfromtxt(BASELINE / "evidence/matched_edges.csv", names=True, delimiter=",")
    records = {}
    rhs = [IDS]
    for method, printed in METHODS.items():
        sample = OUT / (method + "01")
        receipt = terminal(sample)
        assert all(sha(Path(k)) == h for k, h in receipt["inputSHA256"].items())
        log = (sample / "debugger.log").read_text(errors="replace")
        assert "Gradient for upwind reconstruction: " + printed + "." in log
        assert "Gradient for viscous and source terms: inverse-distance weighted Least-Squares." in log
        assert (sample / "bound_state.bin").read_bytes() == (MATCHED / "bound_state.bin").read_bytes()
        rows = np.loadtxt(sample / "face_samples.csv", delimiter=",")
        rows = rows[np.argsort(rows[:, 0])]
        assert rows.shape == (5, len(FACE_COLUMNS)) and np.isfinite(rows).all()
        report, nodes, edges = assembly(sample, 0)
        indices = np.searchsorted(edges["edge"], rows[:, 0])
        for col, field in ((0,"edge"), (1,"left"), (2,"right"), (40,"central"), (42,"roe")):
            np.testing.assert_allclose(rows[:, col], edges[field][indices], rtol=0, atol=1e-15)
        np.testing.assert_array_equal(edges["edge"], baseline_edges["edge"])
        visc = np.max(abs((edges["total"]-edges["roe"]) - (baseline_edges["total"]-baseline_edges["roe"])))
        assert visc < 1e-12
        np.testing.assert_array_equal(nodes[:,0], reference["point"])
        hist = numbers(sample / "history_00001.csv")
        assert hist["Inner_Iter"].tolist() == [0] and hist["Time_Step"].tolist() == [50e-6]
        logs = {k: float(hist[k][0]) for k in RES}
        if method == "weighted":
            np.testing.assert_allclose(nodes[:,-1], reference["implicitRHS"], rtol=0, atol=1e-12)
            previous = numbers(MATCHED / "history_00001.csv")
            assert max(abs(v-previous[k][0]) for k,v in logs.items()) < 1e-8
        rhs.append(nodes[:,-1])
        write_array(folder / (method + "_faces.csv"), rows, FACE_COLUMNS)
        records[method] = dict(globalLogRMS=logs, fluxClosure=report["checks"], maxViscousChange=float(visc),
                               gradient=printed, runSeconds=receipt["elapsedSeconds"], peakRSSBytes=receipt["peakFamilyRSSBytes"],
                               finalFrameSHA256=sha(sample/"restart_00001.csv"), boundStateSHA256=sha(sample/"bound_state.bin"),
                               cfgSHA256=sha(sample/"case.cfg"), historySHA256=sha(sample/"history_00001.csv"))
    write_array(folder / "nodes.csv", np.column_stack(rhs), ["point", *METHODS])
    save(folder / "runs.json", records)
    save(folder / "status.json", compare(folder))


def compare(folder):
    """由原始双侧面值还原中心质量/能量通量，再执行事前三项选择条件。"""
    records = read(folder / "runs.json")
    n = np.genfromtxt(folder / "nodes.csv", names=True, delimiter=",")
    np.testing.assert_array_equal(n["point"], IDS)
    wls_nodes = np.genfromtxt(BASELINE / "evidence/matched_nodes.csv", names=True, delimiter=",")
    np.testing.assert_allclose(n["weighted"], wls_nodes["implicitRHS"], atol=1e-12, rtol=0)
    result = {}
    baseline = None
    for method in METHODS:
        e = np.genfromtxt(folder / (method + "_faces.csv"), names=True, delimiter=",")
        np.testing.assert_array_equal(np.sort(np.minimum(e["left"], e["right"])), np.arange(3825,3830))
        assert np.all(abs(e["right"]-e["left"]) == 1)
        if baseline is None:
            baseline = e
        else:
            for field in ["edge","left","right","nxArea","nrArea", *FACE_COLUMNS[7:19]]:
                np.testing.assert_array_equal(e[field], baseline[field])
        fluxes = {}
        for prefix in ("raw", "face"):
            sides = [prefix + p for p in ("L","R")]
            mass = [e[s+"rho"] * (e[s+"u"] * e["nxArea"] + e[s+"v"] * e["nrArea"]) for s in sides]
            energy = [m*e[s+"H"] for m,s in zip(mass, sides)]
            fluxes[prefix] = (sum(mass)/2, sum(energy)/2)
        errors = [float(max(abs(x-e[k]))) for x,k in zip(fluxes["face"], ("centerMass","centerEnergy"))]
        assert max(errors) < 1e-12
        point = 3828
        sign = (e["left"] == point).astype(int) - (e["right"] == point).astype(int)
        assert sum(abs(sign)) == 2
        j = np.argmax(abs(n[method]))
        result[method] = dict(pointResidualAbs=float(abs(n[method][n["point"] == point][0])),
                              fourColumnRhsNorm=float(np.linalg.norm(n[method])),
                              worstGlobalLogRMS=max(records[method]["globalLogRMS"].values()),
                              localMaxPoint=int(n["point"][j]), localMaxAbs=float(abs(n[method][j])),
                              rawCenterEnergy=float(sign@fluxes["raw"][1]), faceCenterEnergy=float(sign@e["centerEnergy"]),
                              rawCenterMass=float(sign@fluxes["raw"][0]), faceCenterMass=float(sign@e["centerMass"]),
                              maxMassEnergyClosure=errors, maxFirstOrderWeight=float(max(e["firstOrderWeight"])))
    tests = {k: result["gg"][k] < result["weighted"][k] for k in ("pointResidualAbs","fourColumnRhsNorm","worstGlobalLogRMS")}
    s = read(BASELINE / "evidence/status.json")
    return dict(**{k:s[k] for k in ("faces","surfaceSHA256","mach","altitudeM","wallGivenK")},
                targetPoint=point, methods=result, trialChecks=tests, admitStepTrial=all(tests.values()), wholeBodyPowerW=None)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path)
    a = p.parse_args()
    if a.archive:
        archive(a.archive)
    else:
        evidence = DOC / "evidence"
        actual = compare(evidence)
        saved = read(evidence / "status.json")
        # 不同 BLAS 归约顺序的末位误差不影响范数与准入条件。
        for key in METHODS:
            for field in actual["methods"][key]:
                np.testing.assert_allclose(actual["methods"][key][field], saved["methods"][key][field],
                                           rtol=2e-14, atol=1e-25)
        assert {k:v for k,v in actual.items() if k!="methods"} == {k:v for k,v in saved.items() if k!="methods"}
        print(actual["trialChecks"])
