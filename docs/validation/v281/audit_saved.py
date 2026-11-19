"""检查三份原面存档在已登记缺项面上的证据，独立复核区域页导出。"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs/validation"
OUT = ROOT / "outputs"
COHORT = DOC / "v263/evidence/sidewall_temperature/faces.csv"
ENGINEERING = OUT / "v176_validation/main02/engineering_data.json"
RESEARCH = ENGINEERING.with_name("research_data.json")
MAT = OUT / "v104_validation/actual03/current_case.mat"
KINNEY = MAT.with_name("faces.csv")
ENTHALPY = DOC / "v190/evidence/real02/face_comparison.csv"
EXPECTED = {
    COHORT: "621b9d8d8bb29bb3e3b62c4be062c77b1ca1fed0f03939fc83168087b30ca831",
    ENGINEERING: "9518d1d1f44abb83f60664a8c1d072d25797ce16d2182823900cfdf119110b1b",
    RESEARCH: "7b5d3821ade4f768b77bd1ae327659c82d49bacc5b1e307ebbe12263dce88889",
    MAT: "a9e45da1561e763e5d903fc4c0257ca5716721f0e4c350733419e7a7f7f8fc70",
    KINNEY: "786d3e7634807dfb60d091e00657d67186109915a57c366f58ec0a6ebbeff205",
    ENTHALPY: "1d4064eab49a5424a0ad4466d55df101d8d42ea12fdc71a033a38cd88eb87214",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tables(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f)


def field(snapshot, key, ledger=False):
    source = snapshot["evidenceLedger"] if ledger else snapshot
    return next(f for f in source["fields"] if f["key"] == key)


def check_export(path, snapshot, selected):
    ledger = snapshot["evidenceLedger"]
    identity = ledger["identity"]
    measures = [(col, field(snapshot, k), field(snapshot, k, True))
                for col, k in (("pressure", "pressure"), ("heat_flux", "heatFlux"))]
    counts = dict(rows=0, jointEngineeringFaces=0, jointVerifiedFaces=0,
                  subsetFiniteHeat=0, subsetEngineeringHeat=0, subsetResearchHeat=0,
                  subsetVerifiedHeat=0, subsetVerifiedJoint=0)
    for i, row in enumerate(tables(path)):
        assert int(row["face_id"]) == i+1 and row["region"] == snapshot["region"][i]
        assert abs(float(row["area_m2"])-snapshot["faceArea"][i]) < 1e-14
        for key, name in (("geometry", "geometrySHA256"), ("condition", "conditionSHA256"), ("source", "sourceSHA256")):
            assert row[key+"_sha256"] == identity[name]
        assert row["condition"] == snapshot["resultCondition"]["caption"]
        assert int(row["synthetic"]) == identity["synthetic"] and int(row["source_traceable"]) == identity["sourceTraceable"]
        for key, raw, evidence in measures:
            value = raw["values"][i]
            if value is None:
                assert row[key] == ""
            else:
                assert float(row[key]) == value
            assert row[key+"_unit"] == raw["unit"]
            assert row[key+"_method"] == raw["methodNames"][raw["methodIndex"][i]]
            assert row[key+"_availability"] == ledger["availabilityNames"][evidence["availability"][i]]
            assert row[key+"_qualification"] == ledger["qualificationNames"][evidence["qualification"][i]]
        pair = all(row[k+"_availability"] == "engineering" for k in ("pressure", "heat_flux"))
        passed = all(row[k+"_qualification"] == "passed" and row[k+"_availability"] in ("engineering","research")
                     for k in ("pressure", "heat_flux"))
        assert int(row["joint_engineering"]) == pair and int(row["joint_verified"]) == passed
        counts["rows"] += 1
        counts["jointEngineeringFaces"] += pair
        counts["jointVerifiedFaces"] += passed
        if i+1 in selected:
            counts["subsetFiniteHeat"] += row["heat_flux"] != ""
            counts["subsetEngineeringHeat"] += row["heat_flux_availability"] == "engineering"
            counts["subsetResearchHeat"] += row["heat_flux_availability"] == "research"
            counts["subsetVerifiedHeat"] += row["heat_flux_qualification"] == "passed"
            counts["subsetVerifiedJoint"] += passed
    assert counts["rows"] == len(snapshot["faces"])
    return counts


def run(work, dest):
    start = time.monotonic()
    assert all(sha(p) == h for p, h in EXPECTED.items())
    eng = json.loads(ENGINEERING.read_bytes())
    research = json.loads(RESEARCH.read_bytes())
    pts = {int(r["faceID"]): r for r in tables(COHORT)}
    assert len(pts) == 936 and len(eng["faces"]) == 19752
    assert eng["vertices"] == research["vertices"] and eng["faces"] == research["faces"]
    np.testing.assert_allclose(eng["faceArea"], research["faceArea"], rtol=0, atol=1e-14)
    cond = eng["resultCondition"]
    other = research["resultCondition"]
    for key in ("Mach","alphaDeg","betaDeg","freestreamPressurePa","freestreamTemperatureK","wallUniformValueK"):
        assert cond[key] == other[key]
    states = field(eng, "pressure", True)["availability"]
    for i, r in pts.items():
        assert r["heatIntoWall_kW_m2"] == "" and states[i-1] == 1 and eng["region"][i-1] == "sidewall"
        assert abs(float(r["area_m2"])-eng["faceArea"][i-1]) < 1e-14
    assessed = {k: check_export(work / (k+".csv"), snap, pts) for k, snap in (("engineering", eng), ("research", research))}
    frame = research["evidenceLedger"]["sourceFrame"]
    assert frame["nativeCheckpointIdentity"]["confirmed"] and frame["checkpointIteration"] == 199
    assert research["evidenceLedger"]["identity"]["origin"] == "research"

    # 新旧工程方法比较从各自真实记录取数。Kinney 另查原 MAT 的实际接受掩码。
    data = loadmat(str(MAT), variable_names=["current","fs"], struct_as_record=False, squeeze_me=True)
    mesh, fs = data["current"], data["fs"]
    np.testing.assert_array_equal(mesh.F-1, eng["faces"])
    np.testing.assert_allclose(mesh.V, eng["vertices"], rtol=0, atol=1e-15)
    assert (fs.altitude, fs.Mach, fs.alpha, fs.beta, fs.pressure, fs.temperature) == (
        cond["altitudeM"], cond["Mach"], cond["alphaDeg"], cond["betaDeg"],
        cond["freestreamPressurePa"], cond["freestreamTemperatureK"])
    assert np.all(mesh.heating.wallTemperatureByFace == cond["wallUniformValueK"])
    ref_inputs = json.loads((ENTHALPY.parent / "source_before.json").read_bytes())
    assert ref_inputs[r"outputs\v104_validation\baseline02\current_case.mat"] == eng["evidenceIdentity"]["sourceFrame"]["sha256"]
    subset = np.array(sorted(pts))-1
    candidates = {}
    for label, path, id_key, area_key, col in (("kinney", KINNEY, "face", "area_m2", "heat_flux_W_m2"),
            ("referenceEnthalpy", ENTHALPY, "faceID", "areaM2", "enthalpyReferenceHeatFluxWm2")):
        finite = rows = 0
        for j, r in enumerate(tables(path)):
            assert j+1 == int(r[id_key]) and abs(float(r[area_key])-eng["faceArea"][j]) < 1e-14
            if j+1 in pts:
                finite += math.isfinite(float(r[col]))
                if label == "kinney":
                    assert bool(int(r["accepted"])) == bool(mesh.heating.kinneyAcreage.acceptedFace[j])
                    assert np.array_equal(np.isnan(float(r[col])), np.isnan(mesh.heating.heatFlux[j]))
            rows += 1
        assert rows == len(eng["faces"])
        candidates[label] = dict(fullMeshRows=rows, selectedFiniteHeatFaces=finite)
    candidates["kinney"]["selectedAcceptedFaces"] = int(np.count_nonzero(mesh.heating.kinneyAcreage.acceptedFace[subset]))
    assert candidates["kinney"]["selectedAcceptedFaces"] == 0
    assert all(s["selectedFiniteHeatFaces"] == 0 for s in candidates.values())
    assert assessed["engineering"]["subsetFiniteHeat"] == assessed["research"]["subsetEngineeringHeat"] == 0
    assert assessed["research"]["subsetResearchHeat"] == assessed["research"]["subsetFiniteHeat"] == len(pts)
    controls = list(tables(work / "signed.csv"))
    for k in ("pressure", "heat_flux"):
        assert [r[k] for r in controls] == ["0","-2",""]
        assert [r[k+"_availability"] for r in controls] == ["engineering","prescribed","out-of-domain"]
    assert [r["joint_engineering"] for r in controls] == ["1","0","0"]
    qualified = list(tables(work / "qualified_synthetic.csv"))
    assert [int(r["face_id"]) for r in qualified if int(r["joint_verified"])] == [1,2]
    assert all(r["synthetic"] == "1" for r in controls+qualified)
    assert qualified[1]["heat_flux"] == "-2" and qualified[1]["heat_flux_availability"] == "engineering"
    assert qualified[5]["heat_flux"] == "6" and qualified[5]["heat_flux_availability"] == "out-of-domain"
    report = dict(targetFaces=len(pts), targetAreaM2=math.fsum(eng["faceArea"][i] for i in subset),
                  exports=assessed, archivedEngineeringCandidates=candidates, researchCheckpoint=frame["checkpointIteration"],
                  researchAltitudeRecorded=other["altitudeM"] is not None,
                  newJointEngineeringAreaM2=0, newIndependentVerifiedAreaM2=0, physicalSolverRun=False,
                  inputs={p.relative_to(ROOT).as_posix(): h for p,h in EXPECTED.items()},
                  exportedFiles={n:sha(work/n) for n in ("engineering.csv","research.csv","signed.csv","qualified_synthetic.csv")},
                  elapsedSeconds=time.monotonic()-start)
    with dest.open("x",encoding="utf8") as f:
        json.dump(report, f, indent=2, allow_nan=False)
    print(json.dumps(assessed))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_folder", type=Path)
    parser.add_argument("new_report", type=Path)
    args = parser.parse_args()
    run(args.export_folder, args.new_report)
