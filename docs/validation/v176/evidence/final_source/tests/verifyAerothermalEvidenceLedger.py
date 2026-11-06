"""独立核验逐原面证据台账；只读小型表面缓存，不运行任何流场求解。

模块一以坐标叉积和集合重算面积，模块二检查合成资格真值，模块三
将真实面值及CSV逐行绑定到既存MAT。MATLAB JSON哈希不冒称跨语言重算。
"""
from __future__ import print_function
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
SOURCE = {
    "engineering": ("outputs/v104_validation/baseline02/current_case.mat",
                    "b4889e42a28131c253a56f06472a58c9d2cb290d7c9937935fbb4b5d7349488f", "current"),
    "research": ("outputs/v167_validation/main01/result.mat",
                 "ac7ffa2ad8eee390f6f0766c5de99879cb912e129caa865a02aa721373cb2744", "r"),
}
HEADERS = ["face_id", "region", "area_m2", "field_key", "value", "unit",
           "availability_code", "availability", "qualification_code", "qualification",
           "method", "scope", "reason", "evidence_ids"]
AVAILABILITY = ["未知", "工程计算", "研究计算", "给定条件", "超出适用域"]
QUALIFICATION = ["未验证", "验证失败", "验证通过"]
CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    with Path(path).open("r", encoding="utf-8-sig") as stream:
        return json.load(stream)


def seq(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def fields(value):
    return {item["key"]: item for item in seq(value["fields"])}


def equal_number(a, b, label, tolerance=3e-12):
    if a is None or b is None:
        check(a is None and b is None, label + " unknown value")
    else:
        check(math.isfinite(float(a)) and math.isfinite(float(b)), label + " finite")
        check(abs(float(a) - float(b)) <= tolerance * max(abs(float(a)), abs(float(b)), 1e-300), label)


def equal_vector(a, b, label, tolerance=3e-12):
    a, b = seq(a), seq(b)
    check(len(a) == len(b), label + " length")
    for x, y in zip(a, b):
        equal_number(x, y, label, tolerance)


def same_nested(a, b, label):
    if isinstance(a, list):
        check(isinstance(b, list) and len(a) == len(b), label + " shape")
        for x, y in zip(a, b):
            same_nested(x, y, label)
    elif isinstance(a, (float, int)) and not isinstance(a, bool):
        equal_number(a, b, label)
    else:
        check(a == b, label)


def triangle_areas(data):
    vertices, faces = data["vertices"], data["faces"]
    result = []
    for face in faces:
        check(len(face) == 3 and len(set(face)) == 3, "original triangular face")
        check(all(isinstance(x, (int, float)) and int(x) == x and 0 <= x < len(vertices) for x in face), "zero-based original node IDs")
        a, b, c = [vertices[int(x)] for x in face]
        u = [b[j] - a[j] for j in range(3)]
        v = [c[j] - a[j] for j in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        area = math.sqrt(math.fsum(x * x for x in cross)) / 2
        check(math.isfinite(area) and area > 0, "positive independently calculated area")
        result.append(area)
    equal_vector(data["faceArea"], result, "independent triangle area", 2e-10)
    return result


def area_summary(areas, availability, qualification, ids):
    total = math.fsum(areas[i] for i in ids)
    out = {"totalArea": total}
    for name, code in [("unknownArea", 0), ("engineeringArea", 1), ("researchArea", 2), ("prescribedArea", 3), ("outOfDomainArea", 4)]:
        out[name] = math.fsum(areas[i] for i in ids if availability[i] == code)
        out[name + "Fraction"] = out[name] / total
    for name, code in [("qualificationUnknownArea", 0), ("qualificationFailedArea", 1), ("qualifiedArea", 2)]:
        out[name] = math.fsum(areas[i] for i in ids if qualification[i] == code)
        out[name + "Fraction"] = out[name] / total
    return out


def summary_matches(actual, expected, label):
    for name, value in expected.items():
        check(name in actual, label + " missing " + name)
        equal_number(actual[name], value, label + " " + name, 4e-12)


def expected_availability(data, field):
    origin = field.get("evidenceOrigin") or data.get("evidenceIdentity", {}).get("origin", "engineering")
    status = seq(field["status"])
    outside = seq(field.get("outOfDomainFace", []))
    if not outside:
        outside = [False] * len(status)
    check(len(outside) == len(status), "outside mask length")
    return [4 if outside[i] else 3 if s == 2 else (2 if origin == "research" else 1) if s == 1 else 0 for i, s in enumerate(status)]


def qualification_for_fixture(name, field, availability):
    q = [1 if a == 4 else 0 for a in availability]
    if name in ("qualifiedLedger", "duplicateLedger", "implicitLedger"):
        if field == "pressure":
            for i in (0, 1, 2):
                q[i] = 2
        elif field == "heatFlux":
            q[0], q[1], q[2] = 2, 2, 1
    return q


def joint_summary(areas, calculated, ids):
    selected = set(ids)
    p, q = calculated.get("pressure"), calculated.get("heatFlux")
    sets = {"Available": set(), "Qualified": set(), "Engineering": set()}
    if p and q:
        pa = {i for i, a in enumerate(p["a"]) if a in (1, 2)}
        qa = {i for i, a in enumerate(q["a"]) if a in (1, 2)}
        sets["Available"] = pa & qa & selected
        sets["Qualified"] = sets["Available"] & {i for i, a in enumerate(p["q"]) if a == 2} & {i for i, a in enumerate(q["q"]) if a == 2}
        sets["Engineering"] = selected & {i for i, a in enumerate(p["a"]) if a == 1} & {i for i, a in enumerate(q["a"]) if a == 1}
    result = {"totalArea": math.fsum(areas[i] for i in ids)}
    for name, face_ids in sets.items():
        prefix = "pressureHeatFlux" + name
        result[prefix + "Area"] = math.fsum(areas[i] for i in face_ids)
        result[prefix + "Fraction"] = result[prefix + "Area"] / result["totalArea"]
        result[prefix + "FaceIDs"] = sorted(face_ids)
    return result


def joint_matches(actual, expected, label):
    for name, value in expected.items():
        if name.endswith("FaceIDs"):
            if name in actual:
                check(seq(actual[name]) == value, label + " " + name)
        else:
            equal_number(actual[name], value, label + " " + name, 4e-12)


def verify_ledger(data, ledger, name="real"):
    areas = triangle_areas(data)
    n = len(areas)
    check(len(seq(data["region"])) == n, "all original region tags")
    for key in ("vertices", "faces", "faceArea", "region"):
        if key in ledger:
            same_nested(ledger[key], data[key], "ledger preserves " + key)
    equal_number(ledger.get("totalArea", math.fsum(areas)), math.fsum(areas), "full original area")
    check(ledger.get("faceCount", n) == n, "full original face count")
    original, actual = fields(data), fields(ledger)
    check(set(original) == set(actual), "all variables retained")
    calculated = {}
    for key, f in original.items():
        lf = actual[key]
        values, status = seq(f["values"]), seq(f["status"])
        check(len(values) == n and len(status) == n, "per-face source lengths")
        check(all((v is None) == (s == 0) for v, s in zip(values, status)), "legacy status and unknown agreement")
        for member in ("values", "status"):
            if member in lf:
                same_nested(lf[member], f[member], "ledger raw " + member)
        if "legacyStatus" in lf:
            same_nested(lf["legacyStatus"], f["status"], "legacyStatus preserved")
        a = expected_availability(data, f)
        q = qualification_for_fixture(name, key, a)
        check(seq(lf["availability"]) == a, key + " availability independent")
        check(seq(lf["qualification"]) == q, key + " qualification independent")
        summary = area_summary(areas, a, q, list(range(n)))
        summary_matches(lf["areaSummary"], summary, key + " global")
        calculated[key] = {"a": a, "q": q, "summary": summary}
    region_tags = seq(data["region"])
    independent_regions = {tag: [i for i, value in enumerate(region_tags) if value == tag] for tag in set(region_tags)}
    seen = []
    region_result = {}
    for region in seq(ledger["regions"]):
        name_r = region["name"]
        check(name_r in independent_regions, "original region name")
        ids = independent_regions[name_r]
        if "faceIDs" in region:
            check(seq(region["faceIDs"]) == ids, "original disjoint region face IDs")
        if "faceCount" in region:
            check(region["faceCount"] == len(ids), "region face count")
        seen.extend(ids)
        equal_number(region["area"], math.fsum(areas[i] for i in ids), "original region area", 4e-12)
        region_result[name_r] = {"area": math.fsum(areas[i] for i in ids), "faceCount": len(ids), "fields": {}}
        for key, item in calculated.items():
            summary = area_summary(areas, item["a"], item["q"], ids)
            summary_matches(region["fields"][key], summary, key + " region " + name_r)
            region_result[name_r]["fields"][key] = summary
        j = joint_summary(areas, calculated, ids)
        joint_matches(region["joint"], j, "region joint")
        region_result[name_r]["joint"] = j
    check(sorted(seen) == list(range(n)), "region union is exactly all original faces once")
    joint = joint_summary(areas, calculated, list(range(n)))
    joint_matches(ledger["joint"], joint, "global joint")
    return {"area": math.fsum(areas), "faceCount": n, "fields": calculated, "regions": region_result, "joint": joint}


def verify_proofs(fixture):
    identity = fixture["initialLedger"]["identity"]
    lf = fields(fixture["initialLedger"])
    for record in seq(fixture["records"]):
        path = Path(record["proofFile"])
        check(path.is_absolute() and path.stat().st_size <= 1024 ** 2, "actual small proof file")
        check(sha(path) == record["proofSHA256"], "proof exact SHA")
        proof = read(path)
        check(proof["evidenceID"] == record["id"], "proof evidence ID")
        check(proof["bindings"] == record["bindings"], "record/proof bindings")
        for key in ("geometrySHA256", "conditionSHA256", "sourceSHA256"):
            check(record["bindings"][key] == identity[key], "bound " + key)
        check(record["bindings"]["valueSHA256"] == lf[record["variable"]]["valueSHA256"], "bound variable hash")
        check(seq(proof["faceIDs"]) == seq(record["faceIDs"]), "bound original face set")
        decisions = []
        for c in seq(proof["criteria"]):
            x, y = c["observed"], c["threshold"]
            op = c["operator"]
            check(op in ("<=", ">=", "<", ">", "==", "abs<="), "known proof comparison")
            comparisons = {"<=": x <= y, ">=": x >= y, "<": x < y, ">": x > y, "==": x == y, "abs<=": abs(x) <= y}
            decisions.append(comparisons[op])
        check(all(decisions) == (not record["id"].endswith("fail")), "independent numeric proof result")


def verify_fixture(path):
    f = read(path)
    check(f["schema"] == "RAH.EvidenceLedger.IndependentFixture.v1" and f["syntheticOnly"], "explicit synthetic fixture")
    check(f["newQualifiedPhysicalArea"] == 0, "synthetic qualification is not physical qualification")
    equal_vector(triangle_areas(f["data"]), [1, 2, 4, 8, 16, 32], "exact unequal areas", 0)
    for ledger_name, data_name in [("initialLedger", "data"), ("qualifiedLedger", "data"), ("duplicateLedger", "data"),
                                   ("conflictedLedger", "data"), ("implicitLedger", "implicitData"),
                                   ("nanOutsideLedger", "nanOutsideData"), ("protectedLedger", "data"),
                                   ("researchLedger", "researchData"), ("legacyLedger", "legacyData")]:
        verify_ledger(f[data_name], f[ledger_name], ledger_name)
    check(f["duplicateLedger"]["duplicateEvidenceCount"] == 1, "one identical proof deduplicated")
    check(f["duplicateLedger"]["acceptedEvidenceCount"] == 4, "unique proofs counted once")
    check(f["conflictedLedger"]["rejectedEvidenceCount"] == 2, "same ID conflict rejects whole group")
    check(f["qualifiedLedger"]["joint"]["pressureHeatFluxQualifiedArea"] == 3, "intersection 3, not marginal mean 5")
    rejection_names = {item["name"] for item in seq(f["rejections"])}
    required = {"same-values-changed-method", "same-values-changed-status", "same-values-changed-field-origin",
                "stale-condition-data", "stale-source-data", "original-face-order-changed", "generic-passed-no-criteria",
                "synthetic-proof-on-real-data", "heat-without-nearwall-rule"}
    check(required <= rejection_names, "stale provenance and invalid proof cases executed")
    verify_proofs(f)
    return {"faces": 6, "area": 63, "qualifiedSyntheticPressureArea": 7, "qualifiedSyntheticHeatArea": 3,
            "qualifiedSyntheticIntersectionArea": 3, "rejectionCases": len(rejection_names), "newPhysicalArea": 0}


def mat_array(value):
    import numpy as np
    a = np.asarray(value)
    if a.ndim == 0:
        return a.item()
    return a.tolist()


def mat_values(value):
    import numpy as np
    return [None if not math.isfinite(float(x)) else float(x) for x in np.asarray(value).reshape(-1)]


def verify_real(data, origin):
    # 只载入已映射表面MAT；本函数没有体场文件入口。
    from scipy.io import loadmat
    relative, expected_sha, variable = SOURCE[origin]
    path = ROOT / relative
    before = sha(path)
    check(before == expected_sha, "fixed historical surface MAT " + origin)
    result = loadmat(str(path), variable_names=[variable], squeeze_me=True, struct_as_record=False)[variable]
    same_nested(data["vertices"], mat_array(result.V), origin + " original vertices")
    same_nested(data["faces"], (result.F - 1).tolist(), origin + " original face order")
    source_fields = fields(data)
    if origin == "engineering":
        values = {"pressure": mat_values(result.pressure.staticPressure / 1000),
                  "heatFlux": mat_values(result.heating.heatFlux / 1000)}
        units = {"pressure": "kPa", "heatFlux": "kW/m^2"}
        states = {"pressure": mat_values(result.pressure.fieldStatus),
                  "heatFlux": mat_values(result.heating.eligibleFace)}
    else:
        values = {"pressure": mat_values(result.pressureRaw), "heatFlux": mat_values(result.heatFluxRaw)}
        units = {"pressure": "Pa", "heatFlux": "W/m^2"}
        states = {key: [0 if value is None else 1 for value in raw] for key, raw in values.items()}
    for key, original in values.items():
        # 老显示器已将旧status=0拒绝值显示为NaN；只核对该既有显示语义。
        status = seq(source_fields[key]["status"])
        expected_status = [0 if original[i] is None else state for i, state in enumerate(states[key])]
        check(status == expected_status, origin + " historical accepted status " + key)
        displayed = [None if status[i] == 0 else value for i, value in enumerate(original)]
        equal_vector(source_fields[key]["values"], displayed, origin + " historical " + key, 8e-14)
        check(source_fields[key]["unit"].replace("²", "^2") == units[key], origin + " original display unit")
    ledger = data["evidenceLedger"]
    verified = verify_ledger(data, ledger)
    check(verified["faceCount"] == 19752, "actual rocket original 19752 faces")
    equal_number(verified["area"], 1.9646292565715551, "actual rocket fixed area", 2e-12)
    check(verified["joint"]["pressureHeatFluxQualifiedArea"] == 0, "no new physical qualification")
    if origin == "research":
        for key, item in verified["fields"].items():
            check(item["summary"]["engineeringArea"] == 0, "research never counted as engineering " + key)
        check(all(v is None for v in seq(source_fields["wallSolved"]["values"])), "unknown solid wall temperature stays unknown")
    check(sha(path) == before, "historical MAT unchanged after read")
    return verified, {"path": str(path), "sha256": before, "bytes": path.stat().st_size}


def decoded_ids(value):
    if value in (None, "", []):
        return []
    if isinstance(value, str):
        try:
            return seq(json.loads(value))
        except ValueError:
            return [value]
    return seq(value)


def reason_at(field, i):
    if "reason" in field:
        return seq(field["reason"])[i]
    return seq(field["reasonNames"])[int(seq(field["reasonIndex"])[i])]


def verify_csv(path, data, recomputed):
    sidecar = Path(path).with_suffix(".json")
    meta = read(sidecar)
    check(meta["kind"] == "aerothermal-evidence-export", "new evidence CSV identity")
    check(meta["csvSHA256"] == sha(path), "CSV SHA matches sidecar")
    check(meta["identity"] == data["evidenceLedger"]["identity"], "CSV same frozen source")
    check(meta["rows"] == recomputed["faceCount"], "CSV exports all original faces")
    key = meta["fieldKey"]
    check(meta["viewState"]["fieldKey"] == key, "CSV current selected field")
    original = fields(data)[key]
    lf = fields(data["evidenceLedger"])[key]
    methods = seq(original["methodNames"])
    method_index = seq(original["methodIndex"])
    values = seq(original["values"])
    scopes = seq(original.get("scopeNames", []))
    scope_index = seq(original.get("scopeIndex", []))
    a = recomputed["fields"][key]["a"]
    q = recomputed["fields"][key]["q"]
    with Path(path).open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        check(reader.fieldnames == HEADERS, "exact evidence CSV schema")
        count = 0
        for i, row in enumerate(reader):
            check(i < len(values), "no excess CSV faces")
            check(int(row["face_id"]) == i + 1, "ordered one-based original CSV face")
            check(row["region"] == seq(data["region"])[i] and row["field_key"] == key, "CSV face region and variable")
            equal_number(float(row["area_m2"]), seq(data["faceArea"])[i], "CSV original face area", 1e-15)
            parsed = None if row["value"] == "" else float(row["value"])
            equal_number(parsed, values[i], "CSV original value including zero/negative/outside", 1e-15)
            check(row["unit"] == original["unit"], "CSV unit")
            check(int(row["availability_code"]) == a[i] and row["availability"] == AVAILABILITY[a[i]], "CSV availability")
            check(int(row["qualification_code"]) == q[i] and row["qualification"] == QUALIFICATION[q[i]], "CSV qualification")
            check(row["method"] == methods[int(method_index[i])], "CSV actual method")
            scope = scopes[int(scope_index[i])] if scopes else original["scope"]
            check(row["scope"] == scope, "CSV actual method scope")
            check(row["reason"] == reason_at(lf, i), "CSV qualification reason")
            check(decoded_ids(row["evidence_ids"]) == decoded_ids(seq(lf["evidenceIDs"])[i]), "CSV evidence IDs")
            count += 1
    check(count == recomputed["faceCount"], "CSV original faces exactly once")
    compact = meta["evidenceLedger"]
    check(compact["identity"] == data["evidenceLedger"]["identity"], "compact ledger identity")
    for key_f, item in fields(compact).items():
        summary_matches(item["areaSummary"], recomputed["fields"][key_f]["summary"], "CSV global " + key_f)
        check(item["valueSHA256"] == fields(data["evidenceLedger"])[key_f]["valueSHA256"], "CSV variable binding")
    joint_matches(compact["joint"], recomputed["joint"], "CSV p/q intersection")
    for r in seq(compact["regions"]):
        expected = recomputed["regions"][r["name"]]
        equal_number(r["area"], expected["area"], "CSV region area", 4e-12)
        check(r["faceCount"] == expected["faceCount"], "CSV region face count")
        for key_f, summary in r["fields"].items():
            summary_matches(summary, expected["fields"][key_f], "CSV region " + key_f)
        joint_matches(r["joint"], expected["joint"], "CSV region intersection")
    return {"path": str(path), "sha256": meta["csvSHA256"], "rows": count, "fieldKey": key,
            "selectedRegion": meta["viewState"].get("regionName", ""), "allOriginalFaces": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--engineering-json", type=Path, required=True)
    parser.add_argument("--research-json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(not args.output.exists(), "refuse to overwrite evidence output")
    args.output.mkdir(parents=True)
    shutil.copyfile(str(Path(__file__).resolve()), str(args.output / "executed_verify.py"))
    started = time.perf_counter()
    inputs = [args.fixture, args.engineering_json, args.research_json] + args.csv + [p.with_suffix(".json") for p in args.csv]
    inputs += [ROOT / item[0] for item in SOURCE.values()]
    before = [{"path": str(p.resolve()), "bytes": p.stat().st_size, "sha256": sha(p)} for p in inputs]
    report = {"passed": False, "newCFDTasks": 0, "newPhysicalQualifiedArea": 0,
              "qualificationDoesNotUseMATLABSummaryAsOracle": True,
              "crossLanguageMATLABJSONHashReproductionClaimed": False}
    try:
        report["fixture"] = verify_fixture(args.fixture)
        datasets, calculations = {}, {}
        for origin, path in [("engineering", args.engineering_json), ("research", args.research_json)]:
            data = read(path)
            verified, provenance = verify_real(data, origin)
            datasets[origin], calculations[origin] = data, verified
            report[origin] = {"faces": verified["faceCount"], "area": verified["area"],
                              "fields": {key: item["summary"] for key, item in verified["fields"].items()},
                              "joint": {key: val for key, val in verified["joint"].items() if not key.endswith("FaceIDs")},
                              "source": provenance}
        report["exports"] = []
        for path in args.csv:
            meta = read(path.with_suffix(".json"))
            origin = meta["identity"]["origin"]
            check(origin in datasets, "export known source origin")
            report["exports"].append(verify_csv(path, datasets[origin], calculations[origin]))
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
    finally:
        stable = all(sha(Path(item["path"])) == item["sha256"] for item in before)
        report["inputsStable"] = stable
        report["passed"] = report["passed"] and stable
        report["checks"] = CHECKS
        report["elapsedSeconds"] = time.perf_counter() - started
        (args.output / "input_manifest.json").write_text(json.dumps(before, ensure_ascii=False, indent=2), encoding="utf-8")
        (args.output / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps({key: report.get(key) for key in ("passed", "checks", "elapsedSeconds", "error", "inputsStable")}, ensure_ascii=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
