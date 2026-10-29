"""固定范围的旧输出缓存清理：先保存来源/目录，删除后独立核对。"""

from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import os
import subprocess
import sys

from scipy.io import loadmat

from planLegacyOutputCleanup import ROOT, OUT, safe, sha, rec, inventory

DEST = ROOT / "docs/validation/storage_20261029/followup"
OLD = OUT / "v83_validation"
STEP = OUT / "v93_validation/continuous01/step01"
FACTORS = {
    "patchhierarchyreal01": ("patch_hierarchy_factor_evidence.mat", "studyRocketPatchHierarchy"),
    "blockfillreal01": ("block_fill1_factor_evidence.mat", "studyRocketBlockFillLevel1"),
    "overlapFreal01": ("restricted_overlap_factor_evidence.mat", "studyRocketRestrictedOverlap"),
    "floworderfactor01": ("flow_ordered_factor_evidence.mat", "studyRocketFlowOrderingFactor"),
    "spatialreductionreal01": ("spatial_coarse_factor_evidence.mat", "studyRocketSpatialReduction"),
}
MAT_KEYS = ["sourceFile", "sourceBinding", "sourceIdentity", "factorIdentity", "packedFile",
            "rhsFile", "anchorFile", "geometryFile", "inputFile", "factorFile", "graphFile",
            "topologyFile", "baseFile", "assemblyFile"]


def read(path):
    return json.loads(safe(path).read_text(encoding="utf-8-sig"))


def terminal(path, code=0):
    value = read(path)
    assert value["status"] == "terminal" and value["exitCode"] == code, str(path)
    assert all(value["sourceStableAfter"].values()), str(path)
    return value


def mat_paths(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from mat_paths(child)
    elif isinstance(value, str) and value.lower().endswith(".mat"):
        path = Path(value)
        yield safe(path if path.is_absolute() else ROOT / path)


def plan():
    """模块一：只读核对五项旧因子与六项Q40矩阵，保存所需历史源码。"""
    assert not (DEST / "plan.json").exists(), "已有清单，不能覆盖或扩大范围。"
    DEST.mkdir(parents=True, exist_ok=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    protected, source_versions, groups, deleted = {}, {}, [], []

    def protect(path, expected=None):
        item = rec(path, expected)
        protected[item["path"]] = item
        return item

    # 旧构造源码已经在此前清理中保存。按原字节SHA寻找，缺失即拒绝。
    archived = {}
    for parent in [ROOT / "docs/validation/storage_20261028", ROOT / "docs/validation/storage_20261029"]:
        for path in parent.rglob("*.raw"):
            archived.setdefault(sha(path), path)

    def bind_source(rel, expected):
        rel = rel.replace("\\", "/")
        key = (rel, expected)
        if key in source_versions:
            return
        original = safe(ROOT / rel)
        path = original if original.is_file() and sha(original) == expected else archived.get(expected)
        assert path is not None, "找不到原始构造源码：" + rel
        item = dict(rec(path, expected), restorePath=rel)
        if path == original:
            result = subprocess.run(["git", "show", commit + ":" + rel], cwd=ROOT,
                                    capture_output=True, check=False)
            if result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == expected:
                item["sourceCommit"] = commit
            else:
                target = DEST / "historical_sources" / (original.name + "." + expected[:12] + ".raw")
                target.parent.mkdir(exist_ok=True)
                if not target.exists():
                    target.write_bytes(path.read_bytes())
                item["archivedFile"] = rec(target, expected)
        source_versions[key] = item

    for name, (filename, study) in FACTORS.items():
        folder = OLD / name
        assert not (folder / "CACHE_PRUNED.json").exists()
        process_file = OLD / (name + "_process.json")
        process = terminal(process_file)
        process_record = protect(process_file)
        for rel, expected in process["sourcesBefore"].items():
            if rel.endswith(".m"):
                bind_source(rel, expected)
        path = folder / filename
        assert path.stat().st_nlink == 1
        metadata = loadmat(path, variable_names=MAT_KEYS, simplify_cells=True)
        for dep in mat_paths({k: v for k, v in metadata.items() if not k.startswith("__")}):
            protect(dep)
        binding = metadata["sourceBinding"]
        protect(Path(binding["sourceFile"]), binding["sourceFileSHA256"])
        # overlap构造入口还显式核对旧强度因子；该依赖保留。
        if name == "overlapFreal01":
            topology = loadmat(metadata["topologyFile"], variable_names=["sourceStrengthFactorSHA256"],
                               simplify_cells=True)
            protect(OLD / "fullFstrengthreal02/full_F_strength_factor_evidence.mat",
                    topology["sourceStrengthFactorSHA256"])
        item = rec(path)
        deleted.append(item)
        for other in folder.iterdir():
            if other.is_file() and other != path:
                protect(other)
        groups.append(dict(folder=folder.relative_to(ROOT).as_posix(), process=process_record,
                           constructorStudy=study, removedArtifacts=[item],
                           rebuildScope="factor-construction-only-from-preserved-inputs-in-new-folder"))

    # Q41接受、原失败/恢复和独立数值审计均终态；本轮不重新求解。
    for filename, code in [("coordinatefull93_process.json", 1),
                           ("coordinaterecover93_process.json", 0), ("coordinateview93_process.json", 0)]:
        terminal(OUT / "v93_validation" / filename, code)
        protect(OUT / "v93_validation" / filename)
    phase = read(STEP.parent / "active_continuation.json")
    assert phase["lastAcceptedFolder"] == str(STEP / "krylov01")
    step = read(STEP / "step.json")
    assert step["main"]["acceptedNonlinearIterations"] == 1
    assert [step[k] for k in ["newPreparationActualNS", "newProbeActualNS", "newMainActualNS"]] == [1, 2, 4]
    binding = read(OUT / "v93_validation/recovery_binding.json")
    protect(OUT / "v93_validation/recovery_binding.json")
    for key in ["failure", "context", "prefix", "preview", "independentProjectorProof"]:
        protect(ROOT / binding[key]["file"], binding[key]["sha256"])
    for rel in ["continuous01/independent_temperature_active.json", "independent_cell_interface_continuation.json",
                "row_real02/independent_row_conditioning.json"]:
        path = OUT / "v93_validation" / rel
        assert read(path)["status"].startswith("passed-")
        protect(path)
    protect(STEP / "krylov01/wall_velocity_krylov_evidence.mat",
            "dc9fea4e454f137fa4e4ead736bed37ffb5407eed2cbc6b5f5012ca9e4ffd1ce")
    protect(OLD / "krylov41/wall_velocity_krylov_evidence.mat",
            "dc9fea4e454f137fa4e4ead736bed37ffb5407eed2cbc6b5f5012ca9e4ffd1ce")
    for rel in ["active_direction.mat", "active_model.mat", "step.json"]:
        protect(STEP / rel)
    current = STEP / "current"
    preparation_file = current / "peak_system.json"
    preparation = read(preparation_file)
    protect(preparation_file)
    assert preparation["sourceBinding"]["sourceFileSHA256"] == "3ba7f649bb86fc211a1bd30aad9fca31c296e6d321cdb14c039dd4316b26ddb1"
    protect(Path(preparation["sourceBinding"]["sourceFile"]), preparation["sourceBinding"]["sourceFileSHA256"])
    protect(Path(preparation["geometryBinding"]["file"]), preparation["geometryBinding"]["sha256"])
    for identity in [preparation["sourceIdentity"], preparation["methodIdentity"]]:
        for name, expected in identity.items():
            paths = [ROOT / base / (name + ".m") for base in ["internalFlow/cutCell", "internalFlow", "tests", "."]]
            paths = [p for p in paths if p.is_file()]
            matching = [p for p in paths if sha(p) == expected]
            assert matching, "准备模块身份不符：" + name
            bind_source(matching[0].relative_to(ROOT).as_posix(), expected)
    six = []
    assert set(preparation["artifacts"]) == {"cache", "operator", "packed", "factor", "scale", "metric"}
    for value in preparation["artifacts"].values():
        path = safe(value["file"])
        assert path.parent == current and path.stat().st_nlink == 1
        six.append(rec(path, value["sha256"]))
    deleted += six
    groups.append(dict(folder=current.relative_to(ROOT).as_posix(), constructorStudy="prepareRocketCurrentPeakSystem",
                       preparation=rec(preparation_file), removedArtifacts=six,
                       rebuildScope="rebuild-Q40-preparation-only-in-new-folder-next-advancement-requires-new-Q41-system"))
    assert len(deleted) == 11
    assert not {x["path"] for x in deleted}.intersection(protected)
    before = inventory()
    inv = DEST / "inventory_before.json.gz"
    with inv.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            stream.write(json.dumps(before).encode())
    result = dict(kind="five-legacy-factors-and-six-obsolete-Q40-preparation-caches", workspace=str(ROOT),
                  createdUTC=datetime.now(timezone.utc).isoformat(), sourceCommit=commit,
                  deleteFiles=deleted, groups=groups, protectedFiles=list(protected.values()),
                  sourceFiles=list(source_versions.values()), inventoryBefore=rec(inv),
                  plannedBytes=sum(x["bytes"] for x in deleted),
                  uniqueFileLengthBytesBefore=before["uniqueFileLengthBytes"],
                  numericalRebuildExecuted=False, historicalMATFileSHAReproductionClaimed=False,
                  v93PublicationStatus="numerics-tested-UI-delivery-closeout-pending-not-published-by-storage-maintenance")
    (DEST / "plan.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    summary = {k: result[k] for k in ["kind", "plannedBytes", "uniqueFileLengthBytesBefore"]}
    summary.update(files=len(deleted), protected=len(protected), sourceVersions=len(source_versions))
    print(json.dumps(summary))


def audit():
    """模块二：全目录独立比较，只允许固定十一项消失，关键保留文件再验SHA。"""
    data = read(DEST / "plan.json")
    receipt = read(DEST / "execution.json")
    assert receipt["status"] == "completed-eleven-explicit-obsolete-caches-removed"
    with gzip.open(DEST / "inventory_before.json.gz", "rb") as stream:
        before = json.load(stream)
    after = inventory()
    removed = set(before["files"]) - set(after["files"])
    assert removed == {x["path"] for x in data["deleteFiles"]}
    changed = [p for p in before["files"] if p in after["files"] and before["files"][p] != after["files"][p]]
    assert not changed, changed
    for item in data["protectedFiles"] + data["sourceFiles"]:
        rec(ROOT / item["path"], item["sha256"])
        if "archivedFile" in item:
            rec(ROOT / item["archivedFile"]["path"], item["sha256"])
    result = dict(status="passed-independent-full-output-inventory-and-protected-SHA-audit",
                  removed=sorted(removed), unchangedOriginalFiles=before["fileCount"]-len(removed),
                  uniqueFileLengthBytesBefore=before["uniqueFileLengthBytes"],
                  uniqueFileLengthBytesAfter=after["uniqueFileLengthBytes"],
                  actualFreeIncreaseBytes=receipt["actualFreeIncreaseBytes"],
                  newOutputFiles=sorted(set(after["files"])-set(before["files"])),
                  noNumericalComputationExecuted=True, currentQ41AndOriginalEvidenceRetained=True)
    (DEST / "postcheck.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps(result))


if __name__ == "__main__":
    {"plan": plan, "audit": audit}[sys.argv[1]]()
