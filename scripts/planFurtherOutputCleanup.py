"""核对三个旧研究的四项派生缓存；规划、删除后审计分别执行。"""

from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import subprocess
import sys

from scipy.io import loadmat
from planLegacyOutputCleanup import ROOT, OUT, safe, sha, rec, inventory

DEST = ROOT / "docs/validation/storage_20261029/factors"
BASE = OUT / "v83_validation"
GROUPS = {
    "blockilureal01": ("studyRocketBlockILU0", 1, [
        "steady_graph_packed_evidence.mat", "block_ilu_factor_evidence.mat"]),
    "coupledpatchreal01": ("studyRocketCoupledPatches", 0, [
        "coupled_patch_factor_evidence.mat"]),
    "fullFstrengthreal02": ("studyRocketFullFStrength", 0, [
        "full_F_strength_factor_evidence.mat"]),
}
KEYS = ["sourceFile", "sourceBinding", "sourceIdentity", "packedFile", "rhsFile",
        "graphFile", "inputFile", "baseFile", "assemblyFile"]


def read(path):
    return json.loads(safe(path).read_text(encoding="utf-8-sig"))


def paths(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from paths(child)
    elif isinstance(value, str) and value.lower().endswith(".mat"):
        path = Path(value)
        yield safe(path if path.is_absolute() else ROOT / path)


def plan():
    """模块一：固定删除范围；保留旧Q17、当前Q41、实际结果和历史源码。"""
    assert not DEST.exists(), "已有规划目录，不得覆盖历史记录。"
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    deleted_paths = {BASE / group / name for group, (_, _, names) in GROUPS.items() for name in names}
    protected, source_versions, groups, removed = {}, {}, [], []
    pending_archives = {}
    archived = {}
    for directory in [ROOT / "docs/validation/storage_20261028", ROOT / "docs/validation/storage_20261029"]:
        for file in directory.rglob("*.raw"):
            archived.setdefault(sha(file), file)

    def protect(path, expected=None):
        path = safe(path)
        if path in deleted_paths:
            return
        item = rec(path, expected)
        protected[item["path"]] = item

    # 完整锁定每个原进程的MATLAB源码版本；找不到历史原字节则拒绝清理。
    def source(rel, expected):
        rel = rel.replace("\\", "/")
        key = (rel, expected)
        if key in source_versions:
            return
        current = safe(ROOT / rel)
        file = current if current.is_file() and sha(current) == expected else archived.get(expected)
        if file is None:
            saved = subprocess.run(["git", "show", commit + ":" + rel], cwd=ROOT,
                                   capture_output=True, check=False)
            assert saved.returncode == 0 and hashlib.sha256(saved.stdout).hexdigest() == expected, "缺少匹配的历史源码：" + rel
            target = DEST / "historical_sources" / (current.name + "." + expected[:12] + ".raw")
            pending_archives[target] = saved.stdout
            source_versions[key] = dict(path=target.relative_to(ROOT).as_posix(), bytes=len(saved.stdout),
                                        sha256=expected, restorePath=rel, sourceCommit=commit,
                                        restoredFromGitWithoutChangingWorkingFile=True)
            return
        item = dict(rec(file, expected), restorePath=rel)
        if file == current:
            saved = subprocess.run(["git", "show", commit + ":" + rel], cwd=ROOT,
                                   capture_output=True, check=False)
            if saved.returncode == 0 and hashlib.sha256(saved.stdout).hexdigest() == expected:
                item["sourceCommit"] = commit
            else:
                item["requiresRawArchive"] = True
        source_versions[key] = item

    for group, (constructor, exit_code, names) in GROUPS.items():
        folder = BASE / group
        assert not (folder / "CACHE_PRUNED.json").exists()
        process_file = BASE / (group + "_process.json")
        process = read(process_file)
        assert process["status"] == "terminal" and process["exitCode"] == exit_code
        assert process["sourceStableAfter"] and all(process["sourceStableAfter"].values())
        protect(process_file)
        for rel, expected in process["sourcesBefore"].items():
            if rel.endswith(".m"):
                source(rel, expected)
        artifacts = []
        for name in names:
            file = safe(folder / name)
            assert file.stat().st_nlink == 1, "本次只处理独占的旧缓存。"
            metadata = loadmat(str(file), variable_names=KEYS, simplify_cells=True)
            binding = metadata["sourceBinding"]
            protect(Path(binding["sourceFile"]), binding["sourceFileSHA256"])
            for dep in paths(metadata):
                protect(dep)
            if "packedFileSHA256" in binding:
                protect(Path(binding["packedFile"]), binding["packedFileSHA256"])
            item = rec(file)
            artifacts.append(item)
            removed.append(item)
        for file in folder.iterdir():
            if file.is_file() and file not in deleted_paths:
                protect(file)
        groups.append(dict(folder=folder.relative_to(ROOT).as_posix(), constructorStudy=constructor,
                           originalProcess=rec(process_file), originalExitCode=exit_code,
                           removedArtifacts=artifacts))
        print("VERIFIED " + group, flush=True)

    # F强连接覆盖另依赖原分区；只删除因子，保留分区、原图和数值结果。
    protect(BASE / "spatialreductionpreflight01/spatial_reduction_input.mat")
    for relative in ["v93_validation/continuous01/step01/krylov01/wall_velocity_krylov_evidence.mat",
                     "v83_validation/krylov41/wall_velocity_krylov_evidence.mat",
                     "v93_validation/continuous01/step01/active_direction.mat",
                     "v93_validation/continuous01/step01/active_model.mat",
                     "v93_validation/continuous01/step01/step.json",
                     "v80_validation/domain02/shell_domain_evidence.mat"]:
        protect(OUT / relative)
    # 最新v94的全部本地验证与工作台正在交付收口，逐文件保留，绝不清理。
    for file in (OUT / "v94_validation").rglob("*"):
        if file.is_file():
            protect(file)
    assert len(removed) == 4 and not {x["path"] for x in removed}.intersection(protected)

    DEST.mkdir(parents=True)
    for target, content in pending_archives.items():
        target.parent.mkdir(exist_ok=True)
        target.write_bytes(content)
    for item in source_versions.values():
        if item.pop("requiresRawArchive", False):
            file = ROOT / item["path"]
            target = DEST / "historical_sources" / (file.name + "." + item["sha256"][:12] + ".raw")
            target.parent.mkdir(exist_ok=True)
            if not target.exists():
                target.write_bytes(file.read_bytes())
            item["archivedFile"] = rec(target, item["sha256"])
    before = inventory()
    inv = DEST / "inventory_before.json.gz"
    with inv.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            stream.write(json.dumps(before).encode("utf8"))
    data = dict(kind="four-obsolete-Q17-graph-and-factor-caches", workspace=str(ROOT),
                createdUTC=datetime.now(timezone.utc).isoformat(), sourceCommit=commit,
                deleteFiles=removed, groups=groups, protectedFiles=list(protected.values()),
                sourceFiles=list(source_versions.values()), inventoryBefore=rec(inv),
                plannedBytes=sum(x["bytes"] for x in removed),
                uniqueFileLengthBytesBefore=before["uniqueFileLengthBytes"],
                numericalRebuildExecuted=False, historicalMATFileSHAReproductionClaimed=False,
                currentState="Q41 unchanged; v94 conditional-streamline results retained",
                note="blockilureal01原进程exit1；保留原失败及后续失败证据，不将其改记成功。")
    (DEST / "plan.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf8")
    print(json.dumps(dict(files=4, plannedBytes=data["plannedBytes"],
                          protectedFiles=len(protected), sourceVersions=len(source_versions),
                          uniqueBytesBefore=before["uniqueFileLengthBytes"])), flush=True)


def audit():
    """模块二：全目录独立复核，只允许四项缓存消失及三个清理标记新增。"""
    data = read(DEST / "plan.json")
    receipt = read(DEST / "execution.json")
    assert receipt["status"] == "completed-four-obsolete-Q17-caches-removed"
    assert sha(DEST / "plan.json") == receipt["planSHA256"]
    rec(ROOT / data["inventoryBefore"]["path"], data["inventoryBefore"]["sha256"])
    with gzip.open(DEST / "inventory_before.json.gz", "rb") as stream:
        before = json.load(stream)
    after = inventory()
    missing = set(before["files"]) - set(after["files"])
    assert missing == {x["path"] for x in data["deleteFiles"]}
    assert missing == {x["path"] for x in receipt["deleted"]}
    changed = [p for p in before["files"] if p in after["files"] and before["files"][p] != after["files"][p]]
    document_updates = []
    note = DEST / "readme_update.json"
    if note.exists():
        update = read(note)
        rel = update["after"]["path"]
        assert rel == "outputs/README.md" and update["beforeMetadata"] == before["files"][rel]
        assert update["afterMetadata"] == after["files"][rel]
        rec(ROOT / rel, update["after"]["sha256"])
        document_updates.append(rel)
    assert set(changed) == set(document_updates), changed
    markers = {g["folder"] + "/CACHE_PRUNED.json" for g in data["groups"]}
    assert set(after["files"]) - set(before["files"]) == markers
    for item in data["protectedFiles"] + data["sourceFiles"]:
        rec(ROOT / item["path"], item["sha256"])
        if "archivedFile" in item:
            rec(ROOT / item["archivedFile"]["path"], item["sha256"])
    doc_delta = sum(after["files"][p]["bytes"] - before["files"][p]["bytes"] for p in document_updates)
    expected = before["uniqueFileLengthBytes"] - data["plannedBytes"] + sum(after["files"][p]["bytes"] for p in markers) + doc_delta
    assert after["uniqueFileLengthBytes"] == expected
    result = dict(status="passed-independent-output-inventory-and-preserved-SHA-audit",
                  deletedFiles=sorted(missing), unchangedOriginalFiles=before["fileCount"]-4-len(document_updates),
                  intentionalDocumentationUpdates=document_updates,
                  protectedFiles=len(data["protectedFiles"]), sourceVersions=len(data["sourceFiles"]),
                  uniqueFileLengthBytesBefore=before["uniqueFileLengthBytes"],
                  uniqueFileLengthBytesAfter=after["uniqueFileLengthBytes"],
                  actualFreeIncreaseBytes=receipt["actualFreeIncreaseBytes"],
                  noNumericalComputationExecuted=True, currentQ41AndAllV94OutputsRetained=True)
    report_file = DEST / ("final_postcheck.json" if document_updates else "postcheck.json")
    assert not report_file.exists(), "不得覆盖独立复核原记录。"
    report_file.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    {"plan": plan, "audit": audit}[sys.argv[1]]()
