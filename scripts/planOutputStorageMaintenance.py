# -*- coding: utf-8 -*-
"""输出维护：回收三项旧构造缓存，透明压缩文本，完整保留科学结果和原路径。"""

import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from scipy.io import loadmat, whosmat

from planLegacyOutputCleanup import ROOT, OUT, inventory, rec, safe, sha

DEST = ROOT / "docs/validation/storage_20261029/text_and_derived"
TARGETS = {
    "spatialreductionreal01/spatial_reduction_assembly_evidence.mat": (
        "spatialreductionreal01", "studyRocketSpatialReduction.m",
        {"reduction", "inputFile", "factorFile", "packedFile", "rhsFile", "sourceBinding", "sourceIdentity"}),
    "fullFreal01/full_F_graph_evidence.mat": (
        "fullFreal01", "studyRocketFullFGraph.m",
        {"graph", "extractionReport", "packedFile", "assemblyFile", "sourceBinding", "sourceIdentity"}),
    "momentumlinereal03/momentum_line_factor_evidence.mat": (
        "momentumlinereal03", "studyRocketMomentumLines.m",
        {"factor", "assemblyFile", "rhsFile", "sourceBinding", "sourceIdentity"}),
}
TEXT_EXTENSIONS = {".csv", ".js", ".json", ".html", ".txt", ".log", ".xml", ".svg", ".m", ".py", ".md", ".css"}


def read(path):
    return json.loads(safe(path).read_text(encoding="utf-8-sig"))


def allocated(path):
    """模块一：读取NTFS实际分配长度；每个文件ID只计一次。"""
    function = ctypes.WinDLL("kernel32", use_last_error=True).GetCompressedFileSizeW
    function.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
    function.restype = wintypes.DWORD
    high = wintypes.DWORD(0)
    ctypes.set_last_error(0)
    low = function(str(safe(path)), ctypes.byref(high))
    error = ctypes.get_last_error()
    assert low != 0xFFFFFFFF or error == 0, (str(path), error)
    return (high.value << 32) | low


def plan():
    """模块二：固定删除范围、历史来源、保留依赖和压缩前的字节身份。"""
    assert not DEST.exists(), "已有本轮记录，不得覆盖。"
    before = inventory()
    print("INVENTORY files=" + str(before["fileCount"]), flush=True)
    aliases = {}
    for relative, item in before["files"].items():
        aliases.setdefault((item["device"], item["fileID"]), []).append(relative)
    groups, deleted, protected, versions, archives = [], [], {}, {}, {}
    # 查找已保留的历史MATLAB源码，不重新复制整个工程。
    source_index = {}
    for base in (ROOT / "docs/validation/storage_20261028", ROOT / "docs/validation/storage_20261029"):
        for path in base.rglob("*.raw"):
            source_index.setdefault(sha(path), path)
    by_name = {}
    for relative in before["files"]:
        if relative.endswith(".m"):
            by_name.setdefault(Path(relative).name, []).append(ROOT / relative)

    def protect(path, expected=None):
        path = Path(path)
        item = rec(path if path.is_absolute() else ROOT / path, expected)
        protected[item["path"]] = item
        return item

    def source(relative, expected):
        relative = relative.replace("\\", "/")
        key = (relative, expected)
        if key in versions:
            return
        current = safe(ROOT / relative)
        candidate = current if current.is_file() and sha(current) == expected else source_index.get(expected)
        if candidate is None:
            candidate = next((p for p in by_name.get(current.name, []) if sha(p) == expected), None)
        assert candidate is not None, "缺少准确历史源码：" + relative
        if candidate == current or str(candidate).startswith(str(ROOT / "docs")):
            versions[key] = dict(rec(candidate, expected), restorePath=relative)
        else:
            archive = DEST / "historical_sources" / (current.name + "." + expected[:12] + ".raw")
            content = candidate.read_bytes()
            archives[archive] = content
            versions[key] = dict(path=archive.relative_to(ROOT).as_posix(), bytes=len(content),
                                 sha256=expected, restorePath=relative, preservedFrom=rec(candidate))
        source_index[expected] = candidate

    upstream_pruning = []
    for relative in ("docs/validation/storage_20261029/factors/plan.json",
                     "docs/validation/storage_20261029/followup/plan.json"):
        old_plan = read(ROOT / relative)
        upstream_pruning.extend(old_plan["deleteFiles"])
        protect(relative)
    deleted_upstream = {item["path"]: item for item in upstream_pruning}
    for relative, (job, entry, allowed) in TARGETS.items():
        file = OUT / "v83_validation" / relative
        st = file.stat()
        links = aliases[(st.st_dev, st.st_ino)]
        assert len(links) == st.st_nlink and st.st_nlink == 1
        assert {item[0] for item in whosmat(file)} == allowed, "缓存夹带未列明的结果。"
        metadata = loadmat(file, variable_names=list(allowed - {"reduction", "graph", "factor", "extractionReport"}), simplify_cells=True)
        binding = metadata["sourceBinding"]
        protect(binding["sourceFile"], binding["sourceFileSHA256"])
        input_records = []
        for name in ("inputFile", "factorFile", "packedFile", "rhsFile", "assemblyFile"):
            if name not in metadata:
                continue
            path = safe(metadata[name])
            rel = path.relative_to(ROOT).as_posix()
            if path.exists():
                # 同轮计划删除的上游仍先核对身份，避免把它列作保留依赖。
                if rel not in {"outputs/v83_validation/" + p for p in TARGETS}:
                    input_records.append(dict(status="retained", **protect(path)))
                else:
                    input_records.append(dict(status="also-retired-in-this-plan", **rec(path)))
            else:
                assert rel in deleted_upstream, "缺失上游没有历史清理回执：" + rel
                input_records.append(dict(status="previously-pruned-rebuild-required", **deleted_upstream[rel]))
        process_file = OUT / "v83_validation" / (job + "_process.json")
        process = read(process_file)
        expected_exit = 4294967295 if job == "momentumlinereal03" else 0
        assert process["status"] == "terminal" and process["exitCode"] == expected_exit
        assert process["sourceStableAfter"] and all(process["sourceStableAfter"].values())
        protect(process_file)
        if job == "momentumlinereal03":
            # 原进程在巨大报告序列化时中断；已保存数值及后来恢复的摘要必须保留。
            stop = read(file.parent / "agent_reporting_stop_record.json")
            protect(file.parent / "agent_reporting_stop_record.json")
            protect(file.parent / "momentum_line_krylov_evidence.mat", stop["savedNumericalEvidenceSHA256"])
            recovery_dir = OUT / "v83_validation/momentumlinereport02"
            recovery_process = OUT / "v83_validation/momentumlinereport02_process.json"
            recovered = read(recovery_dir / "momentum_lines_real_summary.json")
            recovered_process = read(recovery_process)
            assert recovered_process["status"] == "terminal" and recovered_process["exitCode"] == 0
            assert all(recovered_process["sourceStableAfter"].values())
            assert recovered["sourceBinding"] == binding
            assert recovered["driverReportingInterrupted"] and recovered["originalDriverExitCode"] == expected_exit
            assert recovered["factorAvailable"] and recovered["threeActionsComplete"] and recovered["innerSolveComplete"]
            protect(recovery_process)
            for path in recovery_dir.rglob("*"):
                if path.is_file():
                    protect(path)
            for rel, expected in recovered_process["sourcesBefore"].items():
                if rel.endswith(".m"):
                    source(rel, expected)
        assert any(Path(p).name == entry for p in process["sourcesBefore"])
        for rel, expected in process["sourcesBefore"].items():
            if rel.endswith(".m"):
                source(rel, expected)
        for name, expected in metadata["sourceIdentity"].items():
            matching = [(rel, digest) for rel, digest in process["sourcesBefore"].items()
                        if Path(rel).name == name + ".m" and digest == expected]
            assert matching, "构造来源与实际进程不符：" + name
        artifact = rec(file)
        deleted.append(artifact)
        marker = (file.parent / "DERIVED_CACHE_PRUNED_20261029.json").relative_to(ROOT).as_posix()
        assert not (ROOT / marker).exists()
        groups.append(dict(file=artifact, aliases=links, variables=sorted(allowed),
                           process=rec(process_file), sourceBinding=binding, upstreamInputs=input_records,
                           originalProcessExitCode=expected_exit,
                           rebuildEntry="tests/" + entry, markerPath=marker))
        # 原向量、探针、迭代记录、失败和摘要都保留，不删整个研究目录。
        for path in file.parent.rglob("*"):
            if path.is_file() and path != file:
                protect(path)
        print("CACHE_VERIFIED " + relative, flush=True)

    for relative in ("v93_validation/continuous01/step01/krylov01/wall_velocity_krylov_evidence.mat",
                     "v93_validation/continuous01/step01/direction_context.mat",
                     "v93_validation/continuous01/step01/active_direction.mat",
                     "v80_validation/domain02/shell_domain_evidence.mat"):
        protect(OUT / relative)
    for version in (94, 95, 96):
        for path in (OUT / ("v%d_validation" % version)).rglob("*.mat"):
            protect(path)
    assert not set(protected).intersection(item["path"] for item in deleted)

    # 模块三：透明压缩大文本；不删除JSON、CSV，不改变下载链接或历史内容。
    compression = []
    for key, links in aliases.items():
        relative = links[0]
        item = before["files"][relative]
        if item["bytes"] < 128 * 1024 or Path(relative).suffix.lower() not in TEXT_EXTENSIONS:
            continue
        path = ROOT / relative
        assert len(links) == path.stat().st_nlink, "压缩对象存在工程外硬链接。"
        assert not path.stat().st_file_attributes & 0x800, "本轮只处理尚未压缩的文本。"
        compression.append(dict(rec(path), aliases=links, device=key[0], fileID=key[1],
                                mtimeNS=item["mtimeNS"], allocatedBytesBefore=allocated(path)))
        if len(compression) % 200 == 0:
            print("TEXT_VERIFIED " + str(len(compression)), flush=True)
    DEST.mkdir(parents=True)
    for path, content in archives.items():
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(content)
    inv = DEST / "inventory_before.json.gz"
    with inv.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            stream.write(json.dumps(before).encode("utf-8"))
    data = dict(kind="three-retired-derived-caches-and-transparent-NTFS-text-compression", workspace=str(ROOT),
                createdUTC=datetime.now(timezone.utc).isoformat(),
                sourceCommit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                groups=groups, deleteFiles=deleted, compressFiles=compression,
                protectedFiles=list(protected.values()), sourceFiles=list(versions.values()), inventoryBefore=rec(inv),
                plannedDeletedUniqueBytes=sum(item["bytes"] for item in deleted),
                compressionLogicalBytes=sum(item["bytes"] for item in compression),
                compressionAllocatedBytesBefore=sum(item["allocatedBytesBefore"] for item in compression),
                numericalRebuildExecuted=False,
                note="只删除三个纯构造缓存；旧上游有清理记录。新目录重建全部来源，不承诺历史MAT SHA。文本只用普通NTFS压缩，不使用WOF、ZIP或内容改写。")
    (DEST / "plan.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(deleteFiles=len(deleted), deleteBytes=data["plannedDeletedUniqueBytes"],
                          compressObjects=len(compression), compressionBytes=data["compressionLogicalBytes"],
                          protectedFiles=len(protected), sourceVersions=len(versions))))


def audit():
    """模块四：独立遍历所有输出，核对删除范围、压缩SHA、原身份和实际分配长度。"""
    data = read(DEST / "plan.json")
    receipt = read(DEST / "execution.json")
    assert receipt["status"] == "completed-storage-maintenance"
    assert receipt["planSHA256"] == sha(DEST / "plan.json")
    rec(ROOT / data["inventoryBefore"]["path"], data["inventoryBefore"]["sha256"])
    with gzip.open(DEST / "inventory_before.json.gz", "rb") as stream:
        before = json.load(stream)
    after = inventory()
    deleted = {item["path"] for item in data["deleteFiles"]}
    markers = {group["markerPath"] for group in data["groups"]}
    updates_file = DEST / "readme_update.json"
    updates = read(updates_file) if updates_file.exists() else None
    changed = {"outputs/README.md"} if updates else set()
    assert set(before["files"]) - set(after["files"]) == deleted
    assert set(after["files"]) - set(before["files"]) == markers
    assert all(before["files"][p] == after["files"][p] for p in before["files"] if p not in deleted | changed)
    if updates:
        assert before["files"]["outputs/README.md"] == updates["beforeIdentity"]
        rec(ROOT / updates["savedOriginal"]["path"], updates["beforeFile"]["sha256"])
        rec(ROOT / updates["afterFile"]["path"], updates["afterFile"]["sha256"])
    allocated_after = 0
    for item in data["compressFiles"]:
        rec(ROOT / item["path"], item["sha256"])
        assert (ROOT / item["path"]).stat().st_file_attributes & 0x800, item["path"]
        allocated_after += allocated(ROOT / item["path"])
    for item in data["protectedFiles"] + data["sourceFiles"]:
        rec(ROOT / item["path"], item["sha256"])
    expected = before["uniqueFileLengthBytes"] - data["plannedDeletedUniqueBytes"] + sum(after["files"][p]["bytes"] for p in markers)
    if updates:
        expected += updates["afterFile"]["bytes"] - updates["beforeFile"]["bytes"]
    assert after["uniqueFileLengthBytes"] == expected
    report = dict(status="passed-full-directory-identity-and-content-audit", deletedFiles=len(deleted),
                  unchangedOriginalFiles=before["fileCount"] - len(deleted) - len(changed),
                  separatelyDocumentedREADMEUpdate=bool(updates),
                  compressedUniqueObjects=len(data["compressFiles"]),
                  compressionOriginalBytes=data["compressionLogicalBytes"],
                  compressionAllocatedBytesBefore=data["compressionAllocatedBytesBefore"],
                  compressionAllocatedBytesAfter=allocated_after,
                  compressionSavedAllocatedBytes=data["compressionAllocatedBytesBefore"] - allocated_after,
                  deletedUniqueFileLengthBytes=data["plannedDeletedUniqueBytes"],
                  actualFreeIncreaseBytes=receipt["actualFreeIncreaseBytes"],
                  freeBytesAtInventory=after["freeBytes"],
                  laterFreeIncreaseFromExecutionStartBytes=after["freeBytes"] - receipt["freeBefore"],
                  uniqueFileLengthBytesAfter=after["uniqueFileLengthBytes"],
                  protectedSHARecords=len(data["protectedFiles"]), sourceVersionSHARecords=len(data["sourceFiles"]),
                  currentQ41AndAllUniqueResultsRetained=True, numericalRebuildExecuted=False)
    target = DEST / ("final_postcheck.json" if updates else "postcheck.json")
    assert not target.exists()
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    {"plan": plan, "audit": audit}[sys.argv[1]]()
