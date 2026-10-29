"""旧 QR 构造缓存清理：先固定输入/来源/硬链接清单，再独立核验删除结果。"""

from datetime import datetime, timezone
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import loadmat, whosmat

from planLegacyOutputCleanup import ROOT, OUT, inventory, rec, safe, sha

DEST = ROOT / "docs/validation/storage_20261029/retired_qr"
CONTEXTS = [
    *(f"v87_validation/continuous01/step{i:02}/direction_context.mat" for i in (1, 2, 3)),
    *(f"v88_validation/continuous01/step{i:02}/direction_context.mat" for i in (1, 2)),
    "v88_validation/recovery01/direction_context.mat",
    *(f"v89_validation/continuous01/step{i:02}/direction_context.mat" for i in (1, 2, 3)),
]
TARGETS = ["v85_validation/qr_real01/augmented_qr_factor.mat",
           "v86_validation/comparison01/temperature_projection_system.mat", *CONTEXTS]
PROCESSES = {"qrreal85": 0, "qrisolated85": 0, "temperaturefull86": 0,
             "temperatureisolated86": 0, "activefull87": 0, "activeisolated87": 0,
             "activefull88": 1, "activerecovery88": 0, "activeisolated88": 0,
             "capacityfull89": 0, "activeisolated89": 0}


def read(path):
    return json.loads(safe(path).read_text(encoding="utf-8-sig"))


def equal(a, b):
    """模块一：逐值核对已保存的温度约束行，不用数值容差替代完全相等。"""
    if sparse.issparse(a) or sparse.issparse(b):
        return sparse.issparse(a) and sparse.issparse(b) and a.shape == b.shape and (a != b).nnz == 0
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, np.ndarray) and a.dtype == object:
        return isinstance(b, np.ndarray) and a.shape == b.shape and all(equal(x, y) for x, y in zip(a.flat, b.flat))
    try:
        return np.array_equal(a, b, equal_nan=True)
    except TypeError:
        return np.array_equal(a, b)


def plan():
    """模块二：只处理十一份构造缓存，保留所有原状态、求解结果与检查点。"""
    assert not DEST.exists(), "已有清理记录，禁止覆盖。"
    before = inventory()
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    aliases = {}
    for rel, item in before["files"].items():
        aliases.setdefault((item["device"], item["fileID"]), []).append(rel)
    protected, versions, archives, groups, deleted = {}, {}, {}, [], []
    historical = {}
    for base in (ROOT / "docs/validation/storage_20261028", ROOT / "docs/validation/storage_20261029"):
        for path in base.rglob("*.raw"):
            historical.setdefault(sha(path), path)

    def protect(path, expected=None):
        path = Path(path)
        path = safe(path if path.is_absolute() else ROOT / path)
        item = rec(path, expected)
        protected[item["path"]] = item

    def source(rel, expected):
        rel = rel.replace("\\", "/")
        if (rel, expected) in versions:
            return
        current = safe(ROOT / rel)
        path = current if current.is_file() and sha(current) == expected else historical.get(expected)
        if path is None:
            blob = subprocess.run(["git", "show", commit + ":" + rel], cwd=ROOT, capture_output=True)
            if blob.returncode == 0 and hashlib.sha256(blob.stdout).hexdigest() == expected:
                content = blob.stdout
                origin = {"sourceCommit": commit}
            else:
                # 早期隔离复核保留过不同版本；只归档匹配SHA的小型源码原文。
                candidates = [ROOT / p for p in before["files"] if Path(p).name == current.name]
                original = next((p for p in candidates if sha(p) == expected), None)
                assert original is not None, "缺少历史源码：" + rel
                content = original.read_bytes()
                origin = {"originalFile": rec(original)}
            path = DEST / "historical_sources" / (current.name + "." + expected[:12] + ".raw")
            archives[path] = content
            versions[(rel, expected)] = dict(path=path.relative_to(ROOT).as_posix(), bytes=len(content),
                                            sha256=expected, restorePath=rel, **origin)
        else:
            versions[(rel, expected)] = dict(rec(path, expected), restorePath=rel)

    for job, expected_exit in PROCESSES.items():
        path = OUT / "v83_validation" / (job + "_process.json")
        process = read(path)
        assert process["status"] == "terminal" and process["exitCode"] == expected_exit
        assert process["sourceStableAfter"] and all(process["sourceStableAfter"].values())
        protect(path)
        for rel, expected in process["sourcesBefore"].items():
            if rel.endswith(".m"):
                source(rel, expected)

    for relative in TARGETS:
        file = safe(OUT / relative)
        st = file.stat()
        links = aliases[(st.st_dev, st.st_ino)]
        assert len(links) == st.st_nlink, "发现本轮目录外的硬链接，不能声称可全部回收。"
        assert not (file.parent / "QR_CACHE_PRUNED.json").exists()
        variables = {x[0] for x in whosmat(file)}
        allowed = {"preparationFile", "preparationSHA256", "factor"}
        if relative in CONTEXTS:
            allowed |= {"rows", "p"}
        elif file.name == "temperature_projection_system.mat":
            allowed |= {"rows", "projector"}
        assert variables == allowed, (relative, variables)
        metadata = loadmat(file, variable_names=["preparationFile", "preparationSHA256", "p", "rows"], simplify_cells=True)
        prep = Path(metadata["preparationFile"])
        prep = safe(prep if prep.is_absolute() else ROOT / prep)
        protect(prep, metadata["preparationSHA256"])
        p = read(prep)
        protect(p["sourceBinding"]["sourceFile"], p["sourceBinding"]["sourceFileSHA256"])
        protect(p["geometryBinding"]["file"], p["geometryBinding"]["sha256"])
        protect(prep.parent / "CACHE_PRUNED.json")
        rows_duplicated = None
        if relative in CONTEXTS:
            active = file.parent / "active_direction.mat"
            saved = loadmat(active, variable_names=["rows", "preparationFile", "preparationSHA256"], simplify_cells=True)
            assert metadata["preparationSHA256"] == saved["preparationSHA256"]
            assert equal(metadata["rows"], saved["rows"]), "约束行不是已保留的完全相同副本。"
            assert metadata["p"]["sourceBinding"] == p["sourceBinding"]
            rows_duplicated = True
        # 同目录的方向、实际作用、检查点、原失败和接受状态全部保留。
        for candidate in file.parent.rglob("*"):
            if candidate.is_file() and candidate != file:
                protect(candidate)
        artifacts = [rec(ROOT / alias, sha(file)) for alias in links]
        deleted.extend(artifacts)
        groups.append(dict(file=rec(file), aliases=artifacts, variables=sorted(variables),
                           preparationFile=rec(prep), rowsExactlyRetainedInActiveDirection=rows_duplicated,
                           markerPaths=sorted({str(Path(x).parent / "QR_CACHE_PRUNED.json").replace("\\", "/") for x in links})))
        print("VERIFIED " + relative, flush=True)

    for rel in ["v93_validation/continuous01/step01/direction_context.mat",
                "v93_validation/continuous01/step01/active_direction.mat",
                "v93_validation/continuous01/step01/krylov01/wall_velocity_krylov_evidence.mat",
                "v80_validation/domain02/shell_domain_evidence.mat"]:
        protect(OUT / rel)
    for version in (94, 95):
        for file in (OUT / f"v{version}_validation").rglob("*"):
            if file.is_file():
                protect(file)
    assert len(groups) == 11 and not set(protected).intersection(x["path"] for x in deleted)
    DEST.mkdir(parents=True)
    for path, content in archives.items():
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(content)
    inv = DEST / "inventory_before.json.gz"
    with inv.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            stream.write(json.dumps(before).encode("utf8"))
    plan_data = dict(kind="eleven-retired-QR-construction-cache-objects", workspace=str(ROOT),
                     createdUTC=datetime.now(timezone.utc).isoformat(), sourceCommit=commit,
                     groups=groups, deleteFiles=deleted, protectedFiles=list(protected.values()),
                     sourceFiles=list(versions.values()), inventoryBefore=rec(inv),
                     plannedUniqueBytes=sum(g["file"]["bytes"] for g in groups),
                     uniqueBytesBefore=before["uniqueFileLengthBytes"], numericalRebuildExecuted=False,
                     note="只回收QR分解/投影构造缓存；全部原结果、失败、检查点、Q41及v94/v95保留。旧审计需在新目录重建，不保证历史MAT SHA。")
    (DEST / "plan.json").write_text(json.dumps(plan_data, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    print(json.dumps(dict(objects=len(groups), paths=len(deleted), uniqueBytes=plan_data["plannedUniqueBytes"],
                          protected=len(protected), sourceVersions=len(versions))))


def audit():
    """模块三：独立遍历全目录并复核保留文件，防止将硬链接累计量算成实收。"""
    plan_data = read(DEST / "plan.json")
    receipt = read(DEST / "execution.json")
    assert receipt["status"] == "completed-eleven-retired-QR-cache-objects"
    assert sha(DEST / "plan.json") == receipt["planSHA256"]
    rec(ROOT / plan_data["inventoryBefore"]["path"], plan_data["inventoryBefore"]["sha256"])
    with gzip.open(DEST / "inventory_before.json.gz", "rb") as stream:
        before = json.load(stream)
    after = inventory()
    deleted = {x["path"] for x in plan_data["deleteFiles"]}
    assert set(before["files"]) - set(after["files"]) == deleted
    assert all(before["files"][p] == after["files"][p] for p in before["files"] if p not in deleted)
    markers = {p for g in plan_data["groups"] for p in g["markerPaths"]}
    assert set(after["files"]) - set(before["files"]) == markers
    for item in plan_data["protectedFiles"] + plan_data["sourceFiles"]:
        rec(ROOT / item["path"], item["sha256"])
    expected = before["uniqueFileLengthBytes"] - plan_data["plannedUniqueBytes"] + sum(after["files"][p]["bytes"] for p in markers)
    assert after["uniqueFileLengthBytes"] == expected
    report = dict(status="passed-full-inventory-and-preserved-file-SHA-audit", deletedObjects=11,
                  deletedPaths=len(deleted), unchangedOriginalFiles=before["fileCount"]-len(deleted),
                  protectedFiles=len(plan_data["protectedFiles"]), sourceVersions=len(plan_data["sourceFiles"]),
                  uniqueBytesBefore=before["uniqueFileLengthBytes"], uniqueBytesAfter=after["uniqueFileLengthBytes"],
                  actualFreeIncreaseBytes=receipt["actualFreeIncreaseBytes"], currentQ41AndV94V95Retained=True,
                  originalFailedV88StudyRetained=True, numericalComputationExecuted=False)
    assert not (DEST / "postcheck.json").exists()
    (DEST / "postcheck.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf8")
    print(json.dumps(report))


if __name__ == "__main__":
    {"plan": plan, "audit": audit}[sys.argv[1]]()
