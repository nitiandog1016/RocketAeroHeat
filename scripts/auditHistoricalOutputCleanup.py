"""独立复核旧矩阵清理：仅24文件缺失，其他原文件元数据及保留内容不变。"""

from pathlib import Path
from collections import Counter
import gzip
import hashlib
import json
import os
import shutil

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "docs/validation/storage_20261028/followup"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(8 * 1024**2), b""):
            h.update(data)
    return h.hexdigest()


def main():
    # 模块一：绑定本次计划、执行回执和完整清理前文件目录。
    plan = json.loads((DEST / "plan.json").read_text(encoding="utf8"))
    receipt = json.loads((DEST / "execution.json").read_text(encoding="utf-8-sig"))
    assert receipt["status"] == "completed-sha-verified-historical-cache-cleanup"
    assert sha(DEST / "plan.json") == receipt["planSHA256"]
    inventory = ROOT / plan["inventoryBefore"]["path"]
    assert sha(inventory) == plan["inventoryBefore"]["sha256"]
    before = json.loads(gzip.decompress(inventory.read_bytes()).decode("utf8"))
    deleted = {a["path"] for a in plan["deleteFiles"]}
    assert deleted == {a["path"] for a in receipt["deleted"]} and len(deleted) == 24
    # 模块二：独立遍历；除四份清理标记外没有新增文件，非目标文件均须逐项不变。
    files, unique, groups = {}, {}, Counter()
    for b, ds, fs in os.walk(ROOT / "outputs", followlinks=False):
        for name in ds + fs:
            assert not (getattr((Path(b) / name).lstat(), "st_file_attributes", 0) & 0x400)
        for name in fs:
            p = Path(b) / name
            rel = p.relative_to(ROOT).as_posix()
            s = p.stat()
            files[rel] = {"bytes": s.st_size, "mtimeNS": s.st_mtime_ns,
                          "device": s.st_dev, "fileID": s.st_ino, "links": s.st_nlink}
            key = (s.st_dev, s.st_ino)
            if key not in unique:
                unique[key] = s.st_size
                groups[rel.split("/")[1]] += s.st_size
    old = before["files"]
    documentation = []
    readme_receipt = DEST / "readme_update.json"
    if readme_receipt.exists():
        a = json.loads(readme_receipt.read_text(encoding="utf8"))
        assert a["path"] == "outputs/README.md" and old[a["path"]] == a["beforeMetadata"]
        assert sha(ROOT / a["path"]) == a["afterSHA256"] and files[a["path"]] == a["afterMetadata"]
        documentation.append(a["path"])
    markers = {a["manifest"]["path"].replace("peak_system.json", "CACHE_PRUNED.json")
               for a in plan["cacheSystems"]}
    assert set(old) - set(files) == deleted
    assert set(files) - set(old) == markers
    assert all(files[k] == old[k] for k in set(old) - deleted - set(documentation))
    for a in plan["deleteFiles"]:
        assert old[a["path"]]["links"] == 1 and old[a["path"]]["bytes"] == a["bytes"]
    # 模块三：保留检查点、几何、方向、小域试验及所有构造源码再次完整验SHA。
    for a in plan["protectedFiles"] + plan["sourceFiles"]:
        p = ROOT / a["path"]
        assert p.stat().st_size == a["bytes"] and sha(p) == a["sha256"], a["path"]
    deleted_bytes = sum(a["bytes"] for a in plan["deleteFiles"])
    doc_delta = sum(files[k]["bytes"] - old[k]["bytes"] for k in documentation)
    assert sum(unique.values()) == before["uniqueFileLengthBytes"] - deleted_bytes + sum(files[k]["bytes"] for k in markers) + doc_delta
    report = {"status": "passed-independent-inventory-and-SHA-audit",
              "missingFilesExactlyMatchPlan": True, "deletedFileCount": 24,
              "allOtherOriginalDataFileMetadataUnchanged": True,
              "intentionalDocumentationUpdates": documentation,
              "protectedAndSourceFilesSHAVerified": True,
              "deletedBytes": deleted_bytes,
              "actualFreeIncreaseBytesFromExecution": receipt["actualFreeIncreaseBytes"],
              "uniqueFileLengthBytesBefore": before["uniqueFileLengthBytes"],
              "uniqueFileLengthBytesAfter": sum(unique.values()),
              "logicalPathBytesAfter": sum(a["bytes"] for a in files.values()),
              "currentFreeBytes": shutil.disk_usage(ROOT / "outputs").free,
              "uniqueBytesGroupedByFirstPath": dict(groups.most_common()),
              "currentState": "Q39; unchanged; not steady",
              "noRebuildOrPhysicalComputationExecuted": True}
    (DEST / "postcheck.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    print(json.dumps({k: v for k, v in report.items() if k != "uniqueBytesGroupedByFirstPath"}), flush=True)


if __name__ == "__main__":
    main()
