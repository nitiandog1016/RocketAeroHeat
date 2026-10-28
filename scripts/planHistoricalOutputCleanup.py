"""只读核对四组已过时矩阵，保存历史源码和清理计划；不执行删除。"""

from pathlib import Path
import hashlib
import gzip
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DEST = ROOT / "docs/validation/storage_20261028/followup"
SYSTEMS = (
    "v84_validation/q24system",
    "v84_validation/full02/baseline/system2",
    "v84_validation/full02/minimum/system2",
    "v91_validation/continuous01/step01/current",
)
KINDS = {"cache", "operator", "packed", "factor", "scale", "metric"}
hashes = {}


def safe(path, base=ROOT):
    """模块一：逐级排除重解析点，确认最终绝对路径位于指定目录。"""
    path = Path(os.path.abspath(path))
    path.relative_to(base)
    for item in (path, *path.parents):
        if item == base.parent:
            break
        assert not item.is_symlink(), str(item)
        if item.exists():
            assert not (getattr(item.lstat(), "st_file_attributes", 0) & 0x400), str(item)
    return path


def digest(path):
    path = safe(path)
    s = path.stat()
    key = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)
    if key not in hashes:
        h = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(8 * 1024**2), b""):
                h.update(chunk)
        after = path.stat()
        assert (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) == key
        hashes[key] = h.hexdigest()
    return hashes[key]


def record(path, expected=None):
    path = safe(path)
    s = path.stat()
    sha = digest(path)
    assert expected is None or sha == expected, str(path)
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": s.st_size, "sha256": sha}


def inventory():
    """模块二：记录所有原有输出；按文件ID计算唯一数据长度。"""
    files, unique = {}, {}
    for b, ds, fs in os.walk(OUT, followlinks=False):
        for name in ds:
            safe(Path(b) / name, OUT)
        for name in fs:
            p = safe(Path(b) / name, OUT)
            s = p.stat()
            files[p.relative_to(ROOT).as_posix()] = {
                "bytes": s.st_size, "mtimeNS": s.st_mtime_ns,
                "device": s.st_dev, "fileID": s.st_ino, "links": s.st_nlink,
            }
            unique[(s.st_dev, s.st_ino)] = s.st_size
    return {"files": files, "fileCount": len(files),
            "logicalBytes": sum(s["bytes"] for s in files.values()),
            "uniqueFileLengthBytes": sum(unique.values()),
            "freeBytes": shutil.disk_usage(OUT).free}


def main():
    assert not (DEST / "plan.json").exists(), "已有计划，不得覆盖。"
    DEST.mkdir(parents=True, exist_ok=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    modules = {p.stem: p for p in (ROOT / "internalFlow").rglob("*.m")}
    historical = ROOT / "outputs/v84_validation/delivery/isolated/RocketAeroHeat"
    protected, sources, systems, deleted = {}, {}, [], []
    # 模块三：逐项验证清单、输入和全部绑定源码；缺失的历史版本单独存原字节。
    for folder in SYSTEMS:
        directory = safe(OUT / folder, OUT)
        assert not (directory / "CACHE_PRUNED.json").exists()
        manifest = directory / "peak_system.json"
        d = json.loads(manifest.read_text(encoding="utf8"))
        assert d["kind"] == "current-state-peak-direction-system-v1"
        assert set(d["artifacts"]) == KINDS
        protected[manifest] = record(manifest)
        s = d["sourceBinding"]
        source = safe(s["sourceFile"], OUT)
        protected[source] = record(source, s["sourceFileSHA256"])
        g = d["geometryBinding"]
        geometry = safe(g["file"], OUT)
        protected[geometry] = record(geometry, g["sha256"])
        bindings = []
        for category in ("sourceIdentity", "methodIdentity"):
            for name, sha in d[category].items():
                p = modules[name]
                rel = p.relative_to(ROOT).as_posix()
                git_bytes = subprocess.check_output(["git", "show", f"{commit}:{rel}"], cwd=ROOT)
                if hashlib.sha256(git_bytes).hexdigest() == sha:
                    item = record(p, sha)
                    item.update(restorePath=rel, sourceCommit=commit)
                else:
                    old = safe(historical / rel, OUT)
                    assert digest(old) == sha, f"找不到匹配的历史源码：{name}"
                    archive = DEST / "historical_sources" / f"{name}_{sha[:12]}.m.raw"
                    archive.parent.mkdir(exist_ok=True)
                    if not archive.exists():
                        archive.write_bytes(old.read_bytes())
                    item = record(archive, sha)
                    item.update(restorePath=rel, sourceCommit=None,
                                preservedFrom=old.relative_to(ROOT).as_posix())
                sources[(item["path"], sha)] = item
                bindings.append(item)
        artifacts = []
        for kind, a in d["artifacts"].items():
            p = safe(a["file"], OUT)
            assert p.parent == directory and p.name == f"current_peak_{kind}_evidence.mat"
            assert p.stat().st_nlink == 1, "共享文件需要另行核对，不纳入本次删除。"
            item = record(p, a["sha256"])
            artifacts.append(item)
            deleted.append(item)
        systems.append({"manifest": record(manifest), "source": protected[source],
                        "geometry": protected[geometry], "sourceModules": bindings,
                        "artifacts": artifacts,
                        "rebuildFunction": "prepareRocketCurrentPeakSystem(newFolder, sourceFile, geometryFile, 'domain')",
                        "rebuildRequiresNewManifest": True,
                        "historicalMATFileSHAReproductionClaimed": False})
        print(f"VERIFIED {folder}", flush=True)
    # 模块四：保留当前Q39、最近自由/受限方向和全部新版小域试验，不删除独有结果。
    for rel in ("v83_validation/krylov39/wall_velocity_krylov_evidence.mat",
                "v91_validation/continuous01/step01/active_direction.mat"):
        p = OUT / rel
        protected[p] = record(p)
    for base in (OUT / "v91_validation/continuous01/step01/free",
                 OUT / "v92_validation"):
        if base.exists():
            for p in base.rglob("*"):
                if p.is_file():
                    protected[p] = record(p)
    assert len(deleted) == 24 and not {a["path"] for a in deleted}.intersection(a["path"] for a in protected.values())
    before = inventory()
    inventory_file = DEST / "inventory_before.json.gz"
    with inventory_file.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(json.dumps(before).encode("utf8"))
    plan = {"kind": "historical-output-cache-cleanup-v2", "workspace": str(ROOT),
            "status": "sha-verified-plan-not-executed", "createdUTC": datetime.now(timezone.utc).isoformat(),
            "sourceCommit": commit, "cacheSystems": systems,
            "deleteFiles": deleted, "protectedFiles": list(protected.values()),
            "sourceFiles": list(sources.values()), "plannedDeletedBytes": sum(a["bytes"] for a in deleted),
            "inventoryBefore": record(inventory_file),
            "noUniqueCheckpointOrNumericalResultDeleted": True,
            "reason": "Q39 requires fresh derivatives; Q24/Q25/Q38 caches are superseded. Inputs, source, directions, checkpoints and results retained.",
            "rebuildExecuted": False, "historicalBinarySHAReproductionClaimed": False}
    (DEST / "plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="utf8")
    print(json.dumps({"plannedBytes": plan["plannedDeletedBytes"], "files": len(deleted),
                      "protectedFiles": len(protected), "uniqueBytesBefore": before["uniqueFileLengthBytes"]}), flush=True)


if __name__ == "__main__":
    main()
