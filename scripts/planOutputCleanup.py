"""生成一次性的 outputs 清理清单；只读计算数据，不执行删除或链接替换。"""

from collections import defaultdict
from pathlib import Path
import hashlib
import json
import os
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DEST = ROOT / "docs/validation/storage_20261028/outputs_cleanup_plan.json"
hashes = {}


def relative(path):
    path = path.resolve(strict=True)
    path.relative_to(OUT)
    return path.relative_to(ROOT).as_posix()


def digest(path):
    s = path.stat()
    key = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)
    if key not in hashes:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
                h.update(chunk)
        after = path.stat()
        assert (after.st_ino, after.st_size, after.st_mtime_ns) == key[1:]
        hashes[key] = h.hexdigest()
    return hashes[key]


def record(path, expected=None):
    s = path.stat()
    sha = digest(path)
    if expected is not None:
        assert sha == expected, str(path)
    return {"path": relative(path), "bytes": s.st_size, "sha256": sha}


def main():
    started = time.time()
    assert not DEST.exists(), "已有清理计划，不能覆盖或重复执行。"
    systems, deleted, protected, skipped = [], [], {}, []
    modules = {p.stem: p for p in (ROOT / "internalFlow").rglob("*.m")}
    # 模块一：仅纳入已结束的旧工况缓存，原输入与全部构造模块必须仍然匹配。
    for manifest in sorted(OUT.rglob("peak_system.json")):
        rel = relative(manifest)
        if "/delivery/" in rel or "/isolated/" in rel:
            continue
        if not any(rel.startswith(f"outputs/v{v}_validation/") for v in range(84, 90)):
            continue
        d = json.loads(manifest.read_text(encoding="utf8"))
        source_modules = []
        mismatch = []
        for category in ("sourceIdentity", "methodIdentity"):
            for name, expected in d[category].items():
                p = modules.get(name)
                if p is None or digest(p) != expected:
                    mismatch.append(name)
                else:
                    source_modules.append({"path": p.relative_to(ROOT).as_posix(), "sha256": expected})
        if mismatch:
            skipped.append({"manifest": rel, "reason": "builder-source-changed", "modules": mismatch})
            continue
        source = record(Path(d["sourceBinding"]["sourceFile"]), d["sourceBinding"]["sourceFileSHA256"])
        geometry = record(Path(d["geometryBinding"]["file"]), d["geometryBinding"]["sha256"])
        protected[source["path"]] = source
        protected[geometry["path"]] = geometry
        protected[rel] = record(manifest)
        artifacts = []
        assert set(d["artifacts"]) == {"cache", "operator", "packed", "factor", "scale", "metric"}
        for kind, item in d["artifacts"].items():
            p = Path(item["file"])
            assert p.parent == manifest.parent and p.name == f"current_peak_{kind}_evidence.mat"
            a = record(p, item["sha256"])
            assert p.stat().st_nlink == 1, "缓存存在其他链接，需单独审计。"
            a["reason"] = "rebuildable-obsolete-state-matrix-cache"
            artifacts.append(a)
            deleted.append(a)
        systems.append({"manifest": rel, "source": source, "geometry": geometry,
                        "sourceModules": source_modules, "artifacts": artifacts,
                        "rebuildFunction": "prepareRocketCurrentPeakSystem(newFolder, sourceFile, geometryFile, 'domain')",
                        "rebuildRequiresNewManifest": True,
                        "historicalMATFileSHAReproductionClaimed": False})
        print(f"CACHE_PLAN {rel}", flush=True)
    assert len(systems) == 10, "本次核定范围应为 v86—v89 的十组旧系统。"
    # 模块二：误将数值基向量展开成 1 GB JSON 的失败报告，可由保留 MAT 追溯。
    failed = OUT / "v83_validation/momentumlinereal03"
    original = json.loads((failed / "reporting_artifact_receipt.json").read_text(encoding="utf8"))
    stopped = json.loads((failed / "agent_reporting_stop_record.json").read_text(encoding="utf8"))
    numerical = record(failed / "momentum_line_krylov_evidence.mat", stopped["savedNumericalEvidenceSHA256"])
    protected[numerical["path"]] = numerical
    bad_report = record(Path(original["archivedPath"]), original["sha256"])
    assert bad_report["bytes"] == original["bytes"]
    bad_report["reason"] = "failed-expanded-JSON-export-numerical-MAT-and-receipts-retained"
    deleted.append(bad_report)
    # 模块三：当前状态、几何、方向及其完整准备产物全部保留并核对。
    for p in [OUT / "v83_validation/krylov38/wall_velocity_krylov_evidence.mat",
              OUT / "v90_validation/recovery01/active_direction.mat"]:
        protected[relative(p)] = record(p)
    for p in (OUT / "v90_validation/continuous01/step01/current").iterdir():
        if p.is_file():
            protected[relative(p)] = record(p)
    delete_paths = {a["path"] for a in deleted}
    assert not delete_paths.intersection(protected)
    # 模块四：按文件 ID 扣除现有硬链接；仅把逐字节相同的旧数据列为合并对象。
    physical = {}
    logical_bytes = 0
    for base, dirs, files in os.walk(OUT, followlinks=False):
        dirs[:] = [d for d in dirs if not (Path(base) / d).is_symlink()]
        for name in files:
            p = Path(base) / name
            assert not p.is_symlink(), str(p)
            s = p.stat()
            rel = relative(p)
            logical_bytes += s.st_size
            item = physical.setdefault((s.st_dev, s.st_ino), {"path": p, "bytes": s.st_size, "links": s.st_nlink, "paths": []})
            item["paths"].append(rel)
    by_size = defaultdict(list)
    for item in physical.values():
        if item["bytes"] >= 1024**2 and not delete_paths.intersection(item["paths"]):
            by_size[item["bytes"]].append(item)
    merges = []
    for size, candidates in sorted(by_size.items(), reverse=True):
        if len(candidates) < 2:
            continue
        # 头尾筛查只减少读取；最后仍以全文件 SHA256 确认，不能以大小代替内容核验。
        sampled = defaultdict(list)
        for item in candidates:
            p = item["path"]
            with p.open("rb") as stream:
                head = stream.read(65536)
                stream.seek(max(0, size - 65536))
                tail = stream.read(65536)
            sampled[hashlib.sha256(head + tail).digest()].append(item)
        for group in sampled.values():
            if len(group) < 2:
                continue
            exact = defaultdict(list)
            for item in group:
                exact[digest(item["path"])].append(item)
            for sha, items in exact.items():
                if len(items) < 2:
                    continue
                def rank(item):
                    paths = item["paths"]
                    is_current = any(x.startswith(("outputs/v90_validation/", "outputs/v91_validation/")) or x in protected for x in paths)
                    return (not is_current, min("/delivery/" in x or "/isolated/" in x for x in paths), min(paths))
                items.sort(key=rank)
                canonical = min(items[0]["paths"], key=lambda x: ("/delivery/" in x, x))
                for item in items[1:]:
                    paths = sorted(item["paths"])
                    if item["links"] != len(paths):
                        continue  # 存在 outputs 之外的硬链接，保持不动。
                    if any(x in protected or x.startswith(("outputs/v90_validation/", "outputs/v91_validation/")) for x in paths):
                        continue
                    merges.append({"canonical": canonical, "paths": paths, "bytes": size, "sha256": sha})
    plan = {"kind": "explicit-output-cleanup-v1", "workspace": str(ROOT), "scope": "outputs",
            "status": "sha-verified-plan-not-yet-executed", "cacheSystems": systems,
            "deleteFiles": deleted, "deduplicatePhysicalFiles": merges, "protectedFiles": list(protected.values()),
            "skippedSystems": skipped, "before": {"logicalBytes": logical_bytes, "uniqueBytes": sum(a["bytes"] for a in physical.values())},
            "plannedDeletedBytes": sum(a["bytes"] for a in deleted),
            "plannedDeduplicatedBytes": sum(a["bytes"] for a in merges),
            "noUniqueCheckpointOrNumericalEvidenceDeleted": True, "historicalCachesRequireExplicitRebuild": True,
            "hashesComputed": len(hashes), "planSeconds": time.time() - started}
    DEST.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    print(json.dumps({k: plan[k] for k in ("plannedDeletedBytes", "plannedDeduplicatedBytes", "before", "hashesComputed", "planSeconds")}), flush=True)


if __name__ == "__main__":
    main()
