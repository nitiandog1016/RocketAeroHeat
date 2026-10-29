"""核对两组过时的 v83 派生缓存；plan 只读规划，audit 独立复核实际删除。"""

from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DEST = ROOT / "docs/validation/storage_20261029"
GROUPS = {
    "currentparent01": ("reservoircurrentparent01", [
        "current_operator_evidence.mat", "current_packed_evidence.mat",
        "current_factor_evidence.mat", "current_image_scale_evidence.mat"]),
    "peakweightreal01": ("peakweightreal01", [
        "current_peak_cache_evidence.mat", "current_peak_operator_evidence.mat",
        "current_peak_packed_evidence.mat", "current_peak_factor_evidence.mat",
        "current_peak_scale_evidence.mat"]),
}
EXCLUDED = {"exportRocketWallVelocityStudy.m", "tests/continueRocketWallVelocityKrylov.m",
            "internalFlow/cutCell/solveRocketRegularizedParentDirection.m"}
HASHES = {}


def safe(path):
    """模块一：逐级拒绝重解析点，最终路径必须位于当前工程。"""
    path = Path(os.path.abspath(path))
    path.relative_to(ROOT)
    for p in (path, *path.parents):
        if p == ROOT.parent:
            break
        if p.exists():
            assert not (getattr(p.lstat(), "st_file_attributes", 0) & 0x400), str(p)
    return path


def sha(path):
    path = safe(path)
    st = path.stat()
    key = (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns)
    if key not in HASHES:
        h = hashlib.sha256()
        with path.open("rb") as f:
            for block in iter(lambda: f.read(8 * 1024**2), b""):
                h.update(block)
        now = path.stat()
        assert key == (now.st_dev, now.st_ino, now.st_size, now.st_mtime_ns)
        HASHES[key] = h.hexdigest()
    return HASHES[key]


def rec(path, expected=None):
    path = safe(path)
    value = sha(path)
    assert expected is None or value == expected, str(path)
    return dict(path=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size, sha256=value)


def inventory():
    files, unique = {}, {}
    for b, ds, fs in os.walk(OUT, followlinks=False):
        for n in ds:
            safe(Path(b) / n)
        for n in fs:
            p = safe(Path(b) / n)
            s = p.stat()
            files[p.relative_to(ROOT).as_posix()] = dict(
                bytes=s.st_size, mtimeNS=s.st_mtime_ns, fileID=s.st_ino, device=s.st_dev)
            unique[(s.st_dev, s.st_ino)] = s.st_size
    return dict(files=files, fileCount=len(files), pathBytes=sum(x["bytes"] for x in files.values()),
                uniqueFileLengthBytes=sum(unique.values()), freeBytes=shutil.disk_usage(OUT).free)


def plan():
    """模块二：只选择九个构造缓存，保留原候选、实际模型、结果和当前依赖。"""
    assert not (DEST / "plan.json").exists(), "本次已有计划，不能覆盖。"
    DEST.mkdir(exist_ok=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    groups, deleted, protected, sources = [], [], {}, {}
    base = OUT / "v83_validation"
    for name, (job, names) in GROUPS.items():
        folder = base / name
        assert not (folder / "CACHE_PRUNED.json").exists()
        process = base / (job + "_process.json")
        p = json.loads(process.read_text(encoding="utf8"))
        assert p["status"] == "terminal" and p["exitCode"] == 0
        assert all(p["sourceStableAfter"].values())
        protected[process] = rec(process)
        bindings = []
        # 保存完整历史 MATLAB 源快照中与构造有关的版本；三个非构造模块不纳入。
        for rel, expected in p["sourcesBefore"].items():
            rel = rel.replace("\\", "/")
            if not rel.endswith(".m") or rel in EXCLUDED:
                continue
            source = ROOT / rel
            if sha(source) != expected:
                assert rel == "internalFlow/cutCell/cutCellExactEulerBranchJacobian.m", rel
                source = ROOT / "docs/validation/storage_20261028/followup/historical_sources" / (
                    "cutCellExactEulerBranchJacobian_" + expected[:12] + ".m.raw")
            item = rec(source, expected)
            item["restorePath"] = rel
            sources[(rel, expected)] = item
            bindings.append(dict(restorePath=rel, sha256=expected, preservedPath=item["path"]))
        artifacts = []
        for filename in names:
            path = safe(folder / filename)
            path.relative_to(OUT)
            assert path.stat().st_nlink == 1
            item = rec(path)
            artifacts.append(item)
            deleted.append(item)
        for path in folder.iterdir():
            if path.is_file() and path.name not in names:
                protected[path] = rec(path)
        groups.append(dict(folder=folder.relative_to(ROOT).as_posix(), process=rec(process),
                           sourceModules=bindings, removedArtifacts=artifacts))
    # 模块三：逐字节核对 Git 已保存的构造源码，未包含于 Git 的原字节另存小文件。
    needed = sorted({x["path"] for x in sources.values() if not x["path"].endswith(".raw")})
    batch = subprocess.check_output(["git", "cat-file", "--batch"], cwd=ROOT,
                                   input="".join(commit + ":" + p + "\n" for p in needed).encode())
    blobs = {}; pos = 0
    for rel in needed:
        end = batch.index(b"\n", pos); header = batch[pos:end].split(); pos = end + 1
        if header[-1:] == [b"missing"]:
            continue  # 历史运行锁也包含 outputs 中的小型原始源码，后续单独保存原字节。
        assert len(header) == 3 and header[1] == b"blob", rel
        length = int(header[2]); data = batch[pos:pos+length]; pos += length + 1
        blobs[rel] = (header[0].decode(), hashlib.sha256(data).hexdigest())
    for item in sources.values():
        if item["path"] in blobs and blobs[item["path"]][1] == item["sha256"]:
            item.update(sourceCommit=commit, gitBlob=blobs[item["path"]][0])
        else:
            source = ROOT / item["path"]
            target = DEST / "historical_sources" / (Path(item["restorePath"]).name + "." + item["sha256"][:12] + ".raw")
            target.parent.mkdir(exist_ok=True)
            if not target.exists():
                target.write_bytes(source.read_bytes())
            item.update(archivedFile=rec(target), sourceCommit=None)
    for rel in ["v80_validation/domain02/shell_domain_evidence.mat",
                "v83_validation/krylov20/wall_velocity_krylov_evidence.mat",
                "v83_validation/krylov22/wall_velocity_krylov_evidence.mat",
                "v83_validation/krylov39/wall_velocity_krylov_evidence.mat",
                "v83_validation/krylov40/wall_velocity_krylov_evidence.mat",
                "v83_validation/reservoirjacreal01/riemann_selected_direction_real_evidence.mat",
                "v83_validation/updateradiusreal01/regularized_parent_candidate_evidence.mat",
                "v92_validation/continuous01/step01/active_direction.mat"]:
        path = OUT / rel
        protected[path] = rec(path)
    assert len(deleted) == 9 and not {x["path"] for x in deleted}.intersection(x["path"] for x in protected.values())
    before = inventory()
    inv = DEST / "inventory_before.json.gz"
    with inv.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            stream.write(json.dumps(before).encode())
    result = dict(kind="v83-two-obsolete-construction-cache-groups", workspace=str(ROOT),
                  createdUTC=datetime.now(timezone.utc).isoformat(), sourceCommit=commit,
                  groups=groups, deleteFiles=deleted, protectedFiles=list(protected.values()),
                  sourceFiles=list(sources.values()), excludedNonConstructionModules=sorted(EXCLUDED),
                  inventoryBefore=rec(inv), plannedBytes=sum(x["bytes"] for x in deleted),
                  numericalRebuildExecuted=False, historicalMATFileSHAReproductionClaimed=False,
                  note="Rebuild constructor stages only, in a new study folder with recorded historical sources; retained candidate/model values are historical evidence, not inputs to a fresh Q40 system.")
    (DEST / "plan.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    print(json.dumps(dict(files=9, plannedBytes=result["plannedBytes"], protected=len(protected),
                          sourceVersions=len(sources), uniqueBytesBefore=before["uniqueFileLengthBytes"])))


def audit():
    """模块四：独立比较全目录，只允许清单九项消失，其他原文件身份不变。"""
    plan = json.loads((DEST / "plan.json").read_text(encoding="utf8"))
    receipt = json.loads((DEST / "execution.json").read_text(encoding="utf-8-sig"))
    assert receipt["status"] == "completed-nine-obsolete-caches-removed"
    with gzip.open(DEST / "inventory_before.json.gz", "rb") as f:
        before = json.load(f)
    after = inventory()
    removed = set(before["files"]) - set(after["files"])
    assert removed == {x["path"] for x in plan["deleteFiles"]}
    changed = [p for p in before["files"] if p in after["files"] and before["files"][p] != after["files"][p]]
    assert not changed, changed
    for item in plan["protectedFiles"] + plan["sourceFiles"]:
        rec(ROOT / item["path"], item["sha256"])
        if "archivedFile" in item:
            rec(ROOT / item["archivedFile"]["path"], item["sha256"])
    result = dict(status="passed-independent-inventory-and-SHA-audit", removed=sorted(removed),
                  unchangedOriginalFiles=before["fileCount"]-len(removed),
                  uniqueBytesBefore=before["uniqueFileLengthBytes"], uniqueBytesAfter=after["uniqueFileLengthBytes"],
                  actualFreeIncreaseBytes=receipt["actualFreeIncreaseBytes"], protectedAndSourceSHA=True,
                  newOutputFiles=sorted(set(after["files"])-set(before["files"])),
                  noNumericalComputationExecuted=True)
    (DEST / "postcheck.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps(result))


if __name__ == "__main__":
    {"plan": plan, "audit": audit}[sys.argv[1]]()
