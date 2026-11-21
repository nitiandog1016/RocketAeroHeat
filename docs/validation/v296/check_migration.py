"""检查热敏珠旧安装迁移、跨卷复制中断与再次使用缓存。全过程只准备输入。"""
import argparse
import errno
import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import thermistorStartupCache as cache
from runThermistorStartupJob import preflight


def check_transfer(directory):
    source = directory / "source"
    target = directory / "target"
    source.write_bytes(b"complete checkpoint")
    digest = cache.sha(source)

    def interrupted_copy(src, pending):
        Path(pending).write_bytes(b"partial")
        raise OSError(errno.ENOSPC, "test: destination is full")

    # 模拟两个磁盘之间无法创建硬链接，随后在复制期间耗尽空间。
    with patch.object(os, "link", side_effect=OSError(errno.EXDEV, "test: different disks")):
        with patch.object(cache.shutil, "copyfile", side_effect=interrupted_copy):
            try:
                cache.link_or_copy(source, target)
            except OSError as error:
                assert error.errno == errno.ENOSPC
            else:
                raise AssertionError("未能复现复制中断")
        assert not target.exists() and cache.sha(source) == digest
        cache.link_or_copy(source, target)
        assert cache.sha(target) == digest and not os.path.samefile(target, source)
    assert not cache.available(directory, [(source.name, "0" * 64)])
    assert cache.sha(source) == digest


def check_existing():
    coarse = ROOT / "outputs/v289_validation/window01"
    fine = ROOT / "outputs/v293_validation/grid01"
    # 只在三套历史输入实际存在且内容匹配时执行这项迁移试验，避免误启冷态计算。
    assert cache.available(coarse, [("mesh.su2", cache.COARSE_MESH), ("nodes.npz", cache.COARSE_NODES),
                                    ("body.json", cache.BODY_SHA)])
    for name, key in (("mesh.su2", "gridSHA256"), ("nodes.npz", "nodeMapSHA256"), ("seed_00000.csv", "seedSHA256")):
        assert cache.available(fine, [(name, cache.GRID_SPEC[key])])
    assert cache.available(fine, [("grid.json", "f796d2d3e3aee6682e2196ac625ab439e1098516bc0dac0832260952faecbf0b")])
    assert cache.available(fine.parent / "build01", [("observe_bead.exe", cache.OBSERVER_SPEC["observerSHA256"]),
                                                    ("CIntegration.cpp", cache.OBSERVER_SPEC["sourceSHA256"])])
    preflight()
    for installed, previous in ((cache.COARSE, coarse), (cache.GRID, fine)):
        assert Path(cache.read(installed / "installed.json")["sourceDirectory"]) == previous
    seed_name = "seed_00000.csv"
    assert cache.sha(cache.GRID / seed_name) == cache.sha(fine / seed_name) == cache.GRID_SPEC["seedSHA256"]
    observer = cache.read(cache.BUILD / "build.json")
    assert observer["origin"] == "verified-local-migration"
    before = sorted(p.name for p in cache.CACHE.iterdir())
    preflight()
    assert before == sorted(p.name for p in cache.CACHE.iterdir())
    return dict(coarseSource=str(coarse), fineSource=str(fine), seedSHA256=cache.GRID_SPEC["seedSHA256"],
                observerSource=observer["origin"], observerSHA256=observer["observerSHA256"],
                secondPreparationReusedCache=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("receipt", type=Path)
    args = cli.parse_args()
    assert not args.receipt.exists()
    with tempfile.TemporaryDirectory(dir=ROOT / "outputs/v296_validation") as work:
        check_transfer(Path(work))
    proof = check_existing()
    with args.receipt.open("x", encoding="utf8") as out:
        json.dump(dict(proof, partialTransferHidden=True, simulatedCrossVolumeRetryPassed=True,
                       mismatchedOptionalInputPreserved=True), out, indent=2)
