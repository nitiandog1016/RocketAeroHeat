"""检查磁盘缓存损坏时显式拒收；只写验证目录，不修改已生成的体网。"""
import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from prepareRFZAirVolume import prepare


def verify(ticket, run):
    good = Path(json.loads(ticket.read_bytes())["directory"])
    assert prepare(good.parent) == good
    for problem in ("producer", "asset"):
        folder = run / problem / good.name
        folder.mkdir(parents=True, exist_ok=False)
        info = json.loads((good / "air.json").read_bytes())
        if problem == "producer":
            info["source"]["gmshVersion"] = "unverified"
            os.link(good / "air.mat", folder / "air.mat")
        else:
            (folder / "air.mat").write_bytes(b"incomplete MAT")
        (folder / "air.json").write_text(json.dumps(info), encoding="utf8")
        try:
            prepare(folder.parent)
        except AssertionError:
            print(problem, "rejected")
        else:
            raise ValueError("错误缓存被接收")
    assert prepare(good.parent) == good


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ticket", type=Path)
    parser.add_argument("new_directory", type=Path)
    args = parser.parse_args()
    verify(args.ticket, args.new_directory)
