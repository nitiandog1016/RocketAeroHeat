"""在 B 型珠体的原扇区内求 6 m/s 启动基流，供完整周向网格初始化。"""
import argparse
import os
from pathlib import Path

from runSU2NativeDebugger import observe
from thermistorBead import ENGINE, ENGINE_SHA, ROOT, sha
from thermistorFullRing import PARENT, SOURCE

FULL = ROOT / "outputs/v287_validation/first02"
FULL_CONFIG_SHA = "5fa790e22f7ebb15f2a0bfb650299bdefa0ad2e09d12f18613450ba1af050092"
CONTRACT = ROOT / "docs/validation/v288/PLAN.md"


def matched_config():
    """从固定全周配置读取 6 m/s 的三个输运系数，镜面、方程和原生容差沿用扇区。"""
    assert sha(FULL/"case.cfg") == FULL_CONFIG_SHA and sha(PARENT/"case.cfg") == SOURCE["case.cfg"]
    fields = dict(line.strip().split("= ",1) for line in (FULL/"case.cfg").read_text("ascii").splitlines()
                  if "= " in line)
    replace = {key: fields[key] for key in ("MU_CONSTANT","THERMAL_CONDUCTIVITY_CONSTANT","DIFFUSIVITY_CONSTANT")}
    replace.update(ITER="7000", OUTPUT_FILES="(RESTART_ASCII)", OUTPUT_WRT_FREQ="800")
    lines = []
    for line in (PARENT/"case.cfg").read_text("ascii").splitlines():
        key = line.partition("=")[0].strip()
        if key in replace:
            line = f"{key}= {replace.pop(key)}"
        lines.append(line)
    assert not replace
    return "\n".join(lines)+"\n"


def run(folder):
    """创建独立的稳态求解目录，并锁定实际参与本次启动的几何、初场、程序和配置。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/"outputs")
    assert sha(ENGINE) == ENGINE_SHA
    for name, digest in SOURCE.items():
        assert sha(PARENT/name) == digest, name
    cfg = matched_config()
    folder.mkdir(parents=True, exist_ok=False)
    inputs = [Path(__file__).resolve(), CONTRACT, ENGINE, PARENT/"case.cfg", FULL/"case.cfg",
              ROOT/"scripts/runSU2NativeDebugger.py", ROOT/"scripts/thermistorFullRing.py",
              ROOT/"scripts/thermistorBead.py"]
    for source, dest in (("restart.csv","initial_flow.csv"),("mesh.su2","mesh.su2"),
                         ("nodes.npz","nodes.npz"),("body.json","body.json")):
        os.link(PARENT/source, folder/dest)
        inputs.append(folder/dest)
    (folder/"case.cfg").write_text(cfg,encoding="ascii")
    inputs.append(folder/"case.cfg")
    receipt = observe([str(ENGINE),"case.cfg","--threads","1"], folder,
                      {str(p):sha(p) for p in inputs},600,180*1024**2)
    return receipt["debuggerExitCode"] == 0 and receipt["stopReason"] is None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    raise SystemExit(0 if run(parser.parse_args().output) else 1)
