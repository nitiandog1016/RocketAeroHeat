"""运行 B 型热敏珠固定工况。前驱只提供初值，壁面结果取自强温差三维终场。"""
import argparse
import os
from pathlib import Path

import beadSurfaceResult
import thermistorBead
import thermistorBeadSector
import thermistorReaction
import thermistorTemperature


def run(folder, angular_cells=2):
    if not __debug__:
        raise RuntimeError("运行时需要开启输入和物理检查，请使用默认 Python 模式。")
    folder = folder.resolve()
    folder.relative_to(thermistorBead.ROOT / "outputs")
    folder.mkdir(parents=True, exist_ok=False)
    os.environ["RAH_THERMISTOR_CANCEL_FILE"] = str(folder.parent / "thermistor_cancel.request")
    axis, sector, strong = (folder / name for name in ("axis", "sector", "strong"))
    probes = folder / "probes"
    thermistorBead.prepare(axis, "outer")
    thermistorBead.solve(axis)
    if angular_cells==2:
        thermistorBeadSector.prepare(axis, sector)
    else:
        thermistorBeadSector.prepare(axis, sector, angular_cells)
    thermistorBead.solve(sector)
    thermistorTemperature.prepare(sector, strong)
    thermistorBead.solve(strong)
    thermistorReaction.observe(strong, probes / "reaction0", release=True)
    thermistorReaction.observe(strong, probes / "plus20", release=False)
    thermistorReaction.observe(probes / "plus20", probes / "reaction20", release=True)
    beadSurfaceResult.collect(strong, probes, folder / "result")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_output", type=Path)
    parser.add_argument("--angular-cells",type=int,choices=(2,4),default=2)
    args=parser.parse_args()
    run(args.new_output,args.angular_cells)
