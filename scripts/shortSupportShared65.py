"""从冷态计算 Short Test 65 双支管装置，或对已收敛的任务重新核算并导出五区原面。"""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"docs/validation/v329"))
from build_shared import build
from map_surface import mapped
from run_shared import read, sha, solve, write_json


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    p.add_argument("--completed",type=Path,help="从已有的真实终点和两份释放壁行重新映射")
    a=p.parse_args();job=a.output.resolve();job.relative_to(ROOT/"outputs")
    job.mkdir(parents=True,exist_ok=False)
    if a.completed:
        source=a.completed.resolve();source.relative_to(ROOT/"outputs")
    else:
        driver=ROOT/"outputs/v329_validation/build01"
        if not (driver/"build.json").is_file():driver=job/"driver";build(driver)
        source=job/"run";solve(source,driver)
    report=mapped(source,job/"surface")
    write_json(job/"completed.json",dict(archiveSHA256=sha(job/"surface/surface.zip"),frameSHA256=report["frameSHA256"],
               thermalSHA256=sha(source/"thermal.json"),source=str(source),areaM2=report["areaM2"],faces=report["faces"]))
