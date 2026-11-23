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
    p.add_argument("--near-wall",action="store_true",help="采用五个间隔各自对半的十层近壁网")
    a=p.parse_args();job=a.output.resolve();job.relative_to(ROOT/"outputs")
    job.mkdir(parents=True,exist_ok=False)
    run_args={};map_args={}
    if a.near_wall:
        sys.path.insert(0,str(ROOT/"docs/validation/v330"))
        from check_grid import inspect
        from shortSupportMesh import prepare
        prep=ROOT/"outputs/v330_validation/prepared01"
        if not (prep/"prepared.json").is_file():
            prep=job/"prepared";prepare(prep,split_near=True)
        grid=inspect(prep);write_json(job/"grid.json",grid)
        run_args=dict(prepared=prep,plan=ROOT/"docs/validation/v330/PLAN.md")
        map_args=dict(prepared=prep,expected_mesh=grid["meshSHA256"])
    if a.completed:
        source=a.completed.resolve();source.relative_to(ROOT/"outputs")
    else:
        driver=ROOT/"outputs"/("v330_validation" if a.near_wall else "v329_validation")/"build01"
        if not (driver/"build.json").is_file():
            driver=job/"driver";build(driver,run_args.get("prepared"))
        source=job/"run";solve(source,driver,**run_args)
    report=mapped(source,job/"surface",**map_args)
    write_json(job/"completed.json",dict(archiveSHA256=sha(job/"surface/surface.zip"),frameSHA256=report["frameSHA256"],
               thermalSHA256=sha(source/"thermal.json"),source=str(source),areaM2=report["areaM2"],faces=report["faces"]))
