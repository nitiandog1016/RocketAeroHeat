"""从冷态计算 Short Test 65 双支管装置，或对已收敛的任务重新核算并导出五区原面。"""
import argparse
import os
from pathlib import Path
import sys

# 周向对照的面积点积与存档后处理均用单线程归约，原生求解的四线程由其启动器设置。
if __name__=="__main__" and "--angular" in sys.argv:
    for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
        os.environ[variable]="1"

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"docs/validation/v329"))
from build_shared import build
from map_surface import mapped
from run_shared import read, sha, solve, write_json


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    p.add_argument("--completed",type=Path,help="从已有的真实终点和两份释放壁行重新映射")
    choice=p.add_mutually_exclusive_group()
    choice.add_argument("--near-wall",action="store_true",help="采用五个间隔各自对半的十层近壁网")
    choice.add_argument("--angular",action="store_true",help="采用十层近壁网，每个周向扇区一分为二")
    a=p.parse_args();job=a.output.resolve();job.relative_to(ROOT/"outputs")
    job.mkdir(parents=True,exist_ok=False)
    run_args={};map_args={}
    if a.angular:
        sys.path.insert(0,str(ROOT/"docs/validation/v332"))
        from sector_grid import inspect
        from shortSupportMesh import prepare
        prep=ROOT/"outputs/v332_validation/prepared01"
        if not (prep/"prepared.json").is_file():
            prep=job/"prepared";prepare(prep,split_near=True,split_sectors=True)
        grid=inspect(prep);transfer=ROOT/"cases/short_sphere/supportedAngularTransfer65.npz"
        grid["publicTransferSHA256"]="63c49c77d06cd2c1814d97f68769c62c3401d52748ad9863f0382fb8408a8609"
        assert sha(transfer)==grid["publicTransferSHA256"]
        write_json(job/"grid.json",grid)
        run_args=dict(prepared=prep,plan=ROOT/"docs/validation/v332/PLAN.md")
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
        version="v332" if a.angular else ("v330" if a.near_wall else "v329")
        driver=ROOT/"outputs"/(version+"_validation")/"build01"
        if not (driver/"build.json").is_file():
            driver=job/"driver";build(driver,run_args.get("prepared"))
        source=job/"run";solve(source,driver,**run_args)
    if a.angular:
        from sector_result import process
        report=process(source,prep,grid,transfer,job)
    else:
        report=mapped(source,job/"surface",**map_args)
    write_json(job/"completed.json",dict(archiveSHA256=sha(job/"surface/surface.zip"),frameSHA256=report["frameSHA256"],
               thermalSHA256=sha(source/"thermal.json"),source=str(source),areaM2=report["areaM2"],faces=report["faces"]))
