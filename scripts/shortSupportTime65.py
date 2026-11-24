"""从给定双支管稳态基流生成扰动初场，求解并导出首个同刻的两套全表面结果。"""
import argparse
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):os.environ[name]="1"
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"docs/validation/v333"))
from build_time import BASE_FRAME, build, sha, read, write_json
from physical_step import setup, run_steps
from time_surface import export


def main(out,completed=None):
    out.relative_to(ROOT/"outputs");out.mkdir(parents=True,exist_ok=False)
    # 旧安装保留已验收的冻结入口；没有历史基流的新安装从本项目的几何和冷初场计算。
    legacy=ROOT/"outputs/v332_validation"
    fresh=not all(p.is_file() for p in (legacy/"prepared01/prepared.json",legacy/"run01/thermal.json",legacy/"build01/build.json",
                                        legacy/"run01/steady05/restart_000199.csv"))
    driver=ROOT/"outputs/v333_validation/build02"
    cached_programs=all((driver/name).is_file() for name in ("build.json","physical.exe","release.exe"))
    archived_toolchain=all((ROOT/path).is_file() for path in ("outputs/v317_validation/build01/wire_process.json",
            "tools/runtime/su2-thermal-build/build-v117/build.ninja","tools/runtime/su2-thermal-build/su2.zip"))
    if (completed is None and (fresh or not (cached_programs or archived_toolchain))) or (
            completed is not None and (completed.resolve().parent/"construction.json").is_file()):
        from shortSupportFirstRun import run
        run(out,completed)
        return
    if completed is None:
        if not (driver/"build.json").is_file():
            driver=out/"driver";build(driver)
        work=out/"run";setup(work,driver);run_steps(work)
    else:
        work=completed.resolve();work.relative_to(ROOT/"outputs")
        from energy import endpoint
        from physical_step import inspect_step
        j=read(work/"job.json");assert j["baseFrame"]==BASE_FRAME and j["seed"]==sha(work/"seed.csv")
        for branch,step in (("full",1),("half",1),("half",2)):
            old=work/"seed.csv" if step==1 else work/"half/step1/restart_00001.csv"
            frame=inspect_step(work,branch,step,old,False);inspect_step(work,branch,step,frame,True)
            endpoint(work,branch,step)
    result=export(work,out/"export")
    write_json(out/"completed.json",dict(pairSHA256=sha(out/"export/pair.json"),fullFrame=result["fullFrame"],halfFrame=result["halfFrame"],
        source=str(work),faces=result["faces"],areaM2=result["areaM2"],elapsedS=result["elapsedS"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    p.add_argument("--completed",type=Path,help="重读已通过的同一初场整步和两半步")
    a=p.parse_args();main(a.output.resolve(),a.completed)
