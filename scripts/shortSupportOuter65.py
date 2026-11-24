"""在已通过验收的 20D 双支管银球基流上，重新求出同刻的整步与两半步结果。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/"docs/validation/v335"
sys.path[:0]=[str(HERE),str(HERE.parent/"v333")]
from physical_step import setup,run_steps,inspect_step
from energy import endpoint
from time_surface import export
from build_time import read,sha,write_json


def run(job,base,completed=None):
    spec=read(base/"steady_complete.json");prep=base/"prepared";driver=base/"time_driver"
    grid=read(prep/"extension.json");pair_file=ROOT/"cases/short_sphere/outer20/supportedTime65.json"
    assert sha(pair_file)=="07df1aa77303c11e81d4ebd2cd92861501cd9d24e58df9fd999e4ab45ab576d3";pair=read(pair_file)
    assert grid["shellRadiiD"]==[16,18,20] and spec["mesh"]==pair["meshSHA256"]==sha(prep/"mesh.su2")
    assert spec["thermal"]==sha(base/"steady/thermal.json") and read(base/"steady/thermal.json")["passed"]
    assert spec["frame"]==pair["steadyFrame"] and spec["prepared"]==sha(prep/"prepared.json")
    baseline=read(base/"steady/converged.json");assert baseline["frame"]==spec["frame"]
    job.mkdir(parents=True,exist_ok=False);edges=921290+grid["addedEdges"]
    if completed is None:
        work=job/"time"
        j=setup(work,driver,prepared=prep,prepared_sha=spec["prepared"],steady=base/"steady"/baseline["stage"],
                base_frame=spec["frame"],steady_build=base/"driver")
        assert j["seed"]==pair["seed"] and j["changed"]==414
        run_steps(work,prepared=prep,plan=HERE/"PLAN.md",expected_edges=edges)
    else:
        work=completed.resolve();work.relative_to(ROOT/"outputs");j=read(work/"job.json")
        assert j["baseFrame"]==spec["frame"] and j["seed"]==sha(work/"seed.csv")==pair["seed"]
        assert j["preparedSHA256"]==spec["prepared"] and j["buildSHA256"]==sha(driver/"build.json")
        for branch,step in (("full",1),("half",1),("half",2)):
            old=work/"seed.csv" if step==1 else work/"half/step1/restart_00001.csv"
            frame=inspect_step(work,branch,step,old,False,prepared=prep)
            inspect_step(work,branch,step,frame,True,prepared=prep)
            endpoint(work,branch,step,prepared=prep,expected_edges=edges)
    result=export(work,job/"export",prepared=prep,prepared_sha=spec["prepared"],base_frame=spec["frame"],
                  steady_surface=base/"steady_surface/surface.zip")
    write_json(job/"completed.json",dict(pairSHA256=sha(job/"export/pair.json"),fullFrame=result["fullFrame"],
               halfFrame=result["halfFrame"],source=str(work),faces=result["faces"],areaM2=result["areaM2"],elapsedS=result["elapsedS"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    p.add_argument("--base",type=Path,default=ROOT/"outputs/v335_validation/outer01")
    p.add_argument("--completed",type=Path,help="复核已经完成的 20D 三个物理步，重新导出两种时间步结果")
    a=p.parse_args();out=a.output.resolve();out.relative_to(ROOT/"outputs")
    run(out,a.base.resolve(),a.completed)
