"""从 20D 自身通过验收的稳态帧推进同刻试验，输出两种步长的 32,480 面。"""
import argparse
from pathlib import Path
import sys

from outer_grid import HERE, ROOT, PREP_KEY, read, sha, write_json
sys.path.insert(0,str(HERE.parent/"v333"))
from build_time import build
from physical_step import setup, run_steps
from map_surface import mapped
from time_surface import export


def physical(out):
    prep=out/"prepared";steady=out/"steady";shared=out/"driver";driver=out/"time_driver";work=out/"time"
    rec=read(out/"steady_complete.json");grid=read(prep/"extension.json")
    assert grid["basePrepared"]==PREP_KEY and grid["shellRadiiD"]==[16,18,20]
    assert rec["prepared"]==grid["prepared"]==sha(prep/"prepared.json")
    assert rec["mesh"]==grid["mesh"]==sha(prep/"mesh.su2")
    assert rec["thermal"]==sha(steady/"thermal.json") and read(steady/"thermal.json")["passed"]
    shared_meta=read(shared/"build.json");assert shared_meta["preparedSHA256"]==rec["prepared"] and shared_meta["meshSHA256"]==rec["mesh"]
    # 一次任务只建一套驱动和一套初始扰动。异常退出后单独读回原端点再决定是否接续。
    build(driver,steady_build=shared,expected_build=sha(shared/"build.json"))
    conv=read(steady/"converged.json");base=conv["frame"];assert base==rec["frame"]
    j=setup(work,driver,prepared=prep,prepared_sha=rec["prepared"],steady=steady/conv["stage"],base_frame=base,steady_build=shared)
    assert j["changed"]==414 and abs(j["peakDeltaOverU"]-.01)<1e-14
    # 新增三层的每层径向边和外圈边来自实际表面连接。
    run_steps(work,prepared=prep,plan=HERE/"PLAN.md",expected_edges=921290+grid["addedEdges"])
    publish(out)


def publish(out):
    prep=out/"prepared";steady=out/"steady";grid=read(prep/"extension.json");rec=read(out/"steady_complete.json")
    transfer=ROOT/"cases/short_sphere/supportedAngularTransfer65.npz"
    steady_report=mapped(steady,out/"steady_surface",prepared=prep,expected_mesh=rec["mesh"],transfer_file=transfer,
                         transfer_sha="63c49c77d06cd2c1814d97f68769c62c3401d52748ad9863f0382fb8408a8609")
    assert steady_report["frameSHA256"]==rec["frame"]
    result=export(out/"time",out/"export",prepared=prep,prepared_sha=grid["prepared"],base_frame=rec["frame"],
                  steady_surface=out/"steady_surface/surface.zip")
    write_json(out/"completed.json",dict(pairSHA256=sha(out/"export/pair.json"),fullFrame=result["fullFrame"],
               halfFrame=result["halfFrame"],source=str(out/"time"),faces=result["faces"],areaM2=result["areaM2"],elapsedS=result["elapsedS"]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("output",type=Path)
    args=parser.parse_args();out=args.output.resolve();out.relative_to(ROOT/"outputs");physical(out)
