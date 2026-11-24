"""组织 Short 双支管试件外域对照，所有新增场保存在独立任务目录。"""
import argparse
from pathlib import Path
import sys

from outer_grid import HERE, ROOT, RADII, PREP_KEY, read, sha, write_json
from build_shared import build, invoke
from run_shared import solve


def steady(out):
    out.mkdir(parents=True,exist_ok=False)
    prepared=out/"prepared"
    invoke("prepare",[sys.executable,HERE/"outer_grid.py",prepared],ROOT,out,180,3)
    extension=read(prepared/"extension.json")
    assert extension["basePrepared"]==PREP_KEY and extension["shellRadiiD"]==list(RADII)
    assert extension["prepared"]==sha(prepared/"prepared.json")
    driver=out/"driver"
    build(driver,prepared)
    solve(out/"steady",driver,prepared=prepared,plan=HERE/"PLAN.md")
    accepted=read(out/"steady/thermal.json")
    assert accepted["passed"] and accepted["meshSHA256"]==extension["mesh"]
    write_json(out/"steady_complete.json",dict(frame=accepted["early"]["sourceFrame"],prepared=extension["prepared"],
                                               mesh=extension["mesh"],thermal=sha(out/"steady/thermal.json")))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    a=p.parse_args();target=a.output.resolve();target.relative_to(ROOT/"outputs")
    steady(target)
