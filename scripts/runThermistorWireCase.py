"""从已审核的同网 12 K B 珠均温解出发，计算长引线温度边界下的全表面响应。"""
import argparse
from pathlib import Path
import shutil

import numpy as np

from runThermistorMeasuredCold import exclusive
from thermistorBead import ENGINE_SHA, ROOT, read, sha, write_json
from thermistorWireCase import native, prepare


def run(job,work,build,*,resume,checkpoint=None):
    if not resume:
        work.mkdir(parents=True,exist_ok=False)
        if shutil.disk_usage(work).free<16*1024**3:raise OSError("需预留 16 GiB 计算空间。")
    with exclusive(work):
        prep=work/"prepared"
        if not prep.exists():prepare(job,prep)
        spec=read(prep/"profile.json");base=Path(spec["base"])
        assert base.parent==job and spec["baseFrame"]==read(job/"completed.json")["frameSHA256"]==sha(base/"restart.csv")
        assert all(sha(prep/n)==digest for n,digest in spec["inputs"].items())
        control=work/"uniform_control"
        if checkpoint:checkpoint("uniform_control")
        native(prep,build,control,base/"restart.csv",1,"uniform")
        official=base/"ledger/control";seal=read(official/"checked.json")
        assert seal["binarySHA256"]==ENGINE_SHA and seal["files"]["initial_flow.csv"]==spec["baseFrame"]
        assert all(sha(official/n)==h for n,h in seal["files"].items())
        assert (official/"case.cfg").read_text("ascii")+"MARKER_PYTHON_CUSTOM= (wall)\n"==(control/"case.cfg").read_text("ascii")
        states=[np.loadtxt(p/"restart.csv",delimiter=",",skiprows=1) for p in (official,control)]
        difference=float(np.max(abs(states[0]-states[1])));assert difference<2e-7,difference
        control_check=dict(maxStateDifference=difference,sourceFrame=spec["baseFrame"],officialFrame=seal["files"]["restart.csv"],
                           customFrame=sha(control/"restart.csv"))
        path=work/"control.json"
        if path.exists():assert read(path)==control_check
        else:write_json(path,control_check)
        del states
        start=prep/"initial_flow.csv";windows=[]
        for n in range(1,6):
            dest=work/f"steady{n:02d}"
            if checkpoint:checkpoint(dest.name)
            r=native(prep,build,dest,start,320,"profile")
            windows.append(max(r["last20"].values()));print(dest.name,windows[-1],flush=True)
            if r["passed"]:
                from thermistorWireSurface import collect
                collect(prep,build,dest,work,checkpoint)
                return
            if len(windows)>=3 and windows[-1]>=windows[-2]>=windows[-3]:
                raise RuntimeError("连续两段末窗未改善，已保存原生终态。")
            start=dest/"restart.csv"
        raise RuntimeError("达到本次预算，尚无满足原六分量末窗条件的终态。")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("cold_job","out","build"):p.add_argument(name,type=Path)
    p.add_argument("--resume",action="store_true")
    a=p.parse_args();paths=[v.resolve() for v in (a.cold_job,a.out,a.build)]
    for path in paths:path.relative_to(ROOT/"outputs")
    run(*paths,resume=a.resume)
