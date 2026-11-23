"""由工作台从已核对的银球基流启动尾迹试验，也可从本次已验收的物理端点接续。"""
import argparse
from pathlib import Path
from uuid import uuid4

from buildThermistorWireDriver import build as build_driver
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive
from shortSphere65 import run as steady_run
from shortSphereWake import BASE, CASE, DRIVER, DT, ROOT, SAMPLES, advance, build_observer, read, setup, sha
from shortSphereWakeReadout import readout
from shortSphereWakeResults import compare


def run(work,through,base=None,driver=None,resume=False):
    assert through in SAMPLES
    if not resume:
        if work.exists():raise FileExistsError("请使用新目录，新算会保留既有任务。")
        if driver is None:
            driver=ROOT/"outputs/v317_validation/build01"
            if not (driver/"wire.exe").is_file() or sha(driver/"wire.exe")!=DRIVER["wire"]:
                driver=work.with_name(work.name+"_driver");build_driver(driver)
        if base is None:
            base=ROOT/"outputs/v321_validation/run01"
            if not (base/"endpoints.json").is_file():
                base=work.with_name(work.name+"_base");steady_run(base,driver)
        setup(work,base,driver)
    with exclusive(work):
        if not (work/"observer/build.json").exists():build_observer(work)
        for t in SAMPLES:
            if t>through:break
            for b,dt in DT.items():
                advance(work,b,t);readout(work,b,round(t/dt))
        name=f"result_{round(through*8):03d}";out=work/name;receipt=work/(name+".json")
        if not out.exists():
            # 读数或打包中断时保留当次文件，只有完整目录才会成为正式结果。
            draft=work/(name+"_assembling_"+uuid4().hex)
            compare(work,through,draft)
            manifest=dict(timeDU=through,base=BASE,seed=read(work/"job.json")["seedSHA256"],case=sha(CASE),
                          comparison=sha(draft/"comparison.json"),archives={b:sha(draft/("wake_"+b+".zip")) for b in DT},
                          summary={b:sha(draft/b/"summary.json") for b in DT})
            atomic_write_json(draft/"manifest.json",manifest)
            draft.rename(out)
        if not receipt.exists():
            atomic_write_json(receipt,read(out/"manifest.json"))
        ack=read(receipt);assert ack["timeDU"]==through and ack["base"]==BASE and ack["case"]==sha(CASE)
        assert ack["seed"]==sha(work/"seed.csv") and ack["comparison"]==sha(out/"comparison.json")
        assert all(sha(out/b/"summary.json")==ack["summary"][b] and sha(out/("wake_"+b+".zip"))==ack["archives"][b] for b in DT)
        for b in DT:
            metadata=read(out/("wake_"+b+".json"))
            assert (metadata["archiveSHA256"],metadata["summarySHA256"],metadata["caseSHA256"])==(ack["archives"][b],ack["summary"][b],ack["case"])
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("--through",type=float,default=8.)
    p.add_argument("--base",type=Path);p.add_argument("--driver",type=Path);p.add_argument("--resume",action="store_true")
    a=p.parse_args();work=a.work.resolve();work.relative_to(ROOT/"outputs")
    run(work,a.through,a.base.resolve() if a.base else None,a.driver.resolve() if a.driver else None,a.resume)
