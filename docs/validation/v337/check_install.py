"""在 v336 独立安装的路径下复读真实三步，内存中加载本版入口，保留原安装所有源码。"""
import json
from pathlib import Path
import sys
from types import ModuleType

ROOT=Path(__file__).resolve().parents[3]
APP=ROOT/"outputs/v336_validation/clean01"
JOB=APP/"outputs/short65_time_6dea1a68-073d-4129-867d-3772b804b131"


def updated(name):
    source=ROOT/"scripts"/(name+".py")
    logical=APP/"scripts"/source.name
    module=ModuleType(name)
    module.__file__=str(logical)
    sys.modules[name]=module
    exec(compile(source.read_bytes(),str(logical),"exec"),module.__dict__)
    return module


if __name__=="__main__":
    sys.path[:0]=[str(APP/"scripts"),str(APP/"docs/validation/v333")]
    entry=updated("shortSupportTime65")
    updated("shortSupportStationGradient");fresh=updated("shortSupportFirstRun")
    assert entry.ROOT==fresh.ROOT==APP and (JOB/"construction.json").is_file()
    out=APP/"outputs"/sys.argv[1]
    assert out.parent==APP/"outputs" and not out.exists()
    entry.main(out,JOB/"run")
    receipt=json.loads((out/"completed.json").read_bytes())
    old=json.loads((JOB/"completed.json").read_bytes())
    for k,v in old.items():assert receipt[k]==v
    reference=ROOT/"outputs/v337_validation/reread01/export/stationGradients.json"
    assert (out/"export/stationGradients.json").read_bytes()==reference.read_bytes()
    assert entry.sha(reference)==receipt["stationGradientsSHA256"]
    print(receipt["stationGradientsSHA256"])
