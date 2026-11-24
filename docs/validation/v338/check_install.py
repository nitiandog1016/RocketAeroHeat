"""从 v336 独立安装的实际保存帧复读，保留原始源码与原任务。"""
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
import importlib.util
spec=importlib.util.spec_from_file_location("short337check",HERE.parent/"v337/check_install.py")
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)

if __name__=="__main__":
    app=prior.APP;job=prior.JOB
    sys.path[:0]=[str(app/"scripts"),str(app/"docs/validation/v333")]
    entry=prior.updated("shortSupportTime65")
    prior.updated("shortSupportStationGradient");prior.updated("shortSupportAzimuth")
    first=prior.updated("shortSupportFirstRun")
    assert entry.ROOT==first.ROOT==app
    out=app/"outputs"/sys.argv[1];assert out.parent==app/"outputs" and not out.exists()
    entry.main(out,job/"run")
    receipt=json.loads((out/"completed.json").read_bytes())
    old=json.loads((job/"completed.json").read_bytes())
    assert all(receipt[k]==v for k,v in old.items())
    assert entry.sha(out/"export/stationAzimuth.zip")==receipt["stationAzimuthSHA256"]==entry.sha(ROOT/"cases/short_sphere/stationAzimuth65.zip")
    print(receipt["stationAzimuthSHA256"])
