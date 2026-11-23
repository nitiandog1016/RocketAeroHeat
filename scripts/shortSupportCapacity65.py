"""用登记的 Short Test 65 双支管混合网从冷态推进，自动核查收敛及端点热量。"""
import argparse
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"docs/validation/v328"))
from build import build
from endpoint import observer_build, sample
from run import BASE, PREP, read, run, sha, write_json
from verify import budget, table


def check_endpoints(samples):
    """按两份真实端点作独立反力及逐边复核，给出本阶段的科研资格。"""
    r=read(samples/"samples.json");assert r["early"]["start"]==r["control"]["start"]==r["plus20"]["start"]
    assert r["late"]["start"]==r["plus20"]["end"]
    for name,rec in r.items():
        stage=samples/name
        assert read(stage/"sample.json")==rec and sha(stage/"initial_flow.csv")==rec["start"]
        assert sha(stage/"restart.csv")==rec["end"] and sha(stage/"history.csv")==rec["history"]
        terminal=read(stage/"debugger_process.json")
        assert terminal["status"]=="terminal" and terminal["debuggerExitCode"]==0 and all(terminal["inputsStable"].values())
        assert terminal["inputSHA256"][terminal["command"][0]]==rec["exe"]
    assert np.max(abs(table(samples/"early","restart.csv")-table(samples/"control","restart.csv")))<2e-12
    s=read(PREP/"prepared.json")["condition"]
    scale=s["rho"]*s["speedMps"]*s["cp"]*s["diameterM"]**2*s["scalarUnitK"]
    with np.load(PREP/"mapping.npz") as mesh:wet=mesh["wet"]
    early=budget(samples/"early",samples/"releaseEarly",wet,scale)
    late=budget(samples/"late",samples/"releaseLate",wet,scale)
    change=abs(late["wallW"]/early["wallW"]-1)
    stable=change<.002 and max(r["plus20"]["window"])<-8
    passed=max(abs(v["physicalGapFraction"]) for v in (early,late))<.005
    return dict(early=early,late=late,extra20Change=change,thermalStable=bool(stable),balancePass=bool(passed))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory",type=Path)
    parser.add_argument("--read-samples",action="store_true",help="仅读取已有的端点目录并核算热量")
    args=parser.parse_args();folder=args.directory.resolve()
    folder.relative_to(ROOT/"outputs")
    if args.read_samples:
        result=check_endpoints(folder)
    else:
        folder.mkdir(parents=True,exist_ok=False)
        production=BASE/"build01"
        if not (production/"build.json").is_file():
            production=folder/"driver";build(production)
        work=folder/"solve";run(work,production)
        decision=read(work/"decision.json")
        if decision["reason"]!="residuals-passed":
            sys.exit("稳态阶段停止，详见 "+str(work/"decision.json"))
        observer=folder/"observer";observer_build(observer,production)
        out=folder/"thermal";sample(out,work/f"steady{decision['stages']:02d}",production,observer)
        result=check_endpoints(out);write_json(folder/"thermal_check.json",result)
    for name in ("early","late"):
        v=result[name]
        print(name,"wall W",v["wallW"],"far W",v["scalarFarW"],"gap %",100*v["physicalGapFraction"])
    print("thermalStable",result["thermalStable"],"balancePass",result["balancePass"])
