"""逐项核对主窗口冷启动现场、确认的体网和已发布的全表面参考解。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorReaction import terminal
from beadSurfaceResult import thermal_frame


def audit(folder):
    ui=read(folder/"workflow.json")
    supervisor=read(folder/"debugger_process.json")
    assert supervisor["status"]=="terminal" and supervisor["debuggerExitCode"]==0
    assert not supervisor["stopReason"] and not supervisor["survivingProcessIDs"] and all(supervisor["inputsStable"].values())
    assert ui["realSolve"]
    prepared=Path(ui["grid15"]);native=Path(ui["output"])/"native"
    expected=read(ROOT/"cases/thermistor_b/grids.json")["grids"]["outer"]
    processes={}
    for name in ("axis","sector","strong","probes/reaction0","probes/plus20","probes/reaction20"):
        directory=native/name;r=read(directory/"engine_process.json")
        if name!="axis":
            terminal(directory)
        assert r["status"]=="terminal" and r["exitCode"]==0 and not (r["cancelled"] or r["timedOut"] or r.get("supervisorError"))
        assert all(r["inputsStable"].values()) and sha(directory/"restart.csv")==r["outputSHA256"]["restart.csv"]
        resources=read(directory/"resources.json");assert not resources["resourceExceeded"]
        processes[name]=dict(pid=r["processID"],seconds=r["elapsedSeconds"],savedFrame=sha(directory/"restart.csv"),
                             receiptSHA256=sha(directory/"engine_process.json"),resources=resources)
    sector=read(native/"sector/engine_process.json")
    for name in ("mesh.su2","nodes.npz","prepared.json"):
        assert sha(prepared/name)==sector["inputSHA256"][str(prepared/name)]
    for name in ("mesh.su2","nodes.npz"):
        assert sha(prepared/name)==expected[name]==sha(native/"sector"/name)==sha(native/"strong"/name)
    old=ROOT/"outputs/v282_validation/result_outer01";new=native/"result"
    report=read(new/"summary.json");baseline=read(old/"summary.json")
    for name in ("bodySHA256","sourceFrameSHA256","sourceReactionSHA256","configSHA256","wallTableSHA256"):
        assert report[name]==baseline[name],name
    assert sha(new/"wall.csv")==report["wallTableSHA256"]
    assert report["last20MaxLogRMS"]==thermal_frame(native/"strong")[1]
    rows=np.genfromtxt(new/"wall.csv",delimiter=",",names=True)
    np.testing.assert_array_equal(rows["faceID"],np.arange(5216))
    power=rows["areaM2"]@rows["heatFluxWm2"]
    assert abs(power-ui["powerW"])<1e-12 and abs(power-report["powerW"])<1e-12
    return dict(preparedDirectory=str(prepared),nativeDirectory=str(native),preparedSHA256=sha(prepared/"prepared.json"),
                sectorMeshSHA256=sha(native/"sector/mesh.su2"),frozenReferenceSource=baseline["sourceFrameSHA256"],
                allFaceTableSHA256=report["wallTableSHA256"],faces=len(rows),areaM2=rows["areaM2"].sum(),
                powerW=power,strongTemperatureRangeK=report["temperatureRangeK"],last20MaxLogRMS=report["last20MaxLogRMS"],
                heatChange20=report["heatChange20"],boundaryRelativeDifference=report["boundaryRelativeDifference"],
                actualMatlabSeconds=supervisor["elapsedSeconds"],peakFamilyRSSBytes=supervisor["peakFamilyRSSBytes"],
                newEngineSeconds=sum(v["seconds"] for v in processes.values()),processes=processes)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("folder",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args()
    result=audit(args.folder)
    with args.output.open("x",encoding="utf8") as stream:
        json.dump(result,stream,indent=2);stream.write("\n")
