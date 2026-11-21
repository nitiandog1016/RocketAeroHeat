"""将两网的实际壁压和热流按同一带线原面配对，并保存小型原生回执和发布输入。"""
import json
from pathlib import Path
import shutil
import sys
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha


def save(path,value):
    with path.open("x",encoding="utf8") as f:
        json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False)


def compare():
    proof=HERE/"evidence"
    samples=[]
    with zipfile.ZipFile(proof/"native_receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as records:
        for grid in ("base","fine"):
            inputs=read(proof/(grid+"_inputs.json"))
            folder=ROOT/inputs["path"]
            assert all(sha(folder/n)==v for n,v in inputs["inputs"].items())
            s=read(folder/"surface/summary.json")
            assert s["sourceFrameSHA256"]==sha(folder/"restart.csv") and s["wallSHA256"]==sha(folder/"surface/wall.csv")
            wall=np.genfromtxt(str(folder/"surface/wall.csv"),delimiter=",",names=True)
            assert np.array_equal(wall["faceID"],np.arange(s["faces"]))
            samples.append((s,wall))
            for path in folder.rglob("*.json"):
                if path.name=="body.json":continue
                records.write(path,grid+"/"+path.relative_to(folder).as_posix())
            for path in (folder/"case.cfg",folder/"history.csv",folder/"engine_stdout.log"):
                records.write(path,grid+"/"+path.name)
            save(proof/(grid+"_surface.json"),s)
    (coarse,a),(fine,b)=samples
    for key in ("faceID","regionID","areaM2"):
        assert np.array_equal(a[key],b[key])
    assert coarse["bodySHA256"]==fine["bodySHA256"] and coarse["meshSHA256"]!=fine["meshSHA256"]
    smooth=read(ROOT/"docs/validation/v283/evidence/energy.json")["wallReleaseW"]
    result=dict(regions=[],relativeTotalPower=fine["totalOutwardW"]/coarse["totalOutwardW"]-1,
                relativeBeadPower=fine["beadAndJointsW"]/coarse["beadAndJointsW"]-1,
                relativeMeanHToExperiment=fine["beadMeanHWm2K"]/fine["experimentMeanHWm2K"]-1,
                previousSmoothBeadW=smooth,relativeBeadToPreviousSmooth=fine["beadAndJointsW"]/smooth-1)
    for region in fine["regions"]:
        pick=b["regionID"]==region["regionID"]
        row=dict(regionID=region["regionID"],name=region["name"],faces=int(pick.sum()),
                 relativePower=region["outwardW"]/coarse["regions"][region["regionID"]-1]["outwardW"]-1)
        for key,label in (("heatFluxWm2","HeatFluxWm2"),("pressurePa","PressurePa")):
            delta=b[key][pick]-a[key][pick]
            row["maxAbs"+label]=float(np.max(abs(delta)))
            row["rms"+label]=float(np.sqrt(np.average(delta**2,weights=b["areaM2"][pick])))
        result["regions"].append(row)
    result["withinThreePercentPower"]=max(abs(result["relativeTotalPower"]),abs(result["relativeBeadPower"]))<.03
    save(proof/"comparison.json",result)
    # 归档保留完整逐面表，工作台新算时从独立目录读取。字段与发布材料逐字节绑定。
    selected=ROOT/"outputs/v297_validation/fine01"
    target=ROOT/"cases/thermistor_b/leads3.zip"
    with (selected/"surface/surface.zip").open("rb") as src,target.open("xb") as dst:
        shutil.copyfileobj(src,dst)
    g=read(selected/"geometry.json")
    contract={k:g[k] for k in ("nodes","cells","meshSHA256","bodySHA256","nodeMapSHA256")}
    contract.update(initialSHA256=sha(selected/"initial_flow.csv"),configSHA256=sha(selected/"case.cfg"),
                    caseSHA256=sha(selected/"case.json"),archiveSHA256=sha(target))
    save(target.with_suffix(".json"),contract)


if __name__=="__main__":
    compare()
