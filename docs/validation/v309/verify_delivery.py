"""按 v308 发布的原面表复核两档三维页、热流 CSV 与珠部实验说明。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(path):
    return json.loads(path.read_text(encoding="utf8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plot_intervals(source,out):
    """用原图点心及读图范围画温差轴，虚线表示本次给定的工作温差。"""
    fig,ax=plt.subplots(figsize=(7.6,2.6))
    for row,g in enumerate(source["figures"]):
        y=1-row
        for point in g["points"]:
            x=point["deltaTK"];lo,hi=point["readRangeK"]
            ax.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt="o",capsize=5,color="#175777")
            ax.annotate(f"{x:.2f} K",(x,y),xytext=(0,9),textcoords="offset points",ha="center",fontsize=9)
    ax.axvline(6,color="#ad4d2d",ls="--")
    ax.annotate("6 K",(6,.5),xytext=(5,0),textcoords="offset points",color="#ad4d2d")
    ax.set(yticks=[1,0],yticklabels=["Fig. 7 (linear)","Fig. 8 (log)"],xlim=(1,14),ylim=(-.38,1.55),
           xlabel="Bead-to-air temperature difference (K)",title="Specimen B at 3 m/s: adjacent plotted measurements")
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="x",alpha=.15)
    fig.tight_layout();fig.savefig(out,dpi=160);plt.close(fig)


def audit(folder,k,base,source,earlier):
    receipt=read(folder/"debugger_process.json");ui=read(folder/"ui.json");side=read(folder/"heat.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not (receipt["survivingProcessIDs"] or receipt["stopReason"] or receipt["supervisorError"])
    page_path=folder/"workbench/case-data.js";raw=page_path.read_text(encoding="utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    page=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in page["fields"]}
    obs=ui["comparison"]
    assert obs==page["thermistorLead"]["beadMeanComparison"]==side["beadMeanComparison"]
    assert obs["sourcePDFSHA256"]==source["pdfSHA256"] and obs["figures"]==[7,8] and obs["deltaTK"]==6
    assert obs["status"]=="noMatchingPoint" and all(obs[f] is None for f in ("measuredHWm2K","relativeDifference","reportedPointUncertaintyWm2K"))
    assert ui["heatNote"]==fields["heatFlux"]["stateNote"] and "实验偏差未知" in ui["heatNote"]
    assert ui["sourceFrame"]==earlier["source"] and page["resultCondition"]["wallUniformValueK"]==299.15
    heat=base[k+"_q_W_m2"]
    np.testing.assert_allclose(fields["heatFlux"]["values"],heat/1000,atol=1e-10,rtol=0)
    np.testing.assert_allclose(fields["pressure"]["values"],base[k+"_p_Pa"]/1000,atol=1e-10,rtol=0)
    np.testing.assert_allclose(fields["nearWallCellTemperature"]["values"],base[k+"_first_K"],atol=1e-10,rtol=0)
    np.testing.assert_array_equal(base["face_id"],np.arange(1,10273))
    area=base["area_m2"];bead=np.isin(base["region_id"],[3,4,5])
    power=float(-area[bead]@heat[bead]);h=power/(6*source["nominalBeadAreaM2"])
    assert abs(h-obs["predictedHWm2K"])<1e-10 and abs(obs["beadAndJointsW"]-power)<1e-13
    assert abs(float(-area@heat)-ui["wholeOutwardW"])<1e-13 and ui["faces"]==10272 and ui["chosenFace"]==207
    assert ui["visibleSamples"]>0 and ui["webglError"]==0 and page["thermistorLead"]["jointExperimentVerifiedM2"]==0
    with (folder/"heat.csv").open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==len(area) and digest(folder/"heat.csv")==side["csvSHA256"]
    for i,row in enumerate(rows):
        assert int(row["face_id"])==i+1 and row["availability_code"]=="2" and row["qualification_code"]=="0"
        assert row["field_key"]=="heatFlux" and row["unit"]=="kW/m²"
        assert abs(float(row["value"])*1000-heat[i])<1e-9 and abs(float(row["area_m2"])-area[i])<1e-18
    return dict(faces=ui["faces"],areaM2=float(sum(area)),beadAndJointsW=power,predictedHWm2K=h,
                wholeOutwardW=ui["wholeOutwardW"],matchedExperimentHWm2K=None,relativeDifference=None,
                sourceFrame=ui["sourceFrame"],heatCSV=digest(folder/"heat.csv"),sidecar=digest(folder/"heat.json"),
                pageData=digest(page_path),elapsedSeconds=receipt["elapsedSeconds"],peakFamilyRSSBytes=receipt["peakFamilyRSSBytes"])


def verify(dirs,out):
    sourcepath=ROOT/"cases/thermistor_b/beadExperimentGap.json"
    source=read(sourcepath)
    previous=HERE.parent/"v308/evidence";proof=read(previous/"verification.json");table=previous/"wall299K.zip"
    assert digest(table)==proof["tableSHA256"]
    with ZipFile(table) as z:
        base=np.genfromtxt(io.StringIO(z.read("faces.csv").decode("utf-8-sig")),delimiter=",",names=True)
    checks={k:audit(folder,k,base,source,proof[k]) for k,folder in zip(("medium","fine"),dirs)}
    out.mkdir(parents=True,exist_ok=False)
    (out/"verification.json").write_bytes((json.dumps(dict(calibrationSHA256=digest(sourcepath),**checks),indent=2)+"\n").encode("utf8"))
    plot_intervals(source,out/"temperature_intervals.png")
    with ZipFile(out/"ui_receipts.zip","x",compression=ZIP_DEFLATED) as z:
        for folder in dirs:
            for name in ("debugger_process.json","ui.json"):
                z.write(folder/name,folder.name+"/"+name)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("medium","fine","output"):parser.add_argument(name,type=Path)
    args=parser.parse_args();verify((args.medium,args.fine),args.output)
