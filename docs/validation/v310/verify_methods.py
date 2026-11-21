"""用量纲式（28）、椭球导热和 v308 原面表独立复核工作台的珠部比较。"""
import argparse
import csv
from hashlib import sha256
import io
import json
import math
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from phaseAwareDryAir import PhaseAwareDryAir


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_formula(m):
    """用旋转长椭球的导热形状因子重建极限 Nu，以量纲式的速度项直接求 h。"""
    e=.53;L=.3430e-3;area=3.094e-6;u=3
    major=1.;minor=major*math.sqrt(1-e**2)
    f=2*math.pi*minor**2*(1+major*math.asin(e)/(minor*e))
    volume=4*math.pi*major*minor**2/3
    capacity=4*math.pi*major*e/math.atanh(e)
    nu0=capacity*(2*volume/f)/f
    air=PhaseAwareDryAir();delta=m["wallTemperatureK"]-293.15;tm=293.15+delta/2
    state=air.evaluate(101325,tm);prop=air.transport(101325,tm)
    k=float(prop["conductivityWmK"]);nu=float(prop["dynamicViscosityPas"]/state["densityKgm3"])
    h0=k/L*nu0;hf=.29*k*math.sqrt(u/nu)/math.sqrt(L);h=h0+hf
    gr=9.80665*delta*L**3/(tm*nu**2);re=u*L/nu
    for key,value in (("meanHWm2K",h),("outwardW",h*area*delta),("reynoldsL",re),
                      ("grashofL",gr),("shapeNu0",nu0),("tableMinusFormulaW",(.664-nu0)*k/L*area*delta)):
        assert math.isclose(m[key],value,rel_tol=1e-12),(key,m[key],value)
    assert m["available"] and m["air"]["pressureCeilingPa"]>=101325 and m["filmTemperatureK"]==tm
    assert math.isclose(m["air"]["kinematicViscosityM2s"],nu,rel_tol=1e-12)
    return dict(filmTemperatureK=tm,reynoldsL=re,grashofL=gr,fourGrOverReSquared=4*gr/re**2,
                conductionHWm2K=h0,forcedIncrementHWm2K=hf,meanHWm2K=h,outwardW=h*area*delta,
                shapeNu0=nu0,tableNu0=.664,tableLastDigitMatches=round(nu0,3)==.664,
                tableMinusFormulaW=(.664-nu0)*k/L*area*delta)


def verify(folders,out):
    prior=HERE.parent/"v308/evidence";base=read(prior/"verification.json")
    table=prior/"wall299K.zip";assert digest(table)==base["tableSHA256"]
    with ZipFile(table) as package:
        faces=np.genfromtxt(io.StringIO(package.read("faces.csv").decode("utf-8-sig")),delimiter=",",names=True)
    area=faces["area_m2"];bead=np.isin(faces["region_id"],[3,4,5])
    assert len(area)==10272 and abs(np.sum(area)-3.19665061256752e-6)<1e-17
    checks={}
    for name,folder in zip(("medium","fine"),folders):
        ui=read(folder/"ui.json");run=read(folder/"debugger_process.json");side=read(folder/"heat.json")
        assert run["status"]=="terminal" and run["debuggerExitCode"]==0 and all(run["inputsStable"].values())
        assert not (run["survivingProcessIDs"] or run["supervisorError"] or run["stopReason"])
        raw=(folder/"workbench/case-data.js").read_text(encoding="utf8").strip()
        assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
        page=json.loads(raw[len("window.ROCKET_CASE="):-1]);m=ui["comparison"]
        assert m==side["meanConvection"]==page["thermistorLead"]["meanConvection"]
        assert m["doi"]=="10.1299/kikai1938.25.810" and m["scope"] and m["wallTemperatureK"]==299.15
        calc=check_formula(m);endpoint=check_formula(ui["endpoint"])
        assert not calc["tableLastDigitMatches"] and side["beadMeanComparison"]["measuredHWm2K"] is None
        q=faces[name+"_q_W_m2"];watts=float(-area[bead]@q[bead]);surface_h=watts/(m["areaM2"]*6)
        assert abs(m["surfaceOutwardW"]-watts)<1e-13 and abs(surface_h-m["surfaceMeanHWm2K"])<1e-10
        assert ui["sourceFrame"]==base[name]["source"] and ui["webglError"]==0 and ui["visibleSamples"]>0
        assert ui["shownRows"][2:]==[f'{m["meanHWm2K"]:.2f}',f'{surface_h:.2f}',f'{1000*m["outwardW"]:.4f}',f'{1000*watts:.4f}']
        fields={v["key"]:v for v in page["fields"]};assert ui["probeFace"]==207 and "实验偏差未知" in fields["heatFlux"]["stateNote"]
        for key,values in (("pressure",faces[name+"_p_Pa"]/1000),("heatFlux",q/1000),("nearWallCellTemperature",faces[name+"_first_K"])):
            np.testing.assert_allclose(fields[key]["values"],values,rtol=0,atol=1e-10)
        assert page["thermistorLead"]["jointExperimentVerifiedM2"]==0
        with (folder/"heat.csv").open(encoding="utf-8-sig",newline="") as f:
            rows=list(csv.DictReader(f))
        assert len(rows)==len(q) and digest(folder/"heat.csv")==side["csvSHA256"]
        for j,row in enumerate(rows):
            assert int(row["face_id"])==j+1 and row["field_key"]=="heatFlux" and row["unit"]=="kW/m²"
            assert row["availability_code"]=="2" and row["qualification_code"]=="0" and abs(float(row["value"])*1000-q[j])<1e-9
        checks[name]=dict(formula=calc,endpoint=endpoint,surfaceMeanHWm2K=surface_h,surfaceOutwardW=watts,
                          relativeToFormula=watts/m["outwardW"]-1,printedTableNuMismatch=calc["shapeNu0"]-.664,
                          sourceFrame=ui["sourceFrame"],faces=len(rows),jointExperimentVerifiedM2=0,
                          heatCSV=digest(folder/"heat.csv"),sidecar=digest(folder/"heat.json"),
                          pageData=digest(folder/"workbench/case-data.js"),seconds=run["elapsedSeconds"],peakRSS=run["peakFamilyRSSBytes"])
    assert checks["medium"]["formula"]==checks["fine"]["formula"]
    out.mkdir(parents=True,exist_ok=False)
    (out/"verification.json").write_bytes((json.dumps(checks,indent=2,allow_nan=False)+"\n").encode("utf8"))
    with ZipFile(out/"ui_receipts.zip","x",compression=ZIP_DEFLATED) as package:
        for folder in folders:
            for key in ("ui.json","debugger_process.json"):
                package.write(folder/key,folder.name+"/"+key)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("medium","fine","output"):parser.add_argument(name,type=Path)
    args=parser.parse_args();verify((args.medium,args.fine),args.output)
