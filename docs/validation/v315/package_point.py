"""核验 12 K 名义侧吹预测及独立原生热边界，在试件数据包内记录风洞测点和适用范围。"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];CASES=ROOT/"cases/thermistor_b"
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorReaction import terminal


def package(folder,six,dest):
    source=read(CASES/"sidewindAngular96.json");rec=terminal(folder);s=read(folder/"surface/summary.json");c=read(folder/"case.json")
    assert sha(folder/"restart.csv")==s["sourceFrameSHA256"]==rec["outputSHA256"]["restart.csv"]==s["energy"]["inputFrameSHA256"]
    assert sha(folder/"surface/wall.csv")==s["wallSHA256"] and sha(folder/"case.json")==s["caseSHA256"]
    for n,k in (("mesh.su2","meshSHA256"),("nodes.npz","nodeMapSHA256"),("body.json","bodySHA256")):
        assert sha(folder/n)==s[k]==source[k]
    assert terminal(folder/"reaction0")["outputSHA256"]["restart.csv"]==s["sourceReactionSHA256"]
    assert terminal(folder/"plus20")["outputSHA256"]["restart.csv"]==s["plus20FrameSHA256"]
    assert terminal(folder/"reaction20")["outputSHA256"]["restart.csv"]==s["reaction20SHA256"]
    assert max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8 and s["relativeChange20"]<.002
    assert s["energy"]==read(folder/"ledger/energy.json") and abs(s["energy"]["relativeDifference"])<.005
    assert s["energy"]["temperatureScaleK"]==c["scalarReferenceDeltaK"]==24 and c["wallK"]-c["airK"]==12
    assert sha(CASES/"sidewindAngular96.zip")==source["archiveSHA256"]
    def parse(data):return np.genfromtxt(io.BytesIO(data),delimiter=",",names=True)
    with ZipFile(CASES/"sidewindAngular96.zip") as z:
        previous=parse(z.read("wall.csv"));original=json.loads(z.read("summary.json"))
    base=read(six/"surface/summary.json");control=terminal(six)
    assert control["outputSHA256"]["restart.csv"]==base["sourceFrameSHA256"]==base["energy"]["inputFrameSHA256"]
    assert sha(six/"initial_flow.csv")==original["sourceFrameSHA256"]==c["sourceSixKFrame"]
    assert read(folder/"field_audit.json")["steps"]==read(six/"field_audit.json")["steps"]==40
    assert sha(six/"case.json")==base["caseSHA256"]==source["caseSHA256"] and sha(six/"nodes.npz")==s["nodeMapSHA256"]
    assert sha(six/"surface/wall.csv")==base["wallSHA256"] and base["meshSHA256"]==s["meshSHA256"]
    for path,name in (("reaction0","sourceReactionSHA256"),("plus20","plus20FrameSHA256"),("reaction20","reaction20SHA256")):
        assert terminal(six/path)["outputSHA256"]["restart.csv"]==base[name]
    assert base["energy"]["temperatureScaleK"]==6 and base["relativeChange20"]<.002 and abs(base["energy"]["relativeDifference"])<.005
    assert max(base["residualWindow"].values())<-8 and max(base["plus20Window"].values())<-8
    ref=parse((six/"surface/wall.csv").read_bytes())
    native=parse((folder/"surface/wall.csv").read_bytes())
    for k in ("faceID","regionID","areaM2"):np.testing.assert_array_equal(native[k],ref[k])
    def compare(table,summary):
        predicted=2*table["heatFluxWm2"];delta=native["heatFluxWm2"]-predicted
        return dict(maxHeatWm2=float(max(abs(delta))),maxRelativeHeat=float(max(abs(delta/predicted))),
                    maxPressurePa=float(max(abs(native["pressurePa"]-table["pressurePa"]))),
                    maxFirstK=float(max(abs(native["nearWallCellTemperatureK"]-(c["airK"]+2*(table["nearWallCellTemperatureK"]-c["airK"]))))),
                    relativePower=abs(s["totalOutwardW"]/(2*summary["totalOutwardW"])-1))
    pcheck=compare(ref,base);unmatched=compare(previous,original)
    assert pcheck["maxRelativeHeat"]<.0005 and pcheck["relativePower"]<.0005 and pcheck["maxPressurePa"]<.001 and pcheck["maxFirstK"]<.001,pcheck
    reading=read(HERE.parent/"v282/evidence/comparison.json")["experiment"];gap=read(CASES/"beadExperimentGap.json")
    assert reading["doi"]==gap["doi"] and reading["specimen"]=="B" and reading["speedMps"]==c["speedMps"]
    intervals=[f["points"][1]["readRangeK"] for f in gap["figures"]];assert all(lo<=12<=hi for lo,hi in intervals)
    exp=dict(doi=gap["doi"],pdfSHA256=gap["pdfSHA256"],printedPage=815,figure=7,speedMps=3,nominalDeltaTK=12,
             figure7DeltaTK=reading["deltaTK"],figure8DeltaTK=gap["figures"][1]["points"][1]["deltaTK"],deltaTReadRangesK=intervals,
             hWm2K=reading["hWm2K"],graphReadRangeWm2K=reading["graphReadRangeWm2K"],reportedPointUncertaintyWm2K=None,
             areaM2=gap["nominalBeadAreaM2"],measuredAxisAngleDeg=None,
             mounting="第 813–814 页，水平气流与直立支承导线。名义共轴模型采用侧吹，B 珠长轴实测夹角未知。",
             heatDefinition="电阻推温；总电热扣除辐射及两根铂铑支线导热。仅对照珠部面积平均换热。",
             temperatureDefinition="第 818–819 页附录讨论内部电阻温度与表面温度的差异。均匀壁温为计算近似。",
             readingSHA256=sha(HERE.parent/"v282/evidence/comparison.json"),gapSHA256=sha(CASES/"beadExperimentGap.json"))
    context=dict(experiment=exp,nativeCheck=pcheck,unmatched=unmatched,sixKPressureDriftPa=float(max(abs(ref["pressurePa"]-previous["pressurePa"]))),
                 sixKFrame=base["sourceFrameSHA256"],sixKTable=base["wallSHA256"],
                 checkSeedFrame=c["sourceSixKFrame"],twelveKFrame=s["sourceFrameSHA256"],jointVerifiedM2=0)
    blob=(json.dumps(context,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode("utf8")
    dest.mkdir(parents=True,exist_ok=True);archive=dest/"windPoint12.zip";manifest=archive.with_suffix(".json")
    assert not manifest.exists()
    with ZipFile(archive,"x") as z:
        for n,path in (("body.json",folder/"body.json"),("case.json",folder/"case.json"),("wall.csv",folder/"surface/wall.csv"),
                       ("summary.json",folder/"surface/summary.json"),("six_wall.csv",six/"surface/wall.csv"),
                       ("six_case.json",six/"case.json"),("six_summary.json",six/"surface/summary.json"),("comparison.json",None)):
            z.writestr(ZipInfo(n,(2026,1,1,0,0,0)),blob if path is None else path.read_bytes(),compress_type=ZIP_DEFLATED)
    key={k:s[k] for k in ("sourceFrameSHA256","sourceReactionSHA256","caseSHA256","bodySHA256","meshSHA256","nodeMapSHA256","wallSHA256")}
    key.update(archiveSHA256=sha(archive),comparisonSHA256=hashlib.sha256(blob).hexdigest(),sixKFrame=base["sourceFrameSHA256"],
               sixKWallSHA256=base["wallSHA256"],sixKReaction=base["sourceReactionSHA256"])
    manifest.write_bytes((json.dumps(key,indent=2)+"\n").encode("utf8"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("run","six","dest"):p.add_argument(n,type=Path)
    a=p.parse_args();package(a.run.resolve(),a.six.resolve(),a.dest.resolve())
