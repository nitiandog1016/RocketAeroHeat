"""比较实窗、离线 HTML 表格、CSV 和已发布的节点剖面，并复核取样失败的行为。"""
import argparse
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import sys

import numpy as np
from verify import CASE, ROOT
sys.path.insert(0,str(ROOT/"scripts"))
from shortSupportStationGradient import derivative, profile, Unavailable


class StationTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active=False;self.cell=False;self.cells=[];self.download=None

    def handle_starttag(self,tag,attrs):
        attributes=dict(attrs)
        if tag=="details" and attributes.get("id")=="short65-stations":self.active=True
        if self.active and tag=="td":self.cell=True
        if self.active and tag=="a":self.download=attributes["href"]

    def handle_endtag(self,tag):
        if tag=="td":self.cell=False
        if self.active and tag=="details":self.active=False

    def handle_data(self,data):
        if self.cell:self.cells.append(float(data))


def check(ui):
    source=json.loads((CASE/"stationGradients65.json").read_bytes())
    half=next(g for g in source["sets"] if g["label"]=="half")
    page=ui/"workbench"
    parsed=StationTable();parsed.feed((page/"index.html").read_text("utf8"))
    values=np.array(parsed.cells).reshape(7,4)
    assert parsed.download=="Short65-stations-"+half["frameSHA256"]+".csv"
    csv_file=page/parsed.download
    assert csv_file.read_bytes()==(ui/"tab.csv").read_bytes()
    with csv_file.open(encoding="utf-8-sig",newline="") as f:table=list(csv.DictReader(f))
    assert len(table)==7
    for i,(r,item) in enumerate(zip(half["rows"],table)):
        np.testing.assert_allclose(values[i],[r[k] for k in ("angleDeg","experimentNu","reactionNu","gradientNu")],rtol=0,atol=5.0001e-5)
        exp=dict(gradientInWm2=r["gradientInWm2"],fluxInWm2=r["reactionInWm2"],reportedNu=r["experimentNu"],
                 computedNu=r["reactionNu"],gradientNu=r["gradientNu"],fluidGradientKPerM=r["gradientInWm2"]/half["kFilm"],
                 secondDistanceM=r["samples"][2]["distanceM"],firstSampleK=r["samples"][1]["temperatureK"],
                 secondSampleK=r["samples"][2]["temperatureK"],firstDistanceM=r["firstDistanceM"],fileFaceID=r["faceID"])
        for key,ref in exp.items():np.testing.assert_allclose(float(item[key]),ref,rtol=1e-13,atol=1e-15)
        assert item["gradientFrameSHA256"]==half["frameSHA256"] and item["gradientStatus"]=="已重建"
    # 用不等距解析二次场同时检查升温、降温、零热流和宽模板；断开邻接的测站要报告缺项。
    h=np.array([0.,1.2,2.7,4.,6.3])*1e-5
    for slope in (0.,-75000.,50000.):
        for curvature in (0.,3e8):
            t=350+slope*h+curvature*h*h
            for a,b in ((1,2),(2,4)):assert abs(derivative(h,t,a,b)-slope)<2e-8
    grid=dict(xyz=np.array([[.5,0,-.01],[.5,0,.01]]),sphere={0,1},first={0:2,1:3},adjacent={(0,1):{}},
              spec={"geometry":{"nearWallIntervals":5}})
    try:
        profile(grid,np.zeros(2),dict(angleDeg=0,leftNode=0,rightNode=1,rightWeight=.5),dict(diameterM=.0127,airK=300,scalarUnitK=90))
    except Unavailable as e:missing=str(e)
    else:raise AssertionError("断开的取样层被采用")
    result=json.loads((ui/"ui.json").read_bytes())
    proc=json.loads((ui/"window_process.json").read_bytes())
    assert not (proc["stopReason"] or proc["livePids"] or proc["exitCode"])
    assert result["damagedRejected"] and result["archiveSets"]==5 and result["frame"]==half["frameSHA256"]
    result.update(seconds=proc["elapsedSeconds"],peakRSSBytes=proc["peakSampledRSSBytes"],csvSHA256=hashlib.sha256(csv_file.read_bytes()).hexdigest(),
                  disconnectedLayer=missing,offlineRows=len(table))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("ui",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
    result=check(a.ui)
    with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
