"""逐点核对导出图表的数据、完整 CSV、重排面号与真实工作台回执。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]


def inspect(folder,archived):
    with ZipFile(ROOT/"cases/short_sphere/stationAzimuth65.zip") as z:
        summary=json.loads(z.read("summary.json"));raw=list(csv.DictReader(io.StringIO(z.read("points.csv").decode())))
    view=json.loads((folder/"ui.json").read_bytes());proc=json.loads((folder/"window_process.json").read_bytes())
    assert not (proc["exitCode"] or proc["livePids"] or proc["stopReason"])
    assert view["fullStepIsolated"] and view["archiveMatches"] and view["damagedRejected"]
    assert view["gap"]["values"] is None and view["gap"]["gaps"]==2 and view["pole"]["rows"]==1
    frame=summary["frameSHA256"];base="Short65-azimuth-"+frame
    index=(folder/"workbench/index.html").read_text("utf8")
    assert f'src="{base}.html"' in index and f'href="{base}.csv"' in index
    page=(folder/"workbench"/(base+".html")).read_text("utf8")
    payload=json.loads(page.split("window.SHORT65_DATA=",1)[1].split(";</script>",1)[0])
    assert payload["frame"]==view["rear"]["frame"]==frame
    if archived:
        np.testing.assert_allclose(payload["gaps"],[[90,78.8700181154,101.1299818846],[90,258.8700181154,281.1299818846]],rtol=0,atol=1e-9)
        assert payload["rear"]==summary["rear"] or abs(payload["rear"]["lowerBound"]-summary["rear"]["lowerBound"])<1e-12
    stl=folder/("input.stl" if (folder/"input.stl").exists() else "reversed.stl")
    geom=hashlib.sha256(stl.read_bytes()).hexdigest()
    f=folder/"workbench"/(base+".csv")
    with f.open(encoding="utf-8-sig",newline="") as stream:output=list(csv.DictReader(stream))
    assert len(raw)==len(output)==view["rows"]==len(payload["rows"])
    numerical=("angleDeg","phiDeg","loDeg","hiDeg","xM","yM","zM","qInWm2","wallK","nu","nuLower","nuUpper")
    def number(x):return float(x) if x not in ("",None) else np.nan
    np.testing.assert_allclose([[number(r[k]) for k in numerical] for r in output],[[number(r[k]) for k in numerical] for r in raw],
                               rtol=3e-14,atol=1e-12,equal_nan=True)
    expected=[]
    for r,v in zip(raw,output):
        ref=int(r["faceID"]);target=ref if archived else 32481-ref
        assert v["frameSHA256"]==frame and v["stlSHA256"]==geom
        assert v["status"]==r["status"] and v["kind"]==r["kind"] and int(v["referenceFaceID"])==ref and int(v["fileFaceID"])==target
        expected.append([*[number(r[k]) for k in ("angleDeg","phiDeg","nu","qInWm2","wallK")],target,ref,payload["states"].index(r["status"])+1,int(r["faceCount"])])
    np.testing.assert_allclose(np.array(payload["rows"],dtype=float),expected,rtol=3e-14,atol=1e-12,equal_nan=True)
    assert view["selectedFileFace"]==view["otherAzimuth"]["faces"][0]
    return dict(rows=len(raw),fileFaceOrder="reference" if archived else "reversed",csvSHA256=hashlib.sha256(f.read_bytes()).hexdigest(),
                stlSHA256=geom,window=view,seconds=proc["elapsedSeconds"],peakRSSBytes=proc["peakSampledRSSBytes"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("local",type=Path);p.add_argument("archive",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
    r=dict(local=inspect(a.local,False),archive=inspect(a.archive,True))
    with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(r,f,indent=2);f.write("\n")
