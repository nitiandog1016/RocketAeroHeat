"""从公开压缩表独立核对同面差异、热峰集合及运行回执，可选核对导出的三维页。"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];DATA=ROOT/"cases/short_sphere"


def close(actual,expected):
    assert math.isclose(actual,expected,rel_tol=2e-12,abs_tol=1e-10),(actual,expected)


def records(z,name):
    return list(csv.DictReader(io.StringIO(z.read(name).decode("ascii")),skipinitialspace=True))


def read_case(name,key):
    path=DATA/name;assert hashlib.sha256(path.read_bytes()).hexdigest()==key
    with ZipFile(path) as z:
        s=json.loads(z.read("summary.json"))
        assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in s["files"].items())
        return s,records(z,"wall.csv"),records(z,"stations.csv")


def check_quantiles(values,weights,expected):
    ordered=sorted(zip(values,weights));total=math.fsum(weights)
    close(math.fsum(v*w for v,w in ordered)/total,expected["areaMean"])
    close(ordered[-1][0],expected["maximum"])
    for key,p in (("p50",.5),("p95",.95),("p99",.99)):
        summed=0
        for value,weight in ordered:
            summed+=weight
            if summed>=p*total:break
        close(value,expected[key])


def verify(ui=None):
    pair=json.loads((DATA/"supportedNearWall65.json").read_bytes())
    c,w,cs=read_case("supportedSurface65.zip",pair["coarseArchiveSHA256"])
    f,v,fs=read_case("supportedFine65.zip",pair["fineArchiveSHA256"])
    assert f["thermal"]["passed"] and f["frameSHA256"]==pair["fineFrame"]!=pair["coarseFrame"]==c["frameSHA256"]
    area=[];region=[];q=[];qf=[];pg=[];pf=[]
    for idx,(ra,rb) in enumerate(zip(w,v),1):
        assert int(ra["faceID"])==idx==int(rb["faceID"])
        assert all(float(ra[k])==float(rb[k]) for k in ("areaM2","wallInputK","regionID"))
        area.append(float(ra["areaM2"]));region.append(int(ra["regionID"]))
        q.append(float(ra["fluxInWm2"]));qf.append(float(rb["fluxInWm2"]))
        pg.append(float(ra["gaugePa"]));pf.append(float(rb["gaugePa"]))
    assert len(w)==len(v)==len(area)==pair["faces"]==32480
    close(math.fsum(area),pair["areaM2"])
    dq=[b-a for a,b in zip(q,qf)];dp=[b-a for a,b in zip(pg,pf)]
    for subset,entry in [(list(range(len(q))),pair)]+[([i for i,k in enumerate(region) if k==r["id"]],r) for r in pair["regions"]]:
        weights=[area[i] for i in subset]
        close(-math.fsum(q[i]*area[i] for i in subset),entry["coarseW"])
        close(-math.fsum(qf[i]*area[i] for i in subset),entry["fineW"])
        for differences,key in ((dq,"heatAbsDifference"),(dp,"pressureAbsDifference")):
            check_quantiles([abs(differences[i]) for i in subset],weights,entry[key])
    near=[]
    for values,key in ((q,"coarsePeak"),(qf,"finePeak")):
        minQ=min(values);ids={i for i,vv in enumerate(values) if vv<=.99*minQ};near.append(ids)
        assert sorted(i+1 for i in ids)==pair[key]["nearFaces"] and values[pair[key]["faceID"]-1]==minQ
        close(math.fsum(area[i] for i in ids),pair[key]["nearAreaM2"])
    close(math.fsum(area[i] for i in near[0]&near[1]),pair["sharedNearAreaM2"])
    close(math.fsum(area[i] for i in near[0]|near[1]),pair["unionNearAreaM2"])
    for i,(s,t) in enumerate(zip(cs,fs)):
        assert s["originalFaceID"]==t["originalFaceID"] and s["reportedNu"]==t["reportedNu"]
        close(float(s["computedNu"]),pair["coarseNu"][i]);close(float(t["computedNu"]),pair["fineNu"][i])
    if ui:
        page=json.loads((ui/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
        fields={f["key"]:f for f in page["fields"]};ranges=page["comparisonColorRanges"]
        for key,values,unit in (("heatFlux",qf,1000),("supportCoarseHeat",q,1000),("pressureGauge",pf,1),
                                ("supportCoarseGauge",pg,1),("supportGaugeChange",dp,1),("supportHeatChange",dq,1)):
            vals=fields[key]["values"];assert len(vals)==len(area)
            for a,b in zip(vals,values):close(a,b/unit)
        assert ranges["heatFlux"]==ranges["supportCoarseHeat"]==pair["sharedFluxRangeKWm2"]
        assert ranges["pressureGauge"]==ranges["supportCoarseGauge"]==pair["sharedGaugeRangePa"]
    with ZipFile(HERE/"evidence/receipts.zip") as z:
        m=json.loads(z.read("manifest.json"))
        for name,h in m["included"].items():assert hashlib.sha256(z.read(name)).hexdigest()==h
        grid=json.loads(z.read("grid.json"));assert grid["preparedSHA256"]==json.loads(z.read("build01/build.json"))["preparedSHA256"]
        assert grid["coarseNodes"]+grid["newMidpoints"]==grid["fineNodes"] and grid["meshSHA256"]==pair["fineMesh"]
        for k in ("early","late"):
            rows=records(z,k+"_wall.csv")
            watts=f["thermal"]["scale"]*math.fsum(-float(r["release"])-float(r["full"])+float(r["limited"]) for r in rows)
            far=records(z,k+"_far.csv")
            assert len(rows)==5954 and len(far)==1346 and all(r["flow"]==r["scalar"] for r in far)
            air=f["thermal"]["scale"]*math.fsum(max(0.,float(r["flow"]))*float(r["old"]) for r in far)
            close(watts,f["thermal"][k]["wallW"]);close(air,f["thermal"][k]["scalarFarW"]);assert abs(air/watts-1)<.005
    print("Verified two meshes, 32480 co-located face pairs, five regions, seven stations and both endpoint budgets.")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--ui",type=Path);verify(p.parse_args().ui)
