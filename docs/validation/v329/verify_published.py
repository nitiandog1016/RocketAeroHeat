"""仅用发布包和 Python 标准库复核原面清单、五区功率、远场传递及七站对照。"""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
DATA=HERE.parents[2]/"cases/short_sphere"


def rows(z,name):
    return csv.DictReader(io.StringIO(z.read(name).decode("utf8")),skipinitialspace=True)


def verify():
    with ZipFile(DATA/"supported65.zip") as g:geometry=json.loads(g.read("geometry.json"))
    with ZipFile(DATA/"supportedSurface65.zip") as z:
        s=json.loads(z.read("summary.json"));pieces=[[] for _ in range(5)];areas=[[] for _ in range(5)]
        for name,h in s["files"].items():assert hashlib.sha256(z.read(name)).hexdigest()==h
        for i,r in enumerate(rows(z,"wall.csv")):
            k=int(r["regionID"])-1;p=float(r["absolutePa"]);area=float(r["areaM2"]);q=float(r["fluxInWm2"])
            assert i+1==int(r["faceID"]) and k+1==geometry["regionID"][i] and area==geometry["areaM2"][i]
            assert p>0 and abs(p-s["condition"]["airPa"]-float(r["gaugePa"]))<1e-8
            assert math.isfinite(q) and float(r["firstK"])>0 and float(r["wallInputK"])>0
            areas[k].append(area);pieces[k].append(-q*area)
        regionW=[math.fsum(v) for v in pieces]
        assert sum(map(len,areas))==len(geometry["F"])==s["faces"]==32480
        for k,r in enumerate(s["regions"]):
            assert len(areas[k])==r["faces"] and abs(math.fsum(areas[k])-r["areaM2"])<1e-15 and abs(regionW[k]-r["outW"])<1e-11
        totals=s["thermal"]
        assert abs(sum(regionW)-totals["early"]["wallW"])<1e-10
        case=json.loads((DATA/"test65.json").read_bytes());ref=list(rows(z,"stations.csv"))
        assert [float(r["reportedNu"]) for r in ref]==case["localNuSurface"]
        assert abs(math.sqrt(sum((float(r["computedNu"])-float(r["reportedNu"]))**2 for r in ref)/len(ref))-s["stationRMSE"])<1e-12
    with ZipFile(HERE/"evidence/receipts.zip") as z:
        for name,h in json.loads(z.read("manifest.json"))["included"].items():assert hashlib.sha256(z.read(name)).hexdigest()==h
        for n in ("early","late"):
            far=list(rows(z,n+"_far.csv"))
            assert len({r["id"] for r in far})==1346 and all(r["flowMass"]==r["scalarMass"] for r in far)
            scale=totals["scale"];power=scale*math.fsum(max(0.,float(r["flowMass"]))*float(r["old"]) for r in far)
            wall=scale*math.fsum(-float(r["releaseRHS"])-float(r["fullCorrection"])+float(r["limitedCorrection"]) for r in rows(z,n+"_wall.csv"))
            assert abs(power-totals[n]["scalarFarW"])<1e-10 and abs(wall-totals[n]["wallW"])<1e-10
            assert abs(power/wall-1)<.005
        assert json.loads(z.read("repeat.json"))["gate"] and totals["thermalChange"]<.002
    print("Verified 32480 original faces, five regions, seven stations, both heat balances and the cold GUI run.")


if __name__=="__main__":verify()
