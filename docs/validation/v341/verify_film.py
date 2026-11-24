"""用原生配置、热反力和发布面表独立复核膜温闭合、割线历史及同面物性差异。"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from zipfile import ZipFile

for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[var]="1"
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("short65_power_audit",HERE.parent/"v340/verify.py")
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
CASE=audit.CASE;PROOF=HERE/"evidence";STEM="calorimetryFilm65"


def relative(actual,expected):
    assert np.isfinite(actual) and abs(actual/expected-1)<1e-12,(actual,expected)


def configs(z,baseline,name,p,case):
    # 从试验来流、几何、湿度及当前候选温度重算，实验表温和局部 Nu 不进入输运。
    air=(case["airF"]-32)*5/9+273.15;wall=air+p["excessK"];film=(air+wall)/2
    k=.0241*(film/273.15)**1.5*(273.15+194)/(film+194)
    y=case["waterMassFraction"];pressure=case["pressureLbfFt2"]*4.4482216152605/.3048**2
    rho=pressure/(air*((1-y)*287.05+y*461.52));cp=(1-y)*1006+y*1860
    ud=case["diameterIn"]*.0254*case["speedFtS"]*.3048
    pe=rho*cp*ud/k;diff=1/pe;f=p["filmTransport"]
    c=json.loads(z.read(name+"/prepared/prepared.json"))["condition"]
    for a,b in ((f["wallK"],wall),(f["filmK"],film),(f["kWmK"],k),(f["pe"],pe),(f["diffusionND"],diff),
                (c["kFilm"],k),(c["pe"],pe),(c["rho"],rho),(c["cp"],cp),(c["mu"],rho*ud/case["reynoldsFree"])):
        relative(a,b)
    count=0
    for file in z.namelist():
        if not(file.startswith(name+"/run/") and file.endswith("/case.cfg")):continue
        path=file.rsplit("/",1)[0];label=path.split("/")[-1]
        template="steady01" if label.startswith("steady") else label
        src=baseline.read("high/run/"+template+"/case.cfg").decode("ascii").splitlines()
        raw=z.read(file);lines=raw.decode("ascii").splitlines()
        receipt=json.loads(z.read(path+"/debugger_process.json"))
        hashes=[h for n,h in receipt["inputSHA256"].items() if n.replace("\\","/").rsplit("/",1)[-1]=="case.cfg"]
        assert hashes==[hashlib.sha256(raw).hexdigest()] and len(src)==len(lines)
        changes={}
        for a,b in zip(src,lines):
            if a!=b:
                old_key,_,_=a.partition("=");key,_,val=b.partition("=")
                assert key==old_key and key not in changes;changes[key]=float(val)
        assert set(changes)=={"THERMAL_CONDUCTIVITY_CONSTANT","DIFFUSIVITY_CONSTANT","CONV_STARTITER"}
        relative(changes["THERMAL_CONDUCTIVITY_CONSTANT"],diff);relative(changes["DIFFUSIVITY_CONSTANT"],diff)
        assert changes["CONV_STARTITER"]==int(next(s.split("=")[1] for s in src if s.startswith("ITER=")))+2
        count+=1
    assert count==6+len(json.loads(z.read(name+"/run/converged.json"))["trend"])
    return dict(wallK=wall,filmK=film,kWmK=k,pe=pe,configurations=count)


def read_surface(stem):
    with ZipFile(CASE/(stem+".zip")) as z:
        summary=json.loads(z.read("summary.json"))
        wall=np.genfromtxt(z.open("wall.csv"),names=True,delimiter=",")
        stations=np.genfromtxt(z.open("stations.csv"),names=True,delimiter=",")
    return summary,wall,stations


def verify():
    proof=audit.verify(film=True,proof_dir=PROOF,stem=STEM)
    p=audit.read(CASE/(STEM+".json"));case=audit.read(CASE/"test65.json")
    order=p["order"];goal=p["targetSphereW"];final=p["balanceName"];rows={}
    assert order[:2]==["low","high"] and 3<=len(order)<=5 and order[-1]==final and p["balance"]==p["trials"][final]
    lo=p["trials"]["low"];hi=p["trials"]["high"];prior=p["source"]["frameSHA256"]
    assert (lo["excessK"],hi["excessK"])==(20,40)
    with ZipFile(PROOF/"native.zip") as z,ZipFile(HERE.parent/"v340/evidence/native.zip") as old:
        ledger=json.loads(z.read("acceptedTrials.json"))
        assert ledger["order"]==order and ledger["trials"]==p["trials"]
        assert hashlib.sha256(z.read("contracts/PLAN.md")).hexdigest()==p["source"]["planSHA256"]
        for i,name in enumerate(order):
            r=p["trials"][name];assert r["parentFrame"]==prior;prior=r["earlyFrame"]
            rows[name]=configs(z,old,name,r,case)
            if i<2:continue
            assert name==f"root{i-1:02d}"
            low=lo["earlyRegionW"][0];high=hi["earlyRegionW"][0];assert low<goal<high
            expect=lo["excessK"]+(goal-low)*(hi["excessK"]-lo["excessK"])/(high-low)
            assert abs(expect-r["excessK"])<1e-11
            success=all(abs(r[end+"RegionW"][0]/goal-1)<.0002 for end in ("early","late"))
            assert success==(name==final)
            if r["earlyRegionW"][0]<goal:lo=r
            else:hi=r
    s,w,stations=read_surface(STEM);old,v,st=read_surface("calorimetry65")
    assert s["frameSHA256"]==p["balance"]["earlyFrame"] and s["thermal"]["meshSHA256"]==old["thermal"]["meshSHA256"]
    relative(p["temperatureK"],rows[final]["wallK"])
    relative(s["condition"]["pe"],rows[final]["pe"]);relative(s["condition"]["kFilm"],rows[final]["kWmK"])
    for key in ("faceID","regionID","areaM2"):np.testing.assert_array_equal(w[key],v[key])
    fields={}
    for field in ("gaugePa","fluxInWm2","firstK","wallInputK"):
        delta=w[field]-v[field];i=np.argmax(abs(delta))
        fields[field]=dict(range=[min(w[field]),max(w[field])],deltaRMS=np.sqrt(np.average(delta**2,weights=w["areaM2"])),
                           maxAbsDelta=abs(delta[i]),signedAtMax=delta[i],referenceFaceID=int(w["faceID"][i]))
    base=audit.read(CASE/"calorimetry65.json");assert base["targetSphereW"]==goal
    dz=np.array(p["balance"]["earlyRegionW"])-base["balance"]["earlyRegionW"]
    assert abs(sum(dz)+(w["fluxInWm2"]-v["fluxInWm2"])@w["areaM2"])<1e-10
    np.testing.assert_array_equal(st["angleDeg"],stations["angleDeg"])
    proof.update(transport=rows,final=final,deltaWallK=p["temperatureK"]-base["temperatureK"],regionDeltaW=dz.tolist(),
                 fieldDifferences=fields,fixedNu=st["computedNu"].tolist(),fixedRMSE=old["stationRMSE"],fixedK=old["condition"]["kFilm"])
    return proof


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--publish",type=Path);parser.add_argument("--pack",type=Path)
    parser.add_argument("--record",action="store_true");args=parser.parse_args();PROOF.mkdir(exist_ok=True)
    if args.publish:audit.publish(args.publish.resolve(),STEM)
    if args.pack:audit.pack(args.pack.resolve(),film=True,proof_dir=PROOF)
    if not args.publish or args.pack:
        result=verify()
        if args.record:audit.store(PROOF/"comparison.json",result)
        print("Verified all candidate properties, native configurations, power brackets, surface fields and seven stations.")
