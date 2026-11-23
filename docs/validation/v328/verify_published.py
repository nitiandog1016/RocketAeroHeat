"""仅凭 GitHub 小记录复核两份壁反力、远场载热、全程有界和收敛门槛。"""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent


def verify():
    d=json.loads((HERE.parents[2]/"cases/short_sphere/supportedCapacity65.json").read_bytes());scale=d["coefficientToW"]
    with ZipFile(HERE/"evidence/receipts.zip") as z:
        m=json.loads(z.read("manifest.json"))
        for path,digest in m["included"].items():assert hashlib.sha256(z.read(path)).hexdigest()==digest
        def rows(path):
            return list(csv.DictReader(io.StringIO(z.read(path).decode("ascii")),skipinitialspace=True))
        for name in ("early","late"):
            w=rows(name+"_wall.csv");f=rows(name+"_far.csv")
            assert len(w)==len({r["id"] for r in w})==5954 and len(f)==len({r["id"] for r in f})==1346
            vals=d[name]
            wall=math.fsum(-float(r["releasedRHS"])-float(r["fullCorrection"])+float(r["limitedCorrection"]) for r in w)*scale
            far=math.fsum(max(0.,float(r["scalarMass"]))*float(r["theta"]) for r in f)*scale
            alternate=math.fsum(max(0.,float(r["flowFDSMass"]))*float(r["theta"]) for r in f)*scale
            for a,b in ((wall,vals["wallW"]),(far,vals["scalarFarW"]),(alternate,vals["flowFluxComparisonW"])):
                assert math.isclose(a,b,abs_tol=1e-11,rel_tol=0)
            assert abs(far/wall-1)>.005 and abs(far/wall-1-vals["physicalGapFraction"])<1e-12
        for idx,folder in enumerate(("run01/steady01","run01/steady02","repeat03")):
            path="outputs/v328_validation/"+folder+"/"
            trace=rows(path+"stages.csv");history=rows(path+"history.csv")
            assert len(trace)==400 and len(history)==200
            for k,r in enumerate(trace):
                assert int(r["iteration"])==k//2 and r["stage"]==("actual" if k%2 else "low")
                assert -1e-12<=float(r["min"])<=float(r["max"])<=.4683950617283952+1e-12
            names=[s for s in history[0] if s.strip().startswith("rms[")]
            window=[max(float(r[n]) for r in history[-20:]) for n in names]
            assert window==d["stages"][idx]["window"] and len(window)==6
        assert max(window)<-8 and d["stages"][-1]["rawMax20"]<1e-8
    assert d["thermalStable"] and not d["balancePass"] and d["acceptedAreaM2"]==0
    print("600 updates and both endpoint heat budgets verified")


if __name__=="__main__":verify()
