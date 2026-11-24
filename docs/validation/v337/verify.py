"""用归档节点、实际相邻时间层和独立多项式拟合复核银球测站的热流及误差。"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
from zipfile import ZipFile

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CASE=ROOT/"cases/short_sphere"
ARCHIVES={"steady5":"supportedSurface65.zip","steady10":"supportedFine65.zip",
          "steady128":"supportedAngularSurface65.zip","full":"supportedTime65Full.zip","half":"supportedTime65Half.zip"}


def wall_slope(h,t,indices):
    # 在首个采样点尺度上解范德蒙德系统，独立于生产实现的拉格朗日导数公式。
    pos=h[indices]/h[indices[1]]
    coef=np.linalg.solve(np.vander(pos,increasing=True),t[indices]-t[0])
    return coef[1]/h[indices[1]]


def check(path,volumes=False,reread=None):
    payload=json.loads(path.read_bytes())
    sets={s["label"]:s for s in payload["sets"]}
    assert payload["schema"]=="short65-station-gradient-v1" and set(sets)==set(ARCHIVES)
    test=json.loads((CASE/"test65.json").read_bytes())
    proof={}
    app=ROOT/"outputs/v336_validation/clean01"
    job=app/"outputs/short65_time_6dea1a68-073d-4129-867d-3772b804b131"
    if volumes:
        origin=json.loads((job/"construction.json").read_bytes())
        frames={"steady5":ROOT/"outputs/v329_validation/run01/early/initial_flow.csv",
                "steady10":ROOT/"outputs/v330_validation/run01/early/initial_flow.csv",
                "steady128":app/origin["steady"]/"early/initial_flow.csv",
                "full":job/"run/full/step1/restart_00001.csv","half":job/"run/half/step2/restart_00002.csv"}
    for label,run in sets.items():
        with ZipFile(CASE/ARCHIVES[label]) as z:
            s=json.loads(z.read("summary.json"));raw=z.read("stations.csv")
            stations=np.genfromtxt(io.BytesIO(raw),names=True,delimiter=",")
        assert hashlib.sha256(raw).hexdigest()==run["stationsSHA256"]
        assert run["frameSHA256"]==s["frameSHA256"] and run["caseSHA256"]==hashlib.sha256((CASE/"test65.json").read_bytes()).hexdigest()
        assert run["meshSHA256"]==(s["meshSHA256"] if "energy" in s else s["thermal"]["meshSHA256"])
        assert run["elapsedS"]==s.get("elapsedSeconds") and len(run["rows"])==len(stations)==7
        physical=s["condition"];diam=physical["diameterM"]
        if volumes:
            f=frames[label]
            assert hashlib.sha256(f.read_bytes()).hexdigest()==run["frameSHA256"]
            state=np.loadtxt(f,skiprows=1,delimiter=",",usecols=(0,1,2,3,9))
        errors=[];gap=[];wide=[];linear=[];experiment=[];react=[];heights=[]
        for i,(r,station) in enumerate(zip(run["rows"],stations)):
            assert r["status"]=="available" and r["angleDeg"]==test["angleDeg"][i]
            assert r["experimentNu"]==test["localNuSurface"][i] and r["faceID"]==station["originalFaceID"]
            np.testing.assert_allclose(r["wireGradientMagnitudeKPerM"],test["wireGradientFPerFt"][i]/.3048/1.8)
            assert r["reactionNu"]==station["computedNu"] and r["reactionInWm2"]==station["fluxInWm2"]
            pts=r["samples"];theta=np.deg2rad(r["angleDeg"]);n=np.array([np.cos(theta),np.sin(theta)])
            assert len(pts)==1+r["layerCount"] and pts[0]["ids"]==[station["leftNode"],station["rightNode"]]
            distances=[];temps=[]
            for p in pts:
                xyz=np.array(p["xyzD"]);K=np.array(p["nodeK"])
                if volumes:
                    nodes=state[p["ids"]]
                    np.testing.assert_array_equal(nodes[:,0],p["ids"])
                    np.testing.assert_array_equal(nodes[:,1:4],xyz)
                    np.testing.assert_allclose(physical["airK"]+physical["scalarUnitK"]*nodes[:,4],K,rtol=0,atol=1e-12)
                assert max(abs(xyz[:,1]))<1e-12
                xz=xyz[:,[0,2]]
                radius,w=np.linalg.solve(np.column_stack((n,xz[0]-xz[1])),xz[0])
                assert -1e-12<=w<=1+1e-12
                np.testing.assert_allclose(w,p["weight"],rtol=0,atol=1e-13)
                distances.append(radius*diam);temps.append((1-w)*K[0]+w*K[1])
            d=np.array(distances)-distances[0];t=np.array(temps)
            assert np.all(np.diff(d)>0) and d[-1]<r["joinPlaneClearanceM"]
            np.testing.assert_allclose(d,[p["distanceM"] for p in pts],rtol=0,atol=4e-18)
            np.testing.assert_allclose(t,[p["temperatureK"] for p in pts],rtol=0,atol=1e-12)
            np.testing.assert_allclose(t[0],station["wallK"],rtol=0,atol=1e-12)
            heats=physical["kFilm"]*np.array([wall_slope(d,t,k) for k in ([0,1],[0,1,2],[0,2,4])])
            errs=abs(heats-[r[k] for k in ("linearInWm2","gradientInWm2","wideInWm2")])
            assert max(errs)<1e-8
            kw=.0241*(t[0]/273.15)**1.5*(273.15+194)/(t[0]+194)
            nu=-heats*diam/(kw*(t[0]-physical["airK"]))
            np.testing.assert_allclose(nu,[r[k] for k in ("linearNu","gradientNu","wideNu")],atol=1e-10,rtol=0)
            gap.append(r["gradientNu"]-r["reactionNu"]);wide.append(nu[2]-nu[1]);linear.append(nu[0]-nu[1])
            experiment.append(r["gradientNu"]-r["experimentNu"]);react.append(r["reactionNu"]-r["experimentNu"])
            heights.append(1e6*d[1]);errors.append(max(errs))
        proof[label]=dict(maxReconstructionErrorWm2=max(errors),maxGradientReactionNu=max(abs(np.array(gap))),
                         gradientMinusReactionNu=gap,maxWideTemplateDeltaNu=max(abs(np.array(wide))),
                         maxLinearTemplateDeltaNu=max(abs(np.array(linear))),gradientExperimentErrorNu=experiment,
                         reactionRMSE=float(np.mean(np.square(react))**.5),gradientRMSE=float(np.mean(np.square(experiment))**.5),
                         firstDistanceRangeUm=[min(heights),max(heights)])
    if reread:
        sidecar=reread/"export/stationGradients.json"
        ack=json.loads((reread/"completed.json").read_bytes())
        assert hashlib.sha256(sidecar.read_bytes()).hexdigest()==ack["stationGradientsSHA256"]
        for branch,obs in zip(("full","half"),json.loads(sidecar.read_bytes())["sets"]):
            reference=dict(sets[branch]);reference.pop("label")
            assert obs==reference
    return proof


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--volumes",action="store_true");p.add_argument("--reread",type=Path)
    p.add_argument("--output",type=Path);a=p.parse_args()
    result=check(CASE/"stationGradients65.json",a.volumes,a.reread)
    if a.output:
        with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2);f.write("\n")
    print("Verified temperatures, ray intersections, signed heat fluxes and both Nusselt methods.")
