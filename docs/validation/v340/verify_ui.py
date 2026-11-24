"""逐行核对离线页、原面 CSV 和截面，分别验证球面求解温度与补偿支管给定温度。"""
import argparse
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CASE=ROOT/"cases/short_sphere"


def check(folder):
    ui=json.loads((folder/"ui.json").read_bytes());p=json.loads((CASE/"calorimetry65.json").read_bytes())
    with ZipFile(CASE/"calorimetry65.zip") as z:
        w=np.genfromtxt(z.open("wall.csv"),delimiter=",",names=True);stations=np.genfromtxt(z.open("stations.csv"),delimiter=",",names=True)
    rec=json.loads((folder/"window_process.json").read_bytes())
    assert rec["exitCode"]==0 and not(rec["livePids"] or rec["stopReason"])
    assert ui["frame"]==p["balance"]["earlyFrame"] and ui["faces"]==len(w) and ui["pixels"]>0 and ui["oldBoundaryRestored"]
    rev=ui["mode"]=="local";order=np.arange(len(w))[::-1] if rev else np.arange(len(w))
    dtype=[("n","<f4",(3,)),("vertex","<f4",(3,3)),("tag","<u2")];tri=np.fromfile(folder/"input.stl",dtype=dtype,offset=84)["vertex"].astype(float)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    w=w[order];sphere=w["regionID"]==1;q=w["fluxInWm2"]*(w["areaM2"]/area)
    assert ui["peakFile"]==np.argmin(q)+1 and ui["peakReference"]==w["faceID"][np.argmin(q)]
    with (folder/"wall.csv").open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    assert len(rows)==len(w)
    for i,r in enumerate(rows):
        assert int(r["reference_face_id"])==order[i]+1 and int(r["face_id"])==i+1 and r["joint_verified"]=="0"
        assert r["source_frame_sha256"]==ui["frame"] and r["pressure_method"] and r["heat_flux_method"]
        assert r["wall_temperature_role"]==("球体量热平衡温度" if sphere[i] else "支管或端盖随球温给定")
    def col(k):return np.array([float(r[k]) if r[k].strip() else np.nan for r in rows])
    np.testing.assert_allclose(col("pressure")*1000,w["absolutePa"],atol=1e-8,rtol=0)
    np.testing.assert_allclose(col("heat_flux")*1000,q,atol=1e-8,rtol=0)
    np.testing.assert_allclose(col("first_fluid_K"),w["firstK"],atol=1e-10,rtol=0)
    np.testing.assert_allclose(col("target_sphere_W"),p["targetSphereW"],atol=1e-12,rtol=0)
    solved=np.where(sphere,p["temperatureK"],np.nan);given=np.where(sphere,np.nan,w["wallInputK"])
    np.testing.assert_allclose(col("sphere_equilibrium_wall_K"),solved,atol=1e-10,rtol=0)
    np.testing.assert_allclose(col("wall_input_K"),given,atol=1e-10,rtol=0)
    page=json.loads((folder/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
    fields={f["key"]:f for f in page["fields"]}
    for key,x in (("pressure",w["absolutePa"]/1000),("heatFlux",q/1000),("wallSolved",solved),("wallInput",given),("nearWallCellTemperature",w["firstK"])):
        np.testing.assert_allclose(np.array(fields[key]["values"],float),x,atol=1e-10,rtol=0)
    assert "physicalStepSizeS" not in page["resultCondition"] and "elapsedPhysicalTimeS" not in page["resultCondition"]
    st=next((folder/"workbench").glob("Short65-stations-*.csv"))
    with st.open(encoding="utf-8-sig",newline="") as f:seven=list(csv.DictReader(f))
    assert len(seven)==7 and all(r["powerFrameSHA256"]==ui["frame"] for r in seven)
    np.testing.assert_allclose([float(r["computedNu"]) for r in seven],stations["computedNu"],atol=1e-10,rtol=0)
    with (folder/"section.csv").open(encoding="utf-8-sig",newline="") as f:cut=list(csv.DictReader(f))
    assert len(cut)==72
    for row in cut:
        assert row["field_key"]=="heatFlux" and abs(float(row["value"])-q[int(row["face_id"])-1]/1000)<1e-10
    return dict(window=ui,rows=len(w),solvedWallFaces=int(sum(sphere)),prescribedSupportFaces=int(sum(~sphere)),
                sectionRows=len(cut),sphereIntegralW=-q[sphere]@area[sphere],seconds=rec["elapsedSeconds"],peakBytes=rec["peakSampledRSSBytes"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("archive",type=Path);p.add_argument("local",type=Path);a=p.parse_args()
    result=dict(archive=check(a.archive),local=check(a.local))
    with (HERE/"evidence/workbench.json").open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
