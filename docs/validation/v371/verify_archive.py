"""从原 STL、两份壁节点表、内迭代日志和逐层原生项独立核算外侧网格对照。"""
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np

DOC=Path(__file__).resolve().parent
DATA=DOC/"evidence"
CASE=DOC.parents[2]/"cases/simple_launcher"


def read(p):
    return json.loads(p.read_bytes())


def face_values(tri,table,columns):
    x,r=table["xM"],table["rM"]
    nodes=np.column_stack([table[k] for k in columns])
    base=abs(x)<1e-12
    start=np.flatnonzero(base)[-1]
    flat=(tri[:,:,0]==0).all(axis=1)
    idx=start+np.abs(tri[~flat,:,0,None]-x[start:]).argmin(axis=2)
    np.testing.assert_allclose(x[idx],tri[~flat,:,0],rtol=0,atol=1e-12)
    a=np.empty((len(tri),len(columns)))
    a[~flat]=nodes[idx].mean(axis=1)
    lo,hi=r[:start]/r[start],r[1:start+1]/r[start]
    a[flat]=((hi-lo)*(2*lo+hi)/3)@nodes[:start]+((hi-lo)*(2*hi+lo)/3)@nodes[1:start+1]
    return a


def verify():
    result=read(DATA/"comparison.json");step=read(DATA/"step.json");grid=read(DATA/"grid.json")
    geom=read(CASE/"geometry.json")
    base=DOC.parent/"v369/evidence"
    other=read(base/"surface.json")
    coarse_native=DOC.parent/"v370/evidence/balance.json"
    fine_native=DATA/"native/balance.json"
    assert result["frameSHA256"]==step["frameSHA256"]==read(fine_native)["frameSHA256"]
    assert other["frameSHA256"]==result["coarseFrameSHA256"]==read(coarse_native)["frameSHA256"]
    assert result["nativeBalance"]["fine"]["excessW"]==read(fine_native)["nativeExcessImbalanceW"]
    assert result["elapsedFromNumericalSeedUs"]==step["elapsedUs"]==6.25 and not result["physicalQualification"]
    assert result["wholeBodyPowerW"] is None and result["thermalSpatialComparisonPassed"]
    for p,k in ((fine_native,"fineNativeSHA256"),(coarse_native,"coarseNativeSHA256")):
        assert hashlib.sha256(p.read_bytes()).hexdigest()==result[k]
    knots=np.array(grid["originalLayerCoordinates"])
    np.testing.assert_array_equal(knots,np.sort(np.r_[np.arange(89),86.5,87.5]))
    assert grid["fineNodes"]==9919 and grid["cells"]==9720 and grid["fineSeedSHA256"]==step["seedSHA256"]
    err=abs(np.subtract(grid["newIntegrals"],grid["oldIntegrals"]))/grid["absoluteScales"]
    np.testing.assert_allclose(err,grid["relativeChanges"])
    assert max(err)<1e-12
    with ZipFile(DATA/"history.zip") as z:
        buf=z.read("history_00001.csv")
    assert hashlib.sha256(buf).hexdigest()==step["historySHA256"]
    stream=StringIO(buf.decode("utf-8-sig"));labels=[k.strip() for k in next(csv.reader(stream,skipinitialspace=True))]
    h=np.loadtxt(stream,delimiter=",")
    np.testing.assert_array_equal(h[:,labels.index("Inner_Iter")],np.arange(step["iterations"]))
    assert np.all(h[:,labels.index("Cur_Time")]==6.25e-6) and np.all(h[:,labels.index("Time_Step")]==6.25e-6)
    for k,v in step["last20Max"].items():assert max(h[-20:,labels.index(k)])==v and v < -8
    for key,p in (("coarse",coarse_native),("fine",fine_native)):
        b=read(p);rows=np.genfromtxt(p.parent/"layers.csv",names=True,delimiter=",")
        assert len(rows)==(89 if key=="coarse" else 91)
        phys=np.zeros(2)
        for k in ("farConvective","farViscous","wall","symmetry","storage"):
            phys+=np.array([rows[k+"_kg_s"].sum(),rows[k+"_W"].sum()])
        expected=result["nativeBalance"][key];href=b["referenceEnthalpyJkg"]
        np.testing.assert_allclose(phys,[expected["massKgS"],expected["energyW"]],rtol=1e-11,atol=1e-9)
        np.testing.assert_allclose(phys[1]-href*phys[0],expected["excessW"],atol=1e-5,rtol=0)
        # 粗网原坐标 86 到 87 线性升权，87 及以上权重为 1，即相同父网的末两层基函数。
        weight=np.clip((np.arange(89) if key=="coarse" else knots)-86,0,1)
        part=np.array([weight@rows["geometryMinusUniform_kg_s"],weight@rows["geometryMinusUniform_W"]])
        outer=result["originalOuterBand"][key]
        np.testing.assert_allclose([*part,part[1]-href*part[0]],[outer["massKgS"],outer["energyW"],outer["excessW"]],atol=1e-6,rtol=0)
        assert max(b["checksNativeND"].values())<1e-12
    stl=(CASE/geom["file"]).read_bytes()
    assert hashlib.sha256(stl).hexdigest()==result["surfaceSHA256"]
    dtype=np.dtype([("n","<f4",3),("v","<f4",(3,3)),("attr","<u2")])
    tri=np.frombuffer(stl,dtype=dtype,offset=84)["v"].astype(float)
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    groups_file=DOC.parent/"v365/research/surface.json"
    assert hashlib.sha256(groups_file.read_bytes()).hexdigest()==result["groupSourceSHA256"]
    ids=np.array(read(groups_file)["faceGroup"])-1
    coarse,fine=[face_values(tri,np.genfromtxt(p/"wall.csv",names=True,delimiter=","),result["columns"]) for p in (base,DATA)]
    np.testing.assert_allclose(coarse,np.array(other["groupValues"])[ids],atol=1e-9,rtol=1e-12)
    np.testing.assert_allclose(fine,np.array(result["fineGroupValues"])[ids],atol=1e-9,rtol=1e-12)
    ranges={r["name"]:slice(r["firstFace"]-1,r["lastFace"]) for r in geom["regions"]}
    ranges.update(all=slice(None),contractionEnd=slice(2016,2112),noseEnd=slice(6816,6912))
    for k,s in ranges.items():
        r=result["regions"][k];w=area[s];p,q=fine[s,:2].T
        dp,dq=(fine[s,:2]-coarse[s,:2]).T
        assert len(q)==r["faces"]
        np.testing.assert_allclose([w@q,w@coarse[s,1]],[r["researchIntoWallW"],r["coarseIntoWallW"]],atol=1e-9,rtol=0)
        qrms=np.sqrt((w@(dq*dq))/(w@(q*q)))
        prms=np.sqrt(w@(dp*dp)/sum(w))/result["dynamicPressurePa"]
        np.testing.assert_allclose([qrms,prms,abs(w@dq/(w@q))],[r["heatRmsFraction"],r["pressureRmsOverDynamic"],r["powerDifferenceFraction"]],atol=1e-12)
        np.testing.assert_allclose([abs(dp).max(),abs(dq).max()],r["maxAbsDifference"][:2],atol=1e-8)
        np.testing.assert_allclose(fine[r["peakOutwardFirstFace"]-1,1],min(q),atol=1e-8)
    print("Two matched physical frames, 9919 refined nodes, 9216 original faces and both native balances verified.")


if __name__=="__main__":
    verify()
