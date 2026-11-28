"""只用发布的原生历史、壁节点、远界节点和原 STL 独立核对 v369。"""
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
from far_readback import calculate

DOC=Path(__file__).resolve().parent
DATA=DOC/"evidence"
CASE=DOC.parents[2]/"cases/simple_launcher"


def read(path):
    return json.loads(path.read_bytes())


def verify():
    runs=read(DATA/"runs.json")
    surface=read(DATA/"surface.json")
    assert set(runs)=={"full","half1"} and runs["full"]["sourceSHA256"]==runs["half1"]["sourceSHA256"]
    with ZipFile(DATA/"histories.zip") as z:
        assert set(z.namelist())=={"full.csv","half1.csv"}
        for key,rec in runs.items():
            blob=z.read(key+".csv")
            assert hashlib.sha256(blob).hexdigest()==rec["historySHA256"]
            text=StringIO(blob.decode("utf-8-sig"))
            labels=[s.strip() for s in next(csv.reader(text,skipinitialspace=True))]
            h=np.loadtxt(text,delimiter=",")
            np.testing.assert_array_equal(h[:,labels.index("Inner_Iter")],np.arange(rec["iterations"]))
            dt={"full":12.5,"half1":6.25}[key]
            assert (h[:,labels.index("Time_Step")]*1e6==dt).all() and (h[:,labels.index("Cur_Time")]*1e6==dt).all()
            vals=h[:,[labels.index(k) for k in rec["last20Max"]]]
            np.testing.assert_array_equal(vals[-20:].max(axis=0),list(rec["last20Max"].values()))
            windows=np.array([vals[j:j+20].max() for j in range(len(h)-19)])
            assert windows.min()==rec["bestWindowWorst"] and windows.argmin()==rec["bestWindowFirstIter"]
            assert (windows< -8).sum()==rec["passingWindows"] and bool(windows[-1]< -8)==rec["passed"]
            b=rec["fieldBalance"]
            np.testing.assert_allclose(b["massRateKgS"]+b["farOutMassKgS"],b["massImbalanceKgS"])
            np.testing.assert_allclose(b["storedTotalEnergyRateW"]+b["farOutTotalEnergyW"]+b["wallIntoSolidW"],b["totalEnergyImbalanceW"])
    assert runs["half1"]["passed"] and not runs["full"]["passed"] and not surface["sameTimePairAvailable"]
    assert surface["frameSHA256"]==runs["half1"]["frameSHA256"] and surface["wholeBodyPowerW"] is None
    assert not surface["physicalQualification"] and not surface["nativeFluxClosureProven"]
    geom=read(CASE/"geometry.json")
    stl=(CASE/geom["file"]).read_bytes()
    assert hashlib.sha256(stl).hexdigest()==surface["surfaceSHA256"]==geom["fileSHA256"]
    dtype=np.dtype([("n","<f4",3),("v","<f4",(3,3)),("a","<u2")])
    tri=np.frombuffer(stl,dtype=dtype,offset=84)["v"].astype(float)
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    wall=np.genfromtxt(DATA/"wall.csv",names=True,delimiter=",")
    nodes=np.column_stack([wall[k] for k in surface["columns"]])
    x,r=wall["xM"],wall["rM"]
    on_base=x==0
    start=np.flatnonzero(on_base)[-1]
    bot=(tri[:,:,0]==0).all(axis=1)
    corners=tri[~bot,:,0]
    idx=start+np.abs(corners[:,:,None]-x[start:]).argmin(axis=2)
    np.testing.assert_allclose(x[idx],corners,rtol=0,atol=1e-12)
    # 侧壁线性值取三顶点均值；每个平底扇形在节点半径区间按 r·dr 积分。
    vals=np.empty((len(area),nodes.shape[1]))
    vals[~bot]=nodes[idx].mean(axis=1)
    ra,rb=r[:start]/r[start],r[1:start+1]/r[start]
    vals[bot]=((rb-ra)*(2*ra+rb)/3)@nodes[:start]+((rb-ra)*(ra+2*rb)/3)@nodes[1:start+1]
    group_src=DOC.parent/"v365/research/surface.json"
    assert hashlib.sha256(group_src.read_bytes()).hexdigest()==surface["groupSourceSHA256"]
    group=np.array(read(group_src)["faceGroup"])-1
    np.testing.assert_allclose(vals,np.asarray(surface["groupValues"])[group],rtol=1e-11,atol=1e-9)
    np.testing.assert_allclose(area.sum(),surface["areaM2"],atol=1e-12,rtol=0)
    scopes={s["name"]:slice(s["firstFace"]-1,s["lastFace"]) for s in geom["regions"]}
    scopes.update(all=slice(None),contractionEnd=slice(2016,2112),noseEnd=slice(6816,6912))
    for key,s in scopes.items():
        summary=surface["regions"][key]
        assert len(vals[s])==summary["faces"]
        np.testing.assert_allclose(area[s]@vals[s,1],summary["researchIntoWallW"],atol=1e-9,rtol=0)
        np.testing.assert_allclose([vals[s,0].min(),vals[s,0].max()],summary["pressureRangePa"])
        np.testing.assert_allclose([vals[s,1].min(),vals[s,1].max()],summary["qIntoWallRangeWm2"])
    tw=wall["wallFluidK"]
    k=1.716e-5*(tw/273.15)**1.5*(273.15+110.4)/(tw+110.4)*1.4*287.05287/.4/.72
    np.testing.assert_allclose(k*(wall["firstFluidK"]-surface["wallGivenK"])/wall["normalGapM"],
                               wall["qIntoWallWm2"],rtol=2e-8,atol=2e-5)
    far=read(DATA/"far.json")
    a=np.loadtxt(DATA/"far_nodes.csv",delimiter=",",skiprows=1)
    result=calculate(a[:,:2],a[:,2:],far["condition"])
    assert far["frameSHA256"]==surface["frameSHA256"] and not result["nativeNumericalFluxMeasured"]
    for name in ("fluidGauss","fluidNodal","characteristicNodal","uniformNodal"):
        for label,value in result[name].items():
            np.testing.assert_allclose(value,far[name][label],atol=1e-5 if label=="netAdvectiveW" else 1e-9,rtol=1e-12)
    print(f'{len(tri)} faces, {len(wall)} wall nodes, {len(a)} far nodes, 2 original histories; qualified final field: false')


if __name__=="__main__":
    verify()
