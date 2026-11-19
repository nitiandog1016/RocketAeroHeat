"""用二进制 STL 和独立的面心匹配检验编号、场量、功率及截面导出。"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Imported import instantiate, match_surface
from cm011Tunnel import write_json
from runEmbeddedSU2 import sha


def read(path):
    return json.loads(path.read_bytes())


def geometry_checks(stl):
    dt=np.dtype([("normal","<f4",3),("xyz","<f4",(3,3)),("attr","<u2")])
    raw=stl.read_bytes()
    triangles=np.frombuffer(raw[84:],dtype=dt)["xyz"].astype(float)
    assert len(triangles)==int.from_bytes(raw[80:84],"little")
    v,flat=np.unique(triangles.reshape(-1,3),axis=0,return_inverse=True)
    f=flat.reshape(-1,3)+1
    body,reg=match_surface(v,f)
    rng=np.random.RandomState(275)
    vp,fp=rng.permutation(len(v)),rng.permutation(len(f))
    mapping=np.argsort(vp)
    faces=np.roll(mapping[f[fp]-1]+1,1,axis=1)
    other,check=match_surface(v[vp],faces)
    assert np.array_equal(other["F"],faces) and np.array_equal(other["V"],v[vp])
    assert np.array_equal(np.array(reg["referenceFaceOneBased"])[fp],check["referenceFaceOneBased"])
    assert other["meridian"]==body["meridian"] and np.array_equal(np.array(body["regionCode"])[fp],other["regionCode"])
    damaged=v.copy();damaged[137,2]+=2e-8
    reversed_faces=f.copy();reversed_faces[21]=reversed_faces[21,::-1]
    dup=f.copy();dup[1]=dup[0]
    for verts,tri in ((damaged,f),(v*.001,f),(v,reversed_faces),(v,dup)):
        try:match_surface(verts,tri)
        except ValueError:pass
        else:raise AssertionError("错误几何未被拒绝")
    return triangles,dict(vertices=len(v),faces=len(f),shuffledFacesAndVertices=True,
                          changedCoordinateRejected=True,wrongScaleRejected=True,reversedFaceRejected=True,duplicateFaceRejected=True)


def verify(gui,out):
    stl=ROOT/"cases/CM011-mm-scale-SI.stl"
    tri,tests=geometry_checks(stl)
    run=Path((gui/"run_folder.txt").read_text("utf-8-sig").strip())
    prep=Path((gui/"mesh_folder.txt").read_text("utf-8-sig").strip())
    native=run/"flow"
    geom=read(native/"geometry.json")
    native_body=read(native/"body.json")
    receipt=read(native/"engine_process.json")
    assert receipt["status"]=="terminal" and receipt["exitCode"]==0 and all(receipt["inputsStable"].values())
    assert read(native/"resources.json")["livePids"]==[]
    with tempfile.TemporaryDirectory() as t:
        try:instantiate(prep,Path(t)/"bad","0"*64)
        except ValueError:pass
        else:raise AssertionError("错误体网指纹未被拒绝")
    page=gui/"workbench/case-data.js"
    s=page.read_text("utf8").strip();assert s.startswith("window.ROCKET_CASE=") and s.endswith(";")
    model=json.loads(s[len("window.ROCKET_CASE="):-1])
    v,f=np.array(model["vertices"]),np.array(model["faces"])
    triangles=v[f]
    # STL 文件的三角面顺序必须与工作台编号一致，允许每面循环改变首顶点。
    assert len(tri)==len(f) and np.all(np.min([np.max(abs(tri-np.roll(triangles,k,axis=1)),axis=(1,2))
                                               for k in (0,1,2)],axis=0)==0)
    assert np.array_equal(v,native_body["V"]) and np.array_equal(f+1,native_body["F"])
    prior=ROOT/"outputs/v274_validation/half01"
    ref=read(prior/"body.json")
    ref_tri=np.asarray(ref["V"])[np.asarray(ref["F"])-1]
    dist,idx=cKDTree(ref_tri.mean(axis=1)).query(triangles.mean(axis=1))
    assert max(dist)<5e-10 and len(np.unique(idx))==len(f)
    assert np.array_equal(idx+1,geom["importedSurface"]["referenceFaceOneBased"])
    original=np.genfromtxt(ROOT/"outputs/v274_validation/audit01/wall.csv",delimiter=",",names=True)[idx]
    fields={field["key"]:field for field in model["fields"]}
    diffs={}
    for key,col,scale,tol in (("pressureGauge","gaugePressurePa",1.,1e-5),("heatFlux","heatFluxWm2",1000.,.005),
                             ("wallFluidTemperature","wallFluidTemperatureK",1.,1e-4),("nearWallCellTemperature","nearWallCellTemperatureK",1.,1e-4)):
        field=fields[key]
        diff=float(max(abs(scale*np.array(field["values"])-original[col])))
        assert diff<tol and set(field["status"])=={1} and len(field["values"])==len(f)
        assert all(len(field["methodNames"][i]) for i in field["methodIndex"])
        diffs[key]=diff
    areas=.5*np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)
    power=areas@(np.asarray(fields["heatFlux"]["values"])*1000)
    ref_power=read(ROOT/"docs/validation/v274/evidence/steady.json")["triangularPowerW"]
    assert abs(power/ref_power-1)<.001 and abs(power-model["wallCoverage"]["wholeSurfaceHeatPowerW"])<1e-15
    assert model["resultCondition"]["alphaDeg"]==10 and model["researchIdentity"]["productionQualified"] is False
    assert model["researchIdentity"]["stlSHA256"]==sha(stl) and model["researchIdentity"]["preparedSHA256"]==sha(prep/"prepared.json")
    assert set(fields["wallInput"]["status"])=={2} and set(fields["wallSolved"]["status"])=={0}
    stats=read(run/"audit/summary.json");x=stats["coneRing"]["xM"]
    # 从最终导出所用的实际截面代码生成 CSV；独立验证每个交线段的原面来源。
    out.mkdir(parents=True,exist_ok=False)
    code="""
const fs=require('fs'),path=require('path');const [dir,x,out]=process.argv.slice(1);
global.window={};require(path.join(dir,'case-data.js'));const m=window.ROCKET_CASE;
const api=require(path.join(dir,'section.js')),viewer=require(path.join(dir,'viewer.js'));
const cut=api.geometry(m.vertices,m.faces,Number(x));
if(cut.components.length!==1||!cut.components[0].closed||cut.components[0].invalidAngleSegments)throw Error('section');
for(const k of ['pressureGauge','heatFlux']){const f=m.fields.find(f=>f.key===k);
fs.writeFileSync(path.join(out,k+'.csv'),api.csv(cut,f,'冻结工况',viewer.caseCaption(m,f,0),{model:m}),{flag:'wx'});}
"""
    subprocess.run(["node","-e",code,str(page.parent.resolve()),str(x),str(out.resolve())],check=True)
    counts={}
    for key in ("pressureGauge","heatFlux"):
        field=fields[key]
        with (out/(key+".csv")).open(encoding="utf-8-sig",newline="") as inp:
            rows=list(csv.DictReader(inp))
        for r in rows:
            face=int(r["面元编号"])-1
            assert float(r[field["label"]+" / "+field["unit"]])==field["values"][face]
            assert r["算法"]==field["methodNames"][field["methodIndex"][face]] and "+10°" in r["工况"]
            assert min(tri[face,:,0])-1e-12<=x<=max(tri[face,:,0])+1e-12
        assert len(rows)>0
        counts[key]=len(rows)
    return dict(tests,meshSHA256=sha(native/"mesh.su2"),bodySHA256=sha(native/"body.json"),stlSHA256=sha(stl),
                preparedSHA256=sha(prep/"prepared.json"),sectionSegments=counts,maxDifferenceInSI=diffs,
                fullAreaM2=areas.sum(),powerIntoWallW=power,referencePowerW=ref_power,
                sectionXM=x,maxCentroidDistanceM=max(dist),payloadSHA256=sha(page),originalSTLOrder=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gui",type=Path);parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    write_json(args.new_output/"verification.json",verify(args.gui,args.new_output))
