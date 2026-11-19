"""从两个保存的三维终场在相同 X、方位复读物理量，并独立审核新 STL 的所有显示面。"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens
from cm011Tunnel import CP,D,DT,P0,RHO,T0,U,REGIONS,meridian,write_json
from runEmbeddedSU2 import sha


def read(path):
    return json.loads(path.read_bytes())


def flow_data(gui):
    run=Path((gui/"run_folder.txt").read_text("utf-8-sig").strip())
    flow=run/"flow"
    receipt=read(flow/"engine_process.json")
    assert receipt["status"]=="terminal" and receipt["exitCode"]==0 and not receipt["timedOut"]
    assert all(receipt["inputsStable"].values())
    for n in ("restart.csv","surface.vtk","history.csv"):
        assert sha(flow/n)==receipt["outputSHA256"][n]
    assert read(run/"audit/summary.json")["acceptedNumericalWall"]
    coords,field=vtk_tokens(flow/"surface.vtk")
    # 由体网原始坐标和全局节点号直接取值，独立于工作台的原面映射。
    grid=np.load(flow/"nodes.npz");ids=grid["ids"][:,:49];xyz=grid["xyz"]
    dist,ids_s=cKDTree(coords).query(xyz[ids])
    assert max(dist.ravel())<2e-6 and len(np.unique(ids_s))==len(coords)
    line=xyz[ids[6]][:,:2]*D
    assert np.all(np.diff(line[:,0])>0)
    restart=np.loadtxt(flow/"restart.csv",delimiter=",",skiprows=1)
    assert np.array_equal(restart[:,0],np.arange(len(xyz)))
    samples=dict(pressureGauge=field["Pressure"][:,0][ids_s]*RHO*U**2,
                 heatFlux=field["Heat_Flux"][:,0][ids_s]*RHO*U*CP*DT,
                 wallFluidTemperature=T0+(field["Temperature"][:,0][ids_s]-1)*DT,
                 nearWallCellTemperature=T0+(restart[grid["ids"][:,49:98],8]+297.15)*DT)
    return run,line,samples


def binary_stl(path):
    data=path.read_bytes()
    dt=np.dtype([("normal","<f4",3),("coords","<f4",(3,3)),("attr","<u2")])
    tris=np.frombuffer(data[84:],dtype=dt)["coords"].astype(float)
    assert len(tris)==int.from_bytes(data[80:84],"little")
    return tris


def native_curves(line,samples,x):
    return {k:np.array([np.interp(x,line[:,0],values) for values in mat])
            for k,mat in samples.items()}


def independent_mapping(vertices,faces,line,samples):
    """每个导入顶点与全部 48 个求解线段求二维射线交点，独立复核连续插值站号。"""
    nominal=meridian(1)[0]*D
    origin=np.array([.5*(nominal[0,0]+nominal[-1,0]),0.])
    rel=line-origin;seg=np.diff(rel,axis=0)
    xr=np.column_stack((vertices[:,0],np.hypot(vertices[:,1],vertices[:,2])))-origin
    ray=xr/np.linalg.norm(xr,axis=1)[:,None]
    num=ray[:,0,None]*rel[None,:-1,1]-ray[:,1,None]*rel[None,:-1,0]
    denom=ray[:,1,None]*seg[None,:,0]-ray[:,0,None]*seg[None,:,1]
    with np.errstate(divide="ignore",invalid="ignore"):
        t=num/denom
    valid=(t>=-1e-11)&(t<=1+1e-11)
    # 相反方向的共线延长线不能提供第二个交点。
    hits=rel[None,:-1,:]+np.where(valid,t,0)[:,:,None]*seg
    valid &= np.einsum("ik,ijk->ij",ray,hits)>0
    assert np.all(valid.any(axis=1))
    left=valid.argmax(axis=1);u=left+np.clip(t[np.arange(len(left)),left],0,1)
    angle=np.arctan2(abs(vertices[:,1]),-vertices[:,2])
    pos=np.column_stack((u,angle))
    return {k:RegularGridInterpolator((np.arange(49),np.linspace(0,np.pi,13)),v.T)(pos)[faces].mean(axis=1)
            for k,v in samples.items()}


def section_values(directory, tag, stl, x, output, angles, modules=None):
    """用工作台的截面模块导出，并由实际 Y-Z 端点求物理方位对应的面常值。"""
    program="""
const fs=require('fs'),path=require('path');
const [dir,x,out,tag,modules]=process.argv.slice(1);
global.window={}; require(path.join(dir,'case-data.js'));
const m=window.ROCKET_CASE, api=require(path.join(modules,'section.js'));
const viewer=require(path.join(modules,'viewer.js'));
const cut=api.geometry(m.vertices,m.faces,Number(x));
if(cut.components.length!==1||!cut.components[0].closed||cut.components[0].invalidAngleSegments)
  throw Error('截面必须是完整闭合轮廓');
for(const key of ['pressureGauge','heatFlux']){
  const f=m.fields.find(field=>field.key===key);
  // 和页面载入时一样，使用单一适用范围的字段没有逐面索引。
  if(!f.scopeNames.length){f.scopeNames=null; f.scopeIndex=null;}
  const csv=api.csv(cut,f,'冻结工况',viewer.caseCaption(m,f,0),{model:m});
  fs.writeFileSync(path.join(out,tag+'_'+key+'.csv'),csv,{flag:'wx'});
}
"""
    page=directory/"workbench"
    subprocess.run(["node","-e",program,str(page.resolve()),str(x),str(output.resolve()),
                    tag,str((modules or page).resolve())],check=True)
    model=json.loads((page/"case-data.js").read_text("utf8").strip().split("=",1)[1][:-1])
    tri=binary_stl(stl)
    result={}
    for key in ("pressureGauge","heatFlux"):
        with (output/(tag+"_"+key+".csv")).open(encoding="utf-8-sig",newline="") as stream:
            records=list(csv.DictReader(stream))
        field=next(z for z in model["fields"] if z["key"]==key)
        starts,ends,values=[],[],[]
        for row in records:
            i=int(row["面元编号"])-1
            value=float(row[field["label"]+" / "+field["unit"]])
            assert value==field["values"][i] and row["算法"]==field["methodNames"][field["methodIndex"][i]]
            assert row["适用范围"]==field["scope"] and "+10°" in row["工况"]
            assert min(tri[i,:,0])-1e-12<=x<=max(tri[i,:,0])+1e-12
            starts.append([float(row[k]) for k in ("A_Y / m","A_Z / m")])
            ends.append([float(row[k]) for k in ("B_Y / m","B_Z / m")])
            values.append(value*(1000 if key=="heatFlux" else 1))
        starts=np.asarray(starts)
        edges=np.asarray(ends)-starts
        picked=[]
        for theta in np.deg2rad(angles):
            direction=np.array([np.cos(theta),np.sin(theta)])
            with np.errstate(divide="ignore",invalid="ignore"):
                t=-np.cross(starts,direction)/np.cross(edges,direction)
            hits=(t>=0)&(t<=1)&((starts+np.where(np.isfinite(t),t,0)[:,None]*edges)@direction>0)
            assert hits.sum()==1
            picked.append(values[int(np.flatnonzero(hits)[0])])
        result[key]=dict(segments=len(records),valuesSI=picked)
    return result


def check(gui,out):
    old_gui=ROOT/"outputs/v275_validation/gui01"
    old_run,old_line,old=flow_data(old_gui)
    run,line,new=flow_data(gui)
    stats=read(run/"audit/summary.json")
    page=gui/"workbench/case-data.js"
    js=page.read_text("utf8").strip();assert js.startswith("window.ROCKET_CASE=") and js.endswith(";")
    model=json.loads(js[len("window.ROCKET_CASE="):-1])
    v,f=np.array(model["vertices"]),np.array(model["faces"])
    tri=binary_stl(ROOT/"cases/CM011-alternate-SI.stl")
    assert len(f)==len(tri) and np.all(np.min([abs(tri-np.roll(v[f],i,axis=1)).max(axis=(1,2))
                                                for i in range(3)],axis=0)==0)
    weights=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    np.testing.assert_allclose(model["faceArea"],weights,rtol=2e-13,atol=1e-18)
    mapped=independent_mapping(v,f,line,new)
    fields={z["key"]:z for z in model["fields"]}
    max_diff={}
    for k,expected in mapped.items():
        z=fields[k]
        shown=np.asarray(z["values"])*(1000 if k=="heatFlux" else 1)
        assert np.all(np.array(z["status"])==1) and all(z["methodNames"][i] for i in z["methodIndex"])
        np.testing.assert_allclose(shown,expected,atol=1e-10,rtol=1e-13)
        max_diff[k]=float(max(abs(shown-expected)))
    np.testing.assert_allclose(np.asarray(fields["pressure"]["values"])*1000,P0+mapped["pressureGauge"],atol=1e-8,rtol=0)
    assert set(fields["wallInput"]["status"])=={2}
    assert all(set(fields[k]["status"])=={0} for k in ("wallSolved","edgeTemperature","recoveryTemperature"))
    power=weights@mapped["heatFlux"]; old_stats=read(old_run/"audit/summary.json")
    assert abs(power-model["wallCoverage"]["wholeSurfaceHeatPowerW"])<1e-15
    power_diff=float(power/old_stats["triangularPowerW"]-1)
    assert abs(power_diff)<.01
    comparison=[]
    for i,k in enumerate((3,14,29,41)):
        x=old_line[k,0];a=native_curves(old_line,old,x);b=native_curves(line,new,x)
        dp=max(abs(a["pressureGauge"]-b["pressureGauge"]));dq=max(abs(a["heatFlux"]-b["heatFlux"]))
        assert dp<2e-4 and dq<.15
        comparison.append(dict(region=REGIONS[i],xM=x,referenceRadiusM=old_line[k,1],
             newRadiusM=float(np.interp(x,line[:,0],line[:,1])),maxPressureDifferencePa=dp,maxHeatDifferenceWm2=dq,
             referenceGaugePa=a["pressureGauge"].tolist(),newGaugePa=b["pressureGauge"].tolist(),
             referenceHeatWm2=a["heatFlux"].tolist(),newHeatWm2=b["heatFlux"].tolist()))
    out.mkdir(parents=True,exist_ok=False)
    x=comparison[1]["xM"]
    angles=np.arange(7.5,360.,15.)
    curves={tag:section_values(directory,tag,ROOT/"cases"/name,x,out,angles)
            for tag,directory,name in (("original",old_gui,"CM011-mm-scale-SI.stl"),
                                       ("alternate",gui,"CM011-alternate-SI.stl"))}
    curve_diff={k:float(max(abs(np.array(curves["original"][k]["valuesSI"])-curves["alternate"][k]["valuesSI"])))
                for k in ("pressureGauge","heatFlux")}
    return dict(faces=len(f),referenceFaces=old_stats["originalFaces"],areaM2=weights.sum(),referenceAreaM2=old_stats["areaM2"],
                powerW=power,referencePowerW=old_stats["triangularPowerW"],relativePowerDifference=power_diff,
                commonPositions=comparison,sectionXM=x,sectionAnglesFromPositiveYDeg=angles.tolist(),
                sectionFaceCurves=curves,maxSectionFaceDifferenceSI=curve_diff,maximumIndependentMapDifferenceSI=max_diff,
                originalFaceOrderVerified=True,caseSHA256=sha(page),sourceRestartSHA256=sha(run/"flow/restart.csv"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gui",type=Path);parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    write_json(args.new_output/"comparison.json",check(args.gui,args.new_output))
