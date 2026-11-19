"""在相同原面和物理截面比较 CM-011 两种半周离散，另存局部超差与全壁功率。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"docs/validation/v276"))
from compare import binary_stl, native_curves, section_values
from auditSU2SSTState import vtk_tokens
from cm011Tunnel import CP, D, DT, P0, RHO, T0, U, write_json
from runEmbeddedSU2 import sha


def read(path):
    return json.loads(path.read_bytes())


def saved_result(gui, sectors, document):
    run=Path((gui/"run_folder.txt").read_text("utf-8-sig").strip())
    flow=run/"flow"
    receipt=read(flow/"engine_process.json")
    assert receipt["status"]=="terminal" and receipt["exitCode"]==0 and not receipt["timedOut"]
    assert not receipt["cancelled"] and all(receipt["inputsStable"].values())
    stats=read(run/"audit/summary.json")
    assert stats["normalTerminal"] and stats["acceptedNumericalWall"] and not stats["physicalAccuracyValidated"]
    for n in ("surface.vtk","restart.csv","history.csv"):
        assert sha(flow/n)==receipt["outputSHA256"][n]
    assert sha(flow/"restart.csv")==document["sourceRestartSHA256"] and document["belowHalfPercent"]
    assert sha(flow/"body.json")=="de44b9ac4b3a74ac530f7c339dda8d31c62fb4f872a9183616188f6d6a65d2b2"
    body=read(flow/"body.json")
    page=run/"workbench/case-data.js"
    model=json.loads(page.read_text("utf8").strip().split("=",1)[1][:-1])
    vertices=np.array(model["vertices"])
    faces=np.array(model["faces"])
    np.testing.assert_array_equal(faces,np.array(body["F"])-1)
    np.testing.assert_array_equal(vertices,body["V"])
    tri=binary_stl(ROOT/"cases/CM011-alternate-SI.stl")
    assert np.min([np.max(abs(tri-np.roll(vertices[faces],k,axis=1)),axis=(1,2)) for k in range(3)],axis=0).max()==0
    centers=tri.mean(axis=1)
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    np.testing.assert_allclose(area,model["faceArea"],rtol=2e-13,atol=1e-18)
    grid=np.load(flow/"nodes.npz")
    xyz=grid["xyz"]
    ids=grid["ids"][:,:49]
    wall_xyz=xyz[ids]
    phi=np.arctan2(wall_xyz[:,14,2],wall_xyz[:,14,1])
    np.testing.assert_allclose(phi,np.linspace(-np.pi/2,np.pi/2,sectors+1),rtol=0,atol=1e-14)
    line=np.stack((wall_xyz[0,:,0],np.hypot(wall_xyz[0,:,1],wall_xyz[0,:,2])),axis=1)*D
    coords,fields=vtk_tokens(flow/"surface.vtk")
    distance,order=cKDTree(coords).query(wall_xyz)
    assert distance.max()<2e-6 and len(np.unique(order))==len(coords)
    restart=np.loadtxt(flow/"restart.csv",delimiter=",",skiprows=1)
    np.testing.assert_array_equal(restart[:,0],np.arange(len(xyz)))
    native=dict(pressureGauge=fields["Pressure"][:,0][order]*RHO*U**2,
                heatFlux=fields["Heat_Flux"][:,0][order]*RHO*U*CP*DT,
                wallFluidTemperature=T0+(fields["Temperature"][:,0][order]-1)*DT,
                nearWallCellTemperature=T0+(restart[grid["ids"][:,49:98],8]+297.15)*DT)
    points=np.column_stack((np.arctan2(vertices[:,2],abs(vertices[:,1])),body["wallStation"]))
    by_key={field["key"]:field for field in model["fields"]}
    face_values={}
    for key,samples in native.items():
        values=RegularGridInterpolator((phi,np.arange(49)),samples)(points)[faces].mean(axis=1)
        field=by_key[key]
        shown=np.array(field["values"])*(1000 if key=="heatFlux" else 1)
        assert set(field["status"])=={1} and all(field["methodNames"][i] for i in field["methodIndex"])
        np.testing.assert_allclose(shown,values,rtol=1e-13,atol=1e-10)
        face_values[key]=values
    np.testing.assert_allclose(np.asarray(by_key["pressure"]["values"])*1000,P0+face_values["pressureGauge"],rtol=0,atol=1e-8)
    assert min(P0+face_values["pressureGauge"])>0 and set(by_key["wallInput"]["status"])=={2}
    assert all(set(by_key[key]["status"])=={0} for key in ("wallSolved","edgeTemperature","recoveryTemperature"))
    assert abs(area@face_values["heatFlux"]-model["wallCoverage"]["wholeSurfaceHeatPowerW"])<1e-15
    assert model["researchIdentity"]["productionQualified"] is False and model["resultCondition"]["alphaDeg"]==10
    return dict(run=run,native=native,faces=face_values,line=line,phi=phi,area=area,centers=centers,
                power=float(area@face_values["heatFlux"]),body=body,source=sha(flow/"restart.csv"),payload=sha(page))


def extrema(values, centers, regions):
    result={}
    for name,index in (("minimum",np.argmin(values)),("maximum",np.argmax(values))):
        result[name]=dict(faceOneBased=int(index)+1,valueSI=float(values[index]),
                          xyzM=centers[index].tolist(),region=regions[index])
    return result


def check(base,out):
    old=saved_result(ROOT/"outputs/v276_validation/gui01",12,read(ROOT/"docs/validation/v276/evidence/energy.json"))
    new=saved_result(base/"gui01",24,read(base/"energy.json"))
    for key in ("area","line","centers"):
        np.testing.assert_allclose(old[key],new[key],rtol=0,atol=1e-16)
    assert old["body"]==new["body"] and len(new["area"])==10080
    out.mkdir(exist_ok=False)
    area=new["area"]
    limits={"pressureGauge":2e-4,"heatFlux":.15}
    regions=[new["body"]["regionNames"][i] for i in new["body"]["regionCode"]]
    local={}
    for key,tolerance in limits.items():
        delta=new["faces"][key]-old["faces"][key]
        i=int(np.argmax(abs(delta)))
        local[key]=dict(maximumAbsoluteDifferenceSI=float(abs(delta[i])),areaWeightedRmsSI=float(np.sqrt(area@(delta**2)/area.sum())),
                        maxDifferenceFaceOneBased=i+1,maxDifferenceXYZM=new["centers"][i].tolist(),
                        facesAtOrAboveSectionTolerance=int(np.sum(abs(delta)>=tolerance)),
                        areaAtOrAboveSectionToleranceM2=float(area[abs(delta)>=tolerance].sum()),
                        coarseExtrema=extrema(old["faces"][key],old["centers"],regions),
                        fineExtrema=extrema(new["faces"][key],new["centers"],regions))
    bands=[]
    pressure_difference=new["faces"]["pressureGauge"]-old["faces"]["pressureGauge"]
    heat_difference=new["faces"]["heatFlux"]-old["faces"]["heatFlux"]
    pressure_flag=abs(pressure_difference)>=limits["pressureGauge"]
    heat_flag=abs(heat_difference)>=limits["heatFlux"]
    with (out/"sensitive_faces.csv").open("x",encoding="utf-8-sig",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["原面编号","区域","X / m","Y / m","Z / m","面积 / m²","12单元表压 / Pa",
                         "24单元表压 / Pa","表压差 / Pa","12单元入壁热流 / (W/m²)",
                         "24单元入壁热流 / (W/m²)","热流差 / (W/m²)"])
        for i in np.flatnonzero(pressure_flag|heat_flag):
            writer.writerow([i+1,regions[i],*new["centers"][i],area[i],old["faces"]["pressureGauge"][i],
                             new["faces"]["pressureGauge"][i],pressure_difference[i],old["faces"]["heatFlux"][i],
                             new["faces"]["heatFlux"][i],heat_difference[i]])
    for region in new["body"]["regionNames"]:
        select=np.asarray(regions)==region
        bands.append(dict(region=region,faces=int(select.sum()),areaM2=float(area[select].sum()),
                          pressureFlaggedFaces=int(np.sum(select&pressure_flag)),
                          heatFlaggedFaces=int(np.sum(select&heat_flag)),
                          power12W=float(area[select]@old["faces"]["heatFlux"][select]),
                          power24W=float(area[select]@new["faces"]["heatFlux"][select])))
    positions=read(ROOT/"docs/validation/v276/evidence/comparison.json")["commonPositions"]
    sections=[]
    for station in positions:
        x=station["xM"]
        coarse=native_curves(old["line"],old["native"],x)
        fine=native_curves(new["line"],new["native"],x)
        section=dict(region=station["region"],xM=x,coarseAzimuthDeg=np.rad2deg(old["phi"]).tolist(),
                     fineAzimuthDeg=np.rad2deg(new["phi"]).tolist())
        for key,tolerance in limits.items():
            a,b=coarse[key],fine[key]
            diff=float(np.max(abs(b[::2]-a)))
            section[key]=dict(coarse=a.tolist(),fine=b.tolist(),maxCommonDifferenceSI=diff,
                              withinTolerance=diff<tolerance,coarseWindwardMinusLeeward=float(a[0]-a[-1]),
                              fineWindwardMinusLeeward=float(b[0]-b[-1]))
        sections.append(section)
    x=positions[1]["xM"]
    angles=np.arange(7.5,360.,15.)
    curves={}
    for tag,result in (("coarse",old),("fine",new)):
        curves[tag]=section_values(result["run"],tag,ROOT/"cases/CM011-alternate-SI.stl",x,out,angles)
    face_differences={key:float(max(abs(np.asarray(curves["coarse"][key]["valuesSI"])-curves["fine"][key]["valuesSI"])))
                      for key in limits}
    relative=new["power"]/old["power"]-1
    result=dict(nativeSections=sections,regions=bands,faceComparison=local,originalAreaM2=float(area.sum()),
                power12W=old["power"],power24W=new["power"],relativePowerChange=relative,
                coneSectionXM=x,sectionAnglesFromPositiveYDeg=angles.tolist(),sectionFaceCurves=curves,
                sectionFaceDifferencesSI=face_differences,sourceFrames=[old["source"],new["source"]],
                workbenchPayloads=[old["payload"],new["payload"]],tolerancesSI=limits,
                powerWithinOnePercent=abs(relative)<.01,
                fourSectionsWithinTolerance=all(s[k]["withinTolerance"] for s in sections for k in limits),
                matchedExperimentalAccuracy=False)
    write_json(out/"comparison.json",result)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    check(args.base.resolve(),args.new_output.resolve())
