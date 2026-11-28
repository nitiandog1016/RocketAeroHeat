"""读回同一迭代的原生面/体场，核对热流和方向，并按原 STL 实际面积积分。"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from run_wall import BASE, DOC, geometry, read, save, sha
from plate_check import numbers
from auditSU2SSTState import vtk_tokens
from simpleLauncherTemperature import transfer

RES = ("rms[Rho]", "rms[RhoU]", "rms[RhoV]", "rms[RhoE]")
COLS = ("pressurePa","qIntoWallWm2","shearDownstreamPa","firstDownstreamMps","wallFluidK","firstK","yPlusFromShear")


def terminal(folder):
    run = read(folder/"debugger_process.json")
    assert run["status"]=="terminal" and run["debuggerExitCode"]==0 and not run["stopReason"] and not run["survivingProcessIDs"]
    assert all(run["inputsStable"].values()) and not run["supervisorError"]
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")
    return run


def mapping():
    model,tri,diam,wall=geometry()
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert len(area)==9216 and abs(sum(area)-model["areaM2"])<1e-12
    matrix=transfer(tri,wall,np.arange(len(wall)),diam)
    return model,tri,area,matrix,wall*diam


def diagnostics(folder,i,hist,wall,c):
    """同时刻核读终帧、给定壁温、壁节点流体温度、首层温度和物理正性。"""
    restart=folder/f"restart_{i:06d}.csv"
    surface=folder/f"surface_{i:06d}.vtk"
    a=numbers(restart)
    xyz,fields=vtk_tokens(surface)
    nw=len(wall)
    assert len(xyz)==nw and len(a["PointID"])==np.prod(c["shape"])
    np.testing.assert_array_equal(a["PointID"],np.arange(len(a["x"])))
    xy=np.column_stack((c["lengthM"]-a["x"],a["y"]))
    np.testing.assert_allclose(xy[:nw],wall,atol=5e-15,rtol=0)
    d,order=cKDTree(xyz[:,:2]).query(np.column_stack((a["x"][:nw],a["y"][:nw])))
    assert len(np.unique(order))==nw and max(d)<1e-6
    # 体场和壁面独立写出。比较导出精度后使用 ASCII 体场中较高精度的真实壁节点值。
    for native,volume in (("Pressure","Pressure"),("Heat_Flux","Heat_Flux"),("Y_Plus","Y_Plus"),
                          ("Temperature","Temperature"),("Skin_Friction_Coefficient","Skin_Friction_Coefficient_x")):
        np.testing.assert_allclose(fields[native][order,0],a[volume][:nw],rtol=6e-6,atol=1e-9)
    rho=a["Density"]
    u,v=a["Momentum_x"]/rho,a["Momentum_y"]/rho
    p=.4*(a["Energy"]-.5*rho*(u*u+v*v))
    np.testing.assert_allclose(a["Pressure"],p,rtol=1e-10)
    np.testing.assert_allclose(a["Temperature"],p/rho,rtol=1e-10)
    assert min(p)>0 and min(rho)>0 and min(a["Temperature"])>0
    assert max(abs(u[:nw]))<1e-12 and max(abs(v[:nw]))<1e-12
    temp=c["freeTemperatureK"];pinf=c["freePressurePa"];rhoinf=c["densityKgM3"]
    speed=c["speedMps"];dynamic=.5*rhoinf*speed**2
    np.testing.assert_allclose(a["Pressure_Coefficient"][:nw],(p[:nw]-1)*pinf/dynamic,atol=1e-10)
    edge=np.diff(wall,axis=0)
    normals=np.column_stack((-edge[:,1],edge[:,0]))
    dual=.5*np.vstack((normals[0],normals[:-1]+normals[1:],normals[-1]))
    length=np.linalg.norm(dual,axis=1)
    gap=np.einsum("ij,ij->i",xy[nw:2*nw]-xy[:nw],dual/length[:,None])
    assert min(gap)>0
    wall_temp=a["Temperature"][:nw]*temp
    near_temp=a["Temperature"][nw:2*nw]*temp
    mu=1.716e-5*(wall_temp/273.15)**1.5*(273.15+110.4)/(wall_temp+110.4)
    conductivity=mu*(1.4*287.05287/.4)/.72
    q=a["Heat_Flux"][:nw]
    reconstructed=conductivity*(near_temp-c["givenWallK"])/gap
    np.testing.assert_allclose(q,reconstructed,atol=2e-5,rtol=2e-8)
    cf=np.column_stack((a["Skin_Friction_Coefficient_x"][:nw],a["Skin_Friction_Coefficient_y"][:nw]))
    tau=cf*dynamic
    yplus=gap*np.sqrt(rho[:nw]*rhoinf*np.linalg.norm(tau,axis=1))/mu
    first_u=u[nw:2*nw]*np.sqrt(287.05287*temp)
    # 原生 HF 的权重是二维对偶面长度；另独立计算真实 STL 表面积分。
    line=float(q@length)
    np.testing.assert_allclose(line,hist["HF"][i],rtol=2e-8,atol=1e-5)
    data=np.column_stack((p[:nw]*pinf,q,tau[:,0],first_u,wall_temp,near_temp,yplus))
    assert np.isfinite(data).all()
    audit=dict(iteration=i,sourceRestartSHA256=sha(restart),sourceSurfaceSHA256=sha(surface),
        maximumFluxReadbackDifferenceWm2=float(np.max(abs(q-reconstructed))),
        rawNativeYPlusRange=[float(min(a["Y_Plus"][:nw])),float(max(a["Y_Plus"][:nw]))],
        nativeLineIntegralWm=float(hist["HF"][i]),independentLineIntegralWm=line,
        minPressurePa=float(min(p)*pinf),minDensityKgM3=float(min(rho)*rhoinf),
        minTemperatureK=float(min(a["Temperature"])*temp))
    return data,audit


def relative_rms(change,reference,area):
    return float(np.sqrt((area@(change**2))/(area@(reference**2))))


def inspect(folder,wall,matrix,area,cylinder):
    rc=terminal(folder)
    hist=numbers(folder/"history.csv")
    ids=hist["Inner_Iter"].astype(int)
    np.testing.assert_array_equal(ids,np.arange(len(ids)))
    last=int(ids[-1])
    volume=sorted(int(p.stem.split("_")[-1]) for p in folder.glob("restart_*.csv"))
    assert len(volume)>=3 and volume[-1]==last and sha(folder/"restart.csv")==sha(folder/f"restart_{last:06d}.csv")
    c=read(folder/"condition.json")
    np.testing.assert_allclose(c["meridianM"],wall,atol=1e-12)
    sampled=[]; windows=[]
    persistent_reverse=np.ones(len(area),dtype=bool)
    a=area[cylinder]
    for i in volume[-3:]:
        values,frame=diagnostics(folder,i,hist,wall,c)
        mapped=matrix@values
        reversed_nodes=(values[:,2]<0)&(values[:,3]<0)
        reverse=np.isclose(matrix@reversed_nodes.astype(float),1.,rtol=0,atol=1e-12)
        persistent_reverse &= reverse
        q=mapped[cylinder,1]
        power=a@q
        frame.update(cylinderIntoWallW=float(power),allFaceCandidateIntegralW=float(area@mapped[:,1]),
                     fullyReverseCylinderFaces=int(reverse[cylinder].sum()))
        if sampled:
            previous=sampled[-1][cylinder,1]
            frame["powerChangeFraction"]=float(abs(a@(q-previous))/abs(power))
            frame["heatChangeRmsFraction"]=relative_rms(q-previous,q,a)
        sampled.append(mapped)
        windows.append(frame)
    residual={k:float(max(hist[k][-20:])) for k in RES}
    gates=dict(residual=max(residual.values())< -8,thermalWindow=all(w["powerChangeFraction"]<.01 and
        w["heatChangeRmsFraction"]<.02 for w in windows[1:]),allWallYPlus=float(max(values[:,6]))<=1)
    receipt=dict(terminal=True,iterations=len(ids),savedIterations=volume,seconds=rc["elapsedSeconds"],
                 peakRSSBytes=rc["peakFamilyRSSBytes"],last20Max=residual,windows=windows,gates=gates,
                 wallMaxDerivedYPlus=float(max(values[:,6])),meshSHA256=sha(folder/"mesh.su2"),
                 configSHA256=sha(folder/"case.cfg"),historySHA256=sha(folder/"history.csv"))
    return values,sampled[-1],receipt,persistent_reverse


def write_report(output):
    model,tri,area,matrix,wall=mapping()
    cyl,=[r for r in model["regions"] if r["name"]=="cylinder"]
    sel=slice(cyl["firstFace"]-1,cyl["lastFace"])
    a=area[sel]
    runs={};final=[];nodal=[];persistent=[]
    for cold,first,warm in (("normal4_02","first4","warm4"),("normal2_01","first2","warm2")):
        _,_,runs[cold],_=inspect(BASE/cold,wall,matrix,area,sel)
        terminal(BASE/first)
        h0=numbers(BASE/first/"history.csv");h1=numbers(BASE/cold/"history.csv")
        assert len(h0["Inner_Iter"])==1 and all(h0[k][0]==h1[k][0] for k in h1)
        _,first_frame=diagnostics(BASE/first,0,h0,wall,read(BASE/first/"condition.json"))
        runs[first]=first_frame
        node,face,runs[warm],held=inspect(BASE/warm,wall,matrix,area,sel)
        persistent.append(held)
        assert sha(BASE/warm/"seed.csv")==sha(BASE/cold/"restart_002499.csv")
        nodal.append(node);final.append(face)

    condition=read(BASE/"warm2/condition.json")
    dynamic=.5*condition["densityKgM3"]*condition["speedMps"]**2
    p4,q4=final[0][sel,0:2].T;p2,q2=final[1][sel,0:2].T
    diffs=dict(powerRelative=float(abs(a@(q4-q2))/abs(a@q2)),heatRelativeRms=relative_rms(q4-q2,q2,a),
               pressureRmsOverDynamic=float(np.sqrt(np.average((p4-p2)**2,weights=a))/dynamic))
    diffs["withinGridThresholds"]=diffs["powerRelative"]<.05 and diffs["heatRelativeRms"]<.05 and diffs["pressureRmsOverDynamic"]<.05
    # 用两端节点同时判定全带反向，跨过变号点的三角面不计。面内线性插值不能证明真实分离。
    sides=[]; reverse_sets=[];merged={}
    for g,node,face in zip(("normal4","normal2"),nodal,final):
        rev=(node[:,2]<0)&(node[:,3]<0)
        mapped_rev=np.isclose(matrix@rev.astype(float),1.,rtol=0,atol=1e-12)
        both=np.zeros(len(area),bool);both[sel]=mapped_rev[sel]
        reverse_sets.append(both)
        merged[g]=dict(cylinderPowerW=float(a@face[sel,1]),reverseCylinderFaces=int(both.sum()),
            reverseCylinderAreaM2=float(area[both].sum()),minHeatWm2=float(min(face[sel,1])),
            maxHeatWm2=float(max(face[sel,1])),maxCylinderFaceYPlus=float(max(face[sel,6])))
        sides.append(face)
    joint=reverse_sets[0]&reverse_sets[1]
    repeated=joint & persistent[0] & persistent[1]
    # 按同一母线插值权重保存面组，使侧栏点选能还原两档数值。
    weights,face_group=np.unique(matrix.toarray(),axis=0,return_inverse=True)
    values=[weights@n for n in nodal]
    for j in range(2):
        np.testing.assert_allclose(values[j][face_group],sides[j],rtol=5e-13,atol=1e-10)
    archives=dict(mach=.3,altitudeM=5000,wallGivenK=300,surfaceSHA256=model["fileSHA256"],
        faces=len(area),areaM2=float(area.sum()),columns=list(COLS),faceGroup=(face_group+1).tolist(),
        normal4=values[0].tolist(),normal2=values[1].tolist(),cylinder=merged,gridDifference=diffs,
        residualPassedBoth=all(runs[k]["gates"]["residual"] for k in ("warm4","warm2")),
        heatWindowPassedBoth=all(runs[k]["gates"]["thermalWindow"] for k in ("warm4","warm2")),
        commonReverseCylinderFaces=int(joint.sum()),commonReverseCylinderAreaM2=float(area[joint].sum()),
        threeFrameBothGridReverseFaces=int(repeated.sum()),threeFrameBothGridReverseAreaM2=float(area[repeated].sum()),
        wholeBodyPowerW=None,physicalQualification=False)
    intervals=np.flatnonzero(np.diff(np.r_[False,joint,False]))
    archives["commonReverseCylinderFaceRanges"]=[[int(lo+1),int(hi)] for lo,hi in intervals.reshape(-1,2)]
    if joint.any():
        archives["commonReverseXRangeM"]=[float(tri[joint,:,0].min()),float(tri[joint,:,0].max())]
    output.mkdir(parents=True,exist_ok=False)
    (output/"surface.json").write_bytes((json.dumps(archives,separators=(",",":"),allow_nan=False)+"\n").encode())
    save(output/"runs.json",runs)
    vector=np.column_stack([wall,*nodal])
    names=["xM","rM",*[f"{key}_{g}" for g in ("4um","2um") for key in COLS]]
    np.savetxt(output/"wall.csv",vector,delimiter=",",header=",".join(names),comments="")
    print(json.dumps({"cylinder":merged,"grid":diffs,"sharedReverse":int(joint.sum())},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    write_report(p.parse_args().out)
