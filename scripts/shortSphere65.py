"""计算 Short 第 65 次银球试验的全周稳态参考场。大体场与原生终点保存在单独任务目录。"""
import argparse
import csv
import os
from pathlib import Path

import numpy as np

from cm011TunnelSector import mesh, write_volume
from runSU2NativeDebugger import observe
from thermistorBead import ROOT, read, sha, write_json
from thermistorMeasuredCold import config, HEADER, RESIDUALS

CASE = ROOT/"cases/short_sphere/test65.json"
DRIVER = {"wire":"9721e25124a866652e8d004b531854ec88d2f856e99ca7b80c0ab34d2bbebb55",
          "wire_observer":"0adeb7792046d6525ba4fc8a4c1fecdbe5b5704c43aac5de24662c987dcbb6c4"}


def kelvin(f):
    return (np.asarray(f)-32)*5/9+273.15


def conductivity(t):
    return .0241*(t/273.15)**1.5*(273.15+194)/(t+194)


def condition(c):
    d=c["diameterIn"]*.0254;u=c["speedFtS"]*.3048;t=kelvin(c["airF"]);tw=kelvin(c["meanSurfaceF"])
    p=c["pressureLbfFt2"]*4.4482216152605/.3048**2;y=c["waterMassFraction"]
    rho=p/(t*((1-y)*287.05+y*461.52));cp=(1-y)*1006+y*1860
    k=conductivity((t+tw)/2);re=c["reynoldsFree"]
    return dict(diameterM=d,speedMps=u,airK=float(t),airPa=p,rho=rho,cp=cp,kFilm=k,
                mu=rho*u*d/re,re=re,pe=rho*cp*u*d/k,scalarUnitK=90.)


def prepare(out):
    # 模块一：球形外流域、用户面号、真实首层和实测壁温的极角延拓。
    out.mkdir(parents=True,exist_ok=False);c=read(CASE);s=condition(c)
    m=48;nt=24;radial=32;t=np.linspace(np.pi,0,m+1)
    wall=.5*np.column_stack((np.cos(t),np.sin(t)));wall[[0,-1],1]=0
    _,xyz,_,_,ids,cells,markers,geom=mesh(np.linspace(0,2*np.pi,nt+1),2,radial,wall,14.,s["diameterM"])
    write_volume(out,xyz,cells,markers)
    wet=np.unique(ids[:,:m+1]);wet_lookup={v:i for i,v in enumerate(wet)}
    first=np.empty_like(wet)
    for ring in ids[:-1]:
        for v,near in zip(ring[:m+1],ring[m+1:2*m+2]):first[wet_lookup[v]]=near
    vertices=xyz[wet].tolist();faces=[];centers=[]
    for boundary in markers["wall"]:
        nodes=[wet_lookup[v] for v in reversed(boundary)]
        if len(nodes)==3:faces.append(nodes);continue
        middle=len(vertices);vertices.append(np.mean(xyz[wet[nodes]],axis=0));centers.append(nodes)
        faces.extend([nodes[i],nodes[(i+1)%4],middle] for i in range(4))
    v=np.asarray(vertices)*s["diameterM"];f=np.array(faces);tri=v[f]
    vector=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    assert np.all(np.einsum("ij,ij->i",vector,tri.mean(axis=1))>0)
    area=np.linalg.norm(vector,axis=1)/2;assert abs(area.sum()/geom["chordWallAreaM2"]-1)<1e-12
    centers=np.array(centers);lump=np.bincount(f.ravel(),np.repeat(area/3,3),minlength=len(v))
    lump=lump[:len(wet)]+np.bincount(centers.ravel(),np.repeat(lump[len(wet):]/4,4),minlength=len(wet))
    angles=np.degrees(np.arctan2(np.linalg.norm(xyz[wet,1:],axis=1),xyz[wet,0]))
    temps=np.interp(angles,c["angleDeg"],kelvin(c["surfaceF"]));theta=(temps-s["airK"])/s["scalarUnitK"]
    assert np.all((theta>0)&(theta<.5))
    for name,scalar in [("profile",theta),("release",np.zeros(len(wet)))]:
        np.savetxt(out/(name+".dat"),np.column_stack((wet,1+theta,scalar)),fmt=["%d","%.17g","%.17g"])
    state=np.zeros((len(xyz),6));state[:,1]=-1;state[wet,1]=0;state[:,4]=1-298.15
    state[wet,4]+=theta;state[wet,5]=theta
    np.savetxt(out/"initial.csv",np.column_stack((np.arange(len(xyz)),xyz,state)),delimiter=",",
               fmt=["%d"]+["%.17g"]*9,header=HEADER,comments="")
    line=np.searchsorted(wet,ids[nt//4,:m+1][::-1])
    np.savez_compressed(out/"mapping.npz",xyz=xyz,wet=wet,first=first,V=v,F=f,centers=centers,
                        area=area,lump=lump,angles=angles,line=line,theta=theta)
    files=["mesh.su2","mapping.npz","initial.csv","profile.dat","release.dat"]
    write_json(out/"prepared.json",dict(condition=s,geometry=geom,files={n:sha(out/n) for n in files},
                                        caseSHA256=sha(CASE),faces=len(f),cells=len(cells)))


def configuration(s,steps,release,species_convection="SCALAR_UPWIND"):
    # 模块二：固定求解器与原标量方程，采用本试件的 Re、Pe、轴向来流和逐节点热界。
    if species_convection not in ("SCALAR_UPWIND","BOUNDED_SCALAR"):
        raise ValueError("未知标量对流配置")
    changes={"INC_VELOCITY_INIT":"(-1, 0, 0)","CONV_NUM_METHOD_SPECIES":species_convection,"MU_CONSTANT":1/s["re"],"THERMAL_CONDUCTIVITY_CONSTANT":1/s["pe"],
             "DIFFUSIVITY_CONSTANT":1/s["pe"]}
    if release:changes.update(MARKER_WALL_SPECIES="(wall, FLUX, 0.0)",VOLUME_OUTPUT="(COORDINATES, SOLUTION, RESIDUAL)")
    lines=[]
    for line in config(steps).splitlines():
        key=line.partition("=")[0].strip()
        lines.append(f"{key}= {changes.pop(key)}" if key in changes else line)
    assert not changes
    return "\n".join(lines)+"\nMARKER_PYTHON_CUSTOM= (wall)\n"+("WRT_RESTART_COMPACT= NO\n" if release else "")


def advance(prep,dest,build,start,steps,mode,species_convection="SCALAR_UPWIND"):
    p=read(prep/"prepared.json");s=p["condition"]
    dest.mkdir();release=mode=="release";trace=mode=="observer"
    exe=build/("wire_observer.exe" if trace else "wire.exe")
    assert sha(exe)==DRIVER[exe.stem] and p["caseSHA256"]==sha(CASE)
    for src,name in [(prep/"mesh.su2","mesh.su2"),(prep/("release.dat" if release else "profile.dat"),"wire_boundary.dat"),(start,"initial_flow.csv")]:
        os.link(src,dest/name)
    (dest/"case.cfg").write_text(configuration(s,steps,release,species_convection),encoding="ascii")
    locks={str(q):sha(q) for q in [exe,Path(__file__).resolve(),CASE,prep/"prepared.json",*(dest/n for n in ("case.cfg","mesh.su2","wire_boundary.dat","initial_flow.csv"))]}
    if trace:os.environ["RAH_CM011_LOW_SPEED_LEDGER"]="1"
    else:os.environ.pop("RAH_CM011_LOW_SPEED_LEDGER",None)
    observe([str(exe),"case.cfg","--threads",str(4 if steps>20 else 1)],dest,locks,
            300 if steps>20 else 150,(450 if steps>20 else 300)*1024**2)
    r=read(dest/"debugger_process.json")
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
    with (dest/"history.csv").open() as h:
        keys=[k.strip() for k in next(csv.reader(h,skipinitialspace=True))];history=np.loadtxt(h,delimiter=",",ndmin=2)
    count=len(history);assert 0<count<=steps and np.isfinite(history).all()
    np.testing.assert_array_equal(history[:,keys.index("Inner_Iter")],np.arange(count))
    window={k:float(history[-20:,keys.index(k)].max()) for k in RESIDUALS}
    passed=count>=20 and max(window.values()) < -8
    assert count==steps or passed
    with np.load(prep/"mapping.npz") as z:xyz,wet,theta=z["xyz"],z["wet"],z["theta"]
    assert f"Applied thermistor wall nodes: {len(wet)}" in (dest/"debugger.log").read_text("ascii",errors="replace")
    frames={}
    for n in sorted({count-1,*range(100,count,100)}):
        path=dest/f"restart_{n:06d}.csv";state=np.loadtxt(path,delimiter=",",skiprows=1,usecols=tuple(range(10)))
        assert np.isfinite(state).all() and np.min(s["airPa"]+s["rho"]*s["speedMps"]**2*state[:,4])>0
        np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-12)
        if not release:
            np.testing.assert_allclose(state[wet,9],theta,atol=2e-14,rtol=0)
            assert state[:,9].min()>=-1e-12 and state[:,9].max()<=theta.max()+1e-12
            assert np.max(abs(state[wet,5:8]))<1e-10
        frames[path.name]=sha(path)
    end=sha(dest/"restart.csv");assert end==frames[path.name]
    done=dict(steps=count,passed=passed,window=window,seconds=r["elapsedSeconds"],peakRSS=r["peakFamilyRSSBytes"],
              frames=frames,frame=end,input=sha(start),history=sha(dest/"history.csv"))
    write_json(dest/"native.json",done);return done


def run(work,build):
    work.mkdir(parents=True,exist_ok=False)
    if build is None:
        build=ROOT/"outputs/v317_validation/build01"
        if not all((build/(k+".exe")).is_file() and sha(build/(k+".exe"))==h for k,h in DRIVER.items()):
            from buildThermistorWireDriver import build as build_driver
            build=work/"driver";build_driver(build)
    prep=work/"prepared";prepare(prep)
    p=read(prep/"prepared.json");assert all(sha(prep/n)==h for n,h in p["files"].items())
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build"))
    os.environ["PATH"]=str(kit/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
    start=prep/"initial.csv";trend=[]
    for i in range(1,5):
        stage=work/f"steady{i:02d}";done=advance(prep,stage,build,start,400,"solve");trend.append(max(done["window"].values()))
        print(stage.name,trend[-1],flush=True)
        if done["passed"]:break
        if len(trend)>=3 and trend[-1]>=trend[-2]>=trend[-3]:raise ValueError("连续两段末窗不改善")
        start=stage/"restart.csv"
    else:raise ValueError("基流预算已用完，尚未达到残差门槛")
    end=stage/"restart.csv"
    for name,seed,count,method in [("reaction",end,1,"release"),("plus20",end,20,"solve"),
             ("reaction20",stage/"plus20/restart.csv",1,"release"),("control",end,1,"solve"),("observer",end,1,"observer")]:
        receipt=advance(prep,stage/name,build,seed,count,method)
        if name=="plus20":assert receipt["passed"]
    write_json(work/"endpoints.json",dict(stage=stage.name,frame=sha(end),case=sha(CASE),residuals=trend))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path)
    parser.add_argument("--driver",type=Path)
    a=parser.parse_args();out=a.work.resolve();out.relative_to(ROOT/"outputs")
    run(out,a.driver.resolve() if a.driver else None)
    from shortSphereSurface import recover
    recover(out,out/"surface")
    write_json(out/"completed.json",dict(summary=sha(out/"surface/summary.json"),archive=sha(out/"surface/surface.zip"),
                                       endpoints=sha(out/"endpoints.json"),case=sha(CASE)))
