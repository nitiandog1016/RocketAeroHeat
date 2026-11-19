"""从同一装配阶段的实际增量复算不可压试件的能量账；不启动流场求解。"""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens
from cm011Tunnel import D, RHO, U, CP, DT, GRIDS, meridian, volume
from readCM011Tunnel import evaluate, match_nodes, shell_balance
from runEmbeddedSU2 import sha

POWER=RHO*CP*U*DT*D**2
MASS=RHO*U*D**2
STAGES=("start","convection","diffusion","source","before_far","after_far","before_wall","after_wall","end")


def read(path):
    return json.loads(path.read_bytes())


def native_array(path):
    a=np.genfromtxt(path,delimiter=",",names=True)
    ids=a["id"].astype(int)
    assert np.array_equal(np.sort(ids),np.arange(len(a))) and np.array_equal(a["id"],ids)
    assert all(np.isfinite(a[k]).all() for k in a.dtype.names)
    return a[np.argsort(ids)]


def block_sums(mass,energy,r,href):
    weight=2*np.pi*r
    net=math.fsum(weight*mass)
    return dict(massKgS=net*MASS, rawEnthalpyW=math.fsum(weight*energy)*POWER,
                referenceAdvectionW=href*net*POWER,
                excessEnergyW=math.fsum(weight*(energy-href*mass))*POWER)


def near(a,b,tolerance=1e-12):
    if abs(a-b)>tolerance:
        raise ValueError(f"账项不闭合：{a:.15g} 与 {b:.15g} 相差 {a-b:.3g}")


def verify_sample(folder):
    info=read(folder/"sample.json");receipt=read(folder/"step_process.json")
    assert receipt["exitCode"]==0 and receipt["stopReason"] is None and receipt["livePids"]==[]
    for path,h in info["inputSHA256"].items(): assert sha(Path(path))==h
    for name,h in info["outputSHA256"].items(): assert sha(folder/name)==h
    assert len(list(csv.DictReader((folder/"history.csv").read_text().splitlines())))==1
    return info


def case(base,kind,origin):
    parent=ROOT/"outputs/v271_validation"/origin
    trace=base/(kind+"_trace");control=base/(kind+"_official")
    traced=verify_sample(trace);official=verify_sample(control)
    assert traced["trace"] and not official["trace"] and traced["source"]==official["source"]==origin
    parent_report=read(parent/"engine_process.json")
    for key in ("restart.csv","surface.vtk","vol_solution.vtk","history.csv"):
        assert sha(parent/key)==parent_report["outputSHA256"][key]
    before=sha(parent/"restart.csv")
    assert before==sha(trace/"initial_flow.csv")==sha(control/"initial_flow.csv")
    assert sha(trace/"mesh.su2")==sha(control/"mesh.su2")==sha(parent/"mesh.su2")
    assert (control/"case.cfg").read_bytes()==(trace/"case.cfg").read_bytes()
    if kind=="fine":
        off=base/"fine_off";verify_sample(off)
        for name in ("restart.csv","surface.vtk","vol_solution.vtk","history.csv"):
            assert sha(trace/name)==sha(off/name),name
    def restart(path):
        x=np.loadtxt(path,delimiter=",",skiprows=1)
        assert np.array_equal(x[:,0],np.arange(len(x))) and np.isfinite(x).all()
        return x
    a=restart(parent/"restart.csv")
    inputs=native_array(trace/"input_state.csv")
    replay=np.column_stack([inputs[k] for k in ("id","x","r","p","u","v","h")])
    assert np.max(abs(a[:,1:3]-replay[:,1:3]))<1e-12 and np.array_equal(a[:,3:],replay[:,3:])
    observed=restart(trace/"restart.csv");stock=restart(control/"restart.csv")
    after_difference=float(np.max(abs(observed[:,3:]-stock[:,3:])))
    assert after_difference < 2e-7
    # 保存场热量账直接用 v271 的完整终场，观测只读它作为第一个残差所用的输入。
    saved,_=evaluate(parent)
    body=read(parent/"body.json");level,layers,radius=GRIDS[body["grid"]]
    w,_=meridian(level);xy=volume(w,layers,radius,level)[0]
    xyz,fluid=vtk_tokens(parent/"vol_solution.vtk");order=match_nodes(xyz,xy)
    xy=xy.reshape(layers+1,len(w),2);vel=fluid["Velocity"][order,:2].reshape(xy.shape)
    temperature=fluid["Temperature"].ravel()[order].reshape(layers+1,len(w))
    old=shell_balance(xy,vel,temperature,layers)
    archived=read(ROOT/f"docs/validation/v271/evidence/{origin}.json")
    near(old["outwardHeatW"],archived["shellReconstruction"][-1]["outwardHeatW"])
    near(saved["curvedPowerW"],archived["curvedPowerW"])

    far=np.genfromtxt(trace/"far_terms.csv",delimiter=",",dtype=None,encoding="ascii",names=True)
    conv=far[far["term"]=="conv"];visc=far[far["term"]=="visc"]
    conv,visc=conv[np.argsort(conv["id"])],visc[np.argsort(visc["id"])]
    assert len(conv)==len(visc)==len(w) and np.array_equal(conv["id"],visc["id"])
    assert len(np.unique(conv["id"]))==len(conv) and len(far)==2*len(conv)
    for key in ("x","r","nx","nr","href","theta","u","v"):assert np.array_equal(conv[key],visc[key])
    ids=conv["id"]; href=float(conv["href"][0])
    assert np.all(conv["href"]==href) and np.all(visc["mass"]==0)
    np.testing.assert_allclose(conv["r"],inputs["r"][ids],rtol=0,atol=1e-12)
    np.testing.assert_allclose(conv["theta"],inputs["h"][ids]-href+1,rtol=0,atol=1e-10)
    conv_sum=block_sums(conv["mass"],conv["energy"],conv["r"],href)
    diff_sum=block_sums(-visc["mass"],-visc["energy"],visc["r"],href)
    native_far=conv_sum["excessEnergyW"]+diff_sum["excessEnergyW"]
    point_adv=math.fsum(2*np.pi*conv["r"]*(conv["nx"]*conv["u"]+conv["nr"]*conv["v"])*(conv["theta"]-1))*POWER

    stages={};prev=None
    for name in STAGES:
        cur=native_array(trace/("stage_"+name+".csv"))
        for key in ("id","x","r"):assert np.array_equal(cur[key],inputs[key])
        dm=cur["mass"].copy();de=cur["energy"].copy()
        if prev is not None: dm-=prev["mass"];de-=prev["energy"]
        stages[name]=block_sums(dm,de,cur["r"],href)
        if name=="after_far":
            np.testing.assert_allclose(dm[ids],conv["mass"]-visc["mass"],rtol=1e-11,atol=1e-12)
            np.testing.assert_allclose(de[ids],conv["energy"]-visc["energy"],rtol=1e-11,atol=1e-10)
            other=np.ones(len(cur),dtype=bool);other[ids]=False
            assert np.all(dm[other]==0) and np.all(de[other]==0)
        prev=cur
    total=block_sums(prev["mass"],prev["energy"],prev["r"],href)
    for k in ("massKgS","rawEnthalpyW","referenceAdvectionW","excessEnergyW"):
        near(math.fsum(s[k] for s in stages.values()),total[k])
    near(stages["after_far"]["excessEnergyW"],native_far)
    wall_in=stages["after_wall"]["excessEnergyW"]
    interior=math.fsum(stages[k]["excessEnergyW"] for k in ("convection","diffusion","source"))
    for name in ("start","before_far","before_wall","end"):near(stages[name]["excessEnergyW"],0)
    near(total["excessEnergyW"],interior+native_far+wall_in)
    gap=-saved["curvedPowerW"]-old["outwardHeatW"]
    terms=dict(wallQuadratureW=wall_in-saved["curvedPowerW"],
               nativeSurfaceImbalanceW=-wall_in-native_far,
               farAdvectionChoiceW=conv_sum["excessEnergyW"]-old["outwardSensibleAdvectionW"],
               farDiffusionChoiceW=diff_sum["excessEnergyW"]-old["outwardConductionW"])
    near(sum(terms.values()),gap)
    return dict(schema="cm011-native-energy-v1",grid=body["grid"],sourceRun=origin,
                sourceSHA256={k:sha(parent/k) for k in ("restart.csv","mesh.su2","case.cfg","surface.vtk","vol_solution.vtk")},
                sampledInitialSolutionExact=True, stockOneStepMaximumStateDifference=after_difference,
                observedSingleStepMaximumChange=float(np.max(abs(observed[:,3:]-a[:,3:]))),
                localObservationOffByteEquality=(True if kind=="fine" else None),referenceNativeEnthalpy=href,
                originalWallReleaseW=-saved["curvedPowerW"],reconstructed=old,nativeConvection=conv_sum,
                farAdvectionWithNodeNormalsW=point_adv,nativeDiffusion=diff_sum,nativeFarExcessOutW=native_far,
                nativeWallReleaseW=-wall_in,stageIncrements=stages,weightedAssembledResidual=total,
                innerFacePlusAxisSourceW=interior,reconstructedGapW=gap,gapTerms=terms,
                physicalValidation=False,newSolvedAreaM2=0)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path);parser.add_argument("new_output",type=Path)
    args=parser.parse_args();args.new_output.mkdir(parents=True,exist_ok=False)
    for k,n in (("fine","fine02"),("outer","outer02")):
        result=case(args.base.resolve(),k,n)
        with (args.new_output/(n+".json")).open("x",encoding="utf8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False);stream.write("\n")
