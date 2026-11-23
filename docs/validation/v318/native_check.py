"""用真实修正壁温推进热敏珠标量方程，另外保存对照所用输入与原生终态。"""
import argparse
from math import asin, log, pi, sqrt
import os
from pathlib import Path
import shutil
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha,write_json
from thermistorMeasuredCold import HEADER
from thermistorWireCase import native
from thermistorWireSurface import collect


def appendix(c):
    # 以共焦椭球解析面积反算焦半距。对数取来自壳内径向导热积分的自然对数。
    xi_w,xi_i=1.80,1.25
    length=.3430e-3
    area=c["sourceAreaM2"]
    assert area==3.094e-6
    shape=xi_w**2-1+xi_w**2*sqrt(xi_w**2-1)*asin(1/xi_w)
    focus=sqrt(area/(2*pi*shape))
    reynolds=c["speedMps"]*length*c["rho"]/c["mu"]
    k_s=8.2*4184/3600
    ratio=focus*c["k"]/(4*length*k_s)*(2/3+.3*sqrt(reynolds))*shape*log((xi_i+1)*(xi_w-1)/((xi_i-1)*(xi_w+1)))
    return dict(xiWall=xi_w,xiInternal=xi_i,lengthM=length,areaM2=area,semifocalM=focus,
                reynoldsL=reynolds,solidConductivityWmK=k_s,ratio=ratio,factor=1/(1+ratio))


def run(base,build,out):
    if shutil.disk_usage(out.parent).free<6*1024**3:raise OSError("预留 6 GiB 数值检查空间")
    done=read(base/"completed.json");source=base/done["stage"]
    assert sha(source/"restart.csv")==done["frame"]=="a98cfbb17a2061e87dccd69fe0c1ebd712933e6f34b510ecb56ad2f11f0df5db"
    src=base/"prepared";c=read(src/"case.json");s=read(base/"surface/summary.json")
    assert s["sourceFrameSHA256"]==done["frame"] and sha(base/"surface/wall.csv")==s["wallSHA256"]
    assert max(s["window"].values())<-8 and abs(s["energy"]["relativeDifference"])<.005
    model=appendix(c);factor=model["factor"];air=c["airK"];wall=c["wallK"]
    assert air==293.15 and wall==305.15 and c["scalarReferenceDeltaK"]==24
    out.mkdir(parents=True,exist_ok=False);prep=out/"prepared";prep.mkdir()
    for n in ("mesh.su2","nodes.npz","body.json"):os.link(src/n,prep/n)
    c["wallK"]=air+(wall-air)*factor
    write_json(prep/"case.json",c)
    profile=read(src/"profile.json")
    for n in ("uniform","profile","release"):
        table=np.loadtxt(src/(n+".dat"))
        table[:,1]=1+factor*(table[:,1]-1)
        table[:,2]*=factor
        np.savetxt(prep/(n+".dat"),table,fmt=["%d","%.17g","%.17g"])
    frame=np.loadtxt(source/"restart.csv",delimiter=",",skiprows=1)
    frame[:,9]*=factor
    env_h=1-298.15
    frame[:,8]=env_h+factor*(frame[:,8]-env_h)
    np.savetxt(prep/"initial_flow.csv",frame,delimiter=",",header=HEADER,comments="",fmt=["%d"]+["%.17g"]*9)
    profile["endK"]=air+factor*(profile["endK"]-air)
    for k in ("oneLeadRootW","oneLeadThroughCutW","oneLeadSideW"):profile[k]*=factor
    profile["inputs"]={n:sha(prep/n) for n in profile["inputs"]}
    profile["temperatureReconstruction"]=model
    write_json(prep/"profile.json",profile)
    write_json(out/"request.json",dict(model=model,internalK=wall,surfaceK=c["wallK"],ambientK=air,wireFrame=done["frame"],
               uniformControl=profile["baseFrame"],initialSHA256=profile["inputs"]["initial_flow.csv"],
               sourceFiles={p.relative_to(ROOT).as_posix():sha(p) for p in (HERE/"PLAN.md",Path(__file__).resolve(),ROOT/"scripts/thermistorWireCase.py",
               ROOT/"scripts/thermistorWireSurface.py",base/"completed.json",source/"native.json",build/"build.json")}))
    stage=out/"steady01"
    record=native(prep,build,stage,prep/"initial_flow.csv",20,"profile")
    assert record["passed"] and record["steps"]==20
    collect(prep,build,stage,out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("base","build","out"):p.add_argument(k,type=Path)
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs");run(a.base.resolve(),a.build.resolve(),out)
