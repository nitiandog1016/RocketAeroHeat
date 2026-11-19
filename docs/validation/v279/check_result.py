"""独立复算给定边界的工程分布和功率，并逐行核对工作台导出。"""
import argparse
import csv
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def values(angle_deg, width_deg, c):
    # 不使用 MATLAB 的多项式原函数，直接在实际角区间独立积分。
    def cp(a):
        x=np.deg2rad(a)/2
        v=x*(3.6314-x*x*(2.1709+1.5144*x*x))
        return 1-v*v
    def nu(a):
        x=np.deg2rad(a)/2
        return np.sqrt(177000)*(.9449-.5100*x*x-.5956*x**4)
    bounds=[angle_deg-width_deg/2,angle_deg+width_deg/2]
    p=c["referencePressurePa"]+c["dynamicPressurePa"]*quad(cp,*bounds)[0]/width_deg
    q=quad(nu,*bounds)[0]/width_deg*c["conductivityWmK"]/.066*(c["adiabaticWallK"]-c["wallTemperatureK"])
    return p,q


def check(gui,out,modules=None):
    out.mkdir(exist_ok=False,parents=True)
    page=gui/"workbench"
    model=json.loads((page/"case-data.js").read_text("utf8").strip().split("=",1)[1][:-1])
    rec=json.loads((gui/"workflow.json").read_bytes());c=rec["condition"]
    # 用温度为未知量闭合质量通量方程，独立于生产代码以速度为未知量的算法。
    rho=101325/(287.05*293.15)
    def speed(t):
        mu=1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
        return 177000*mu/(rho*.066)
    def reference(t):
        tad=293.15+np.sqrt(.7)*speed(t)**2/(2*1005)
        return .5*303.15+.28*293.15+.22*tad
    t=brentq(lambda t:reference(t)-t,290,310)
    assert abs(t-c["referenceTemperatureK"])<1e-10 and abs(speed(t)-c["velocityMps"])<1e-10
    fields={f["key"]:f for f in model["fields"]}
    v=np.array(model["vertices"]);tri=v[model["faces"]]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    angle=np.rad2deg(np.unwrap(np.arctan2(tri[:,:,2],tri[:,:,1]),axis=1))
    theta=(angle.min(axis=1)+.5+180)%360-180
    select=(np.abs(theta)<55)&(np.max(abs(tri[:,:,0]),axis=1)<=.033)
    expected=np.array([values(a,1.,c) for a in theta[select]])
    diffs={}
    for field,key in enumerate(("pressure","heatFlux")):
        f=fields[key]
        assert np.array_equal(np.array(f["status"])==1,select) and all(f["values"][i] is None for i in np.flatnonzero(~select))
        shown=np.asarray([f["values"][i]*1000 for i in np.flatnonzero(select)])
        diffs[key]=float(max(abs(shown-expected[:,field])))
        np.testing.assert_allclose(shown,expected[:,field],rtol=0,atol=1e-8)
    power=float(area[select]@expected[:,1])
    assert abs(power-rec["hot"]["calculatedHeatPowerW"])<1e-10 and sum(select)==220
    assert model["wallCoverage"]["wholeSurfaceHeatPowerW"] is None and not model["wallCoverage"]["fullSurfaceComputed"]
    assert set(fields["wallInput"]["status"])=={2} and set(fields["wallSolved"]["status"])=={0}
    assert model["resultCondition"]["wallUniformValueK"]==303.15 and rec["cold"]["calculatedHeatPowerW"]>0
    for key in ("trialPressure","trialHeatFlux"):
        assert fields[key]["computedStateLabel"]=="未收敛试算" and set(fields[key]["status"])=={1}
    trial=np.asarray(fields["trialHeatFlux"]["values"])*1000
    stations=np.genfromtxt(gui/"stations.csv",delimiter=",",names=True,encoding="utf-8-sig")
    check_stations=np.array([values(row["angleDeg"],row["upperDeg"]-row["lowerDeg"],c) for row in stations])
    np.testing.assert_allclose(stations["pressurePa"],check_stations[:,0],atol=1e-8,rtol=0)
    np.testing.assert_allclose(stations["heatFluxWm2"],check_stations[:,1],atol=1e-8,rtol=0)
    # 计算和诊断共用相同截面，CSV 中两种结果有各自的状态、适用说明和原面编号。
    code="""
const fs=require('fs'),p=require('path'),[page,mods,out]=process.argv.slice(1);
global.window={};require(p.join(page,'case-data.js'));const m=window.ROCKET_CASE;
const s=require(p.join(mods,'section.js')),v=require(p.join(mods,'viewer.js')),cut=s.geometry(m.vertices,m.faces,0);
if(cut.components.length!==1||!cut.components[0].closed||cut.components[0].invalidAngleSegments)throw Error('section');
for(const key of ['pressure','heatFlux','trialHeatFlux']){
  const f=m.fields.find(k=>k.key===key);if(!f.scopeNames.length)f.scopeNames=null;
  fs.writeFileSync(p.join(out,key+'.csv'),s.csv(cut,f,'冻结工况',v.caseCaption(m,f,0),{model:m}),{flag:'wx'});
}
"""
    subprocess.run(["node","-e",code,str(page.resolve()),str((modules or page).resolve()),str(out.resolve())],check=True)
    counts={}
    for key in ("pressure","heatFlux","trialHeatFlux"):
        with (out/(key+".csv")).open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        f=fields[key];assert len(rows)==len({r["面元编号"] for r in rows})==720
        known=0
        for row in rows:
            i=int(row["面元编号"])-1;value=row[f["label"]+" / "+f["unit"]]
            assert row["算法"]==f["methodNames"][f["methodIndex"][i]] and "303.15" in row["工况"]
            assert min(tri[i,:,0])<0<max(tri[i,:,0])
            if f["status"][i]:
                assert float(value)==f["values"][i];known+=1
            else:assert value==""
            if key.startswith("trial"):assert row["状态"]=="未收敛试算"
        counts[key]=known
    assert counts==dict(pressure=220,heatFlux=220,trialHeatFlux=720)
    ans=dict(knownFaceCount=int(sum(select)),areaM2=float(area[select].sum()),partialPowerW=power,
             projectedTrialPowerW=float(area@trial),wholePhysicalPowerW=None,sectionKnownSegments=counts,
             maxIndependentDifferenceSI=diffs,pressureRangePa=[expected[:,0].min(),expected[:,0].max()],
             heatRangeWm2=[expected[:,1].min(),expected[:,1].max()])
    (out/"checks.json").write_text(json.dumps(ans,indent=2),encoding="utf8")
    plot(stations,out/"forward.png")


def plot(stations,path):
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
    fig,axes=plt.subplots(1,2,figsize=(8.2,3.75))
    for ax,key,unit,color in zip(axes,("pressurePa","heatFluxWm2"),("绝对壁压 / kPa","入壁热流 / kW/m²"),("#245980","#a15c35")):
        ax.plot(stations["angleDeg"],stations[key]/1000,"o-",color=color)
        ax.set_xlabel("测带中心距驻点的角度 / °");ax.set_ylabel(unit);ax.set_xticks(np.arange(0,51,10))
        ax.grid(axis="y",color="#dce1e1");ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
    fig.suptitle("给定边界的 66 mm 圆柱 · 中央六条测带的工程预测")
    fig.tight_layout(rect=(.01,.01,.99,.94));fig.savefig(path,dpi=155);plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gui",type=Path);parser.add_argument("new_output",type=Path)
    parser.add_argument("--web-modules",type=Path)
    args=parser.parse_args();check(args.gui,args.new_output,args.web_modules)
