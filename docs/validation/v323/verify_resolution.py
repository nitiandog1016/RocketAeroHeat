"""以球心距离和已知角度重建交点，独立复核银球诊断、支撑节点及工作台显示。"""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from shortSphereWakeResolution import axis_result, recovery_band, section


def read(path):
    return json.loads(path.read_bytes())


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-11)


def manufactured():
    # 人工速度 = 2x+3z；射线交点必须准确复现一次场，并采用 −x 下游约定。
    xyz = np.array([[-.5,0,0],[-2,0,0],[-.5,0,.2],[-2,0,.8],[-.5,0,.3],[-2,0,1.2]])
    u = np.zeros_like(xyz)
    u[:,0] = 2*xyz[:,0]+3*xyz[:,2]
    cut = section(-1,xyz,np.arange(6).reshape(3,2),u)
    close(cut["radiusD"],[0,.4,.6])
    close(cut["downstreamSpeedOverU"],[2,.8,.2])
    rad = np.array([0,.2,.4,.6,1])
    band = recovery_band(rad,rad,np.ones(5)*.25)
    close([band["r10D"],band["r90D"],band["widthD"]],[.1,.9,.8])
    assert band["innerRays"]==3 and band["sparse"]
    assert recovery_band(rad,np.ones(5),rad)["status"]=="no-deficit"
    assert recovery_band(rad,np.zeros(5),rad)["status"]=="unbracketed-recovery"
    xyz = np.zeros((4,3));xyz[:,0] = [-.5,-.6,-.7,-.8]
    u = np.zeros_like(xyz);u[:,0]=[0,.2,.1,-.2]
    close(axis_result(xyz,np.arange(4),u)["lengthD"],.2+.1/3)
    assert axis_result(xyz,np.arange(4),-u)["status"]=="no-wall-connected-reversal"


def verify_branch(d, branch, work):
    n = {"coarse":32, "fine":64}[branch]
    frame = work/branch/f"step{n:03d}"/f"restart_{n:05d}.csv"
    assert sha(frame)==d["frameSHA256"]
    a = np.loadtxt(frame, skiprows=1, delimiter=",", usecols=(1,2,3,5,6,7))
    xyz, u = a[:,:3], a[:,3:]
    tree = cKDTree(xyz)
    axis = np.flatnonzero((xyz[:,0]<0)&(np.linalg.norm(xyz[:,1:],axis=1)<1e-13))
    axis = axis[np.argsort(-xyz[axis,0])]
    close(axis,d["axis"]["nodeIDs"])
    radius = -xyz[axis,0]
    bubble = np.flatnonzero((u[axis[:-1],0]>0)&(u[axis[1:],0]<=0))[0]
    bracket = radius[bubble:bubble+2]
    root = brentq(lambda r: np.interp(r,radius,u[axis,0]),*bracket)
    close(bracket-.5,d["axis"]["bracketD"])
    close(root-.5,d["axis"]["lengthD"])
    close(d["axis"]["downstreamSpeedOverU"],-u[axis,0])
    assert d["axis"]["reverseNodes"]==bubble and np.all(u[axis[1:bubble+1],0]>0)
    # 从独立的等角母线及轴线实测球壳半径出发，三维空间寻找每个原始网格节点。
    angles = np.linspace(0,np.pi,49)
    directions = np.c_[-np.cos(angles),np.zeros(49),np.sin(angles)]
    distances,all_ids = tree.query(directions[:,None,:]*radius[None,:,None])
    assert distances.max()<1e-12
    comparison=[]
    for s in d["sections"]:
        x=s["xD"]
        rays=np.flatnonzero(directions[:,0]<x/radius[-1])
        target=x/directions[rays,0]
        j=np.searchsorted(radius,target)
        ids=all_ids[rays]
        np.testing.assert_array_equal(ids[np.arange(len(ids)),j-1],s["sourceNode0"])
        np.testing.assert_array_equal(ids[np.arange(len(ids)),j],s["sourceNode1"])
        r=target*directions[rays,2]
        speed=np.array([-np.interp(t,radius,u[row,0]) for t,row in zip(target,ids)])
        dx=(radius[j]-radius[j-1])*(-directions[rays,0])
        close(r,s["radiusD"]);close(speed,s["downstreamSpeedOverU"]);close(dx,s["axialSpanD"])
        band=s["band"]
        core=np.flatnonzero(r<=.6);i=core[np.argmin(speed[core])]
        levels=speed[i]+np.array([.1,.9])*(1-speed[i])
        bounds=[]
        for level in levels:
            k=i+np.flatnonzero(speed[i:]>=level)[0]
            bounds.append(brentq(lambda v: np.interp(v,r,speed)-level,r[k-1],r[k]))
        r10,r90=bounds;inside=(r>r10)&(r<r90)
        lo=np.searchsorted(r,r10)-1;hi=np.searchsorted(r,r90)
        gaps=np.diff(r[lo:hi+1])
        close([r10,r90,min(gaps),max(gaps),min(dx[lo:hi+1]),max(dx[lo:hi+1])],
              [band[k] for k in ("r10D","r90D","minRadialGapD","maxRadialGapD","minAxialSpanD","maxAxialSpanD")])
        close(band["widthD"],r90-r10)
        ratio=(r90-r10)/max(gaps)
        close(ratio,band["widthOverMaxGap"])
        assert int(sum(inside))==band["innerRays"] and bool(sum(inside)<4 or ratio<4)==band["sparse"]
        comparison.append(dict(xD=x, innerRays=int(sum(inside)), widthOverMaxGap=ratio, sparse=band["sparse"]))
    polar=np.pi-angles[::-1]
    first=all_ids[::-1,1]
    vt=np.sum(u[first]*np.c_[-np.sin(polar),np.zeros(len(polar)),np.cos(polar)],axis=1)
    for station in d["stations"]:
        ang=np.radians(station["angleDeg"])
        idx=np.searchsorted(polar,ang)
        np.testing.assert_array_equal(first[idx-1:idx+1],station["adjacentFirstNodeIDs"])
        close(np.interp(ang,polar,vt),station["firstTangentialSpeedOverU"])
    return dict(frame=d["frameSHA256"],axisLengthD=root-.5,axisBracketD=(bracket-.5).tolist(),sections=comparison)


def plot_resolution(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.spines.top":False,"axes.spines.right":False})
    fig,(a,b)=plt.subplots(1,2,figsize=(9,3.45),constrained_layout=True)
    axis=data["axis"];mm=12.7
    a.axvspan(*(mm*np.array(axis["bracketD"])),alpha=.16,color="#b46a32",label="零速点所在网格间隔")
    a.plot(mm*np.array(axis["distanceD"]),axis["downstreamSpeedOverU"],"-o",ms=3,color="#2e5579",label="原轴线节点")
    a.axhline(0,color="#51565b",lw=.7)
    a.plot(mm*axis["lengthD"],0,"s",ms=5,color="#b46a32",label="线性插值零点")
    a.set(xlim=(-.05,2.4),ylim=(-.12,.12),xlabel="距球尾的轴向距离 / mm",ylabel="下游速度 / 来流速度",title="球后轴线   ·   细步 82.10 ms")
    for index,color in [(0,"#2e5579"),(2,"#b46a32"),(5,"#547562")]:
        s=data["sections"][index];band=s["band"]
        b.plot(s["radiusD"],s["downstreamSpeedOverU"],"-o",ms=2.5,color=color,label=f"x/D = {s['xD']:g}")
        b.plot([band["r10D"],band["r90D"]],band["minSpeedOverU"]+np.array([.1,.9])*(1-band["minSpeedOverU"]),
               "s",ms=5,mfc="white",mec=color)
    b.set(xlim=(0,1.55),ylim=(-.1,1.15),xlabel="距轴线 r/D",ylabel="下游速度 / 来流速度",title="三个截面   ·   方框标出 10% 和 90%")
    for ax in (a,b):
        ax.legend(frameon=False,fontsize=8)
        ax.grid(axis="y",alpha=.18)
    fig.savefig(output,dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    evidence = HERE/"evidence";evidence.mkdir(exist_ok=True)
    actual,ui = map(Path,sys.argv[1:3])
    target=ROOT/"cases/short_sphere/wakeResolution.json"
    assert sha(target)==sha(actual/"wakeResolution.json")==read(actual/"inputs.json")["output"]
    src=read(target);work=ROOT/"outputs/v322_validation/run01"
    assert sha(work/"fine/step064/mesh.su2")==src["meshSHA256"]
    manufactured()
    proof={b:verify_branch(src[b],b,work) for b in ("coarse","fine")}
    proof["deltaLengthD"]=proof["fine"]["axisLengthD"]-proof["coarse"]["axisLengthD"]
    proof["maxWidthDifferenceD"]=max(abs(a["band"]["widthD"]-b["band"]["widthD"])
                                     for a,b in zip(src["fine"]["sections"],src["coarse"]["sections"]))
    proof["maxSameRaySpeedDifferenceOverU"]=max(max(abs(np.subtract(a["downstreamSpeedOverU"],b["downstreamSpeedOverU"])))
                                     for a,b in zip(src["fine"]["sections"],src["coarse"]["sections"]))
    proc=read(actual/"resolution_process.json");matlab=read(ui/"debugger_process.json");check=read(ui/"ui.json")
    assert proc["exitCode"]==0 and not (proc["stopReason"] or proc["livePids"])
    assert matlab["status"]=="terminal" and matlab["debuggerExitCode"]==0 and all(matlab["inputsStable"].values())
    assert not (matlab["stopReason"] or matlab["supervisorError"] or matlab["survivingProcessIDs"])
    for mode,branch in (("wakeCoarse","coarse"),("wakeFine","fine")):
        assert check[mode]["frame"]==proof[branch]["frame"] and check[mode]["sections"]==6
    assert check["steadyExcluded"]
    proof.update(manufacturedLinearCase=True,diagnosticSHA256=sha(target),seconds=proc["elapsedSeconds"],
                 maxRSSBytes=proc["peakSampledRSSBytes"],uiSeconds=matlab["elapsedSeconds"],
                 uiMaxRSSBytes=matlab["peakFamilyRSSBytes"],workbench=check)
    with (evidence/"verification.json").open("x",encoding="utf8") as f:
        json.dump(proof,f,indent=2);f.write("\n")
    plot_resolution(src["fine"],evidence/"profiles.png")
