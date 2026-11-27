"""独立积分流体储热与六区功率，分别比较整步及两半步的周向扰动信号。"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from flux_transfer import DOC,ROOT,build_model,project,read
from runEmbeddedSU2 import sha
from simpleLauncherTunnel import CASE,geometry


def audit():
    out=ROOT/"outputs/v354_validation"
    heat=out/"heat02";flow=out/"flow01"
    defn,tri,d,_=geometry();air=read(CASE/"tunnel.json")
    surface=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    model=build_model(ROOT/"outputs/launcher_6639759c-6b80-41b7-bb07-55437bcd3a4a",ROOT/"outputs/v352_validation/prepared01")
    pts=model["xy"][model["tri"]];ar=np.cross(pts[:,1]-pts[:,0],pts[:,2]-pts[:,0])/2
    studied=read(out/"release01/study.json")
    tables={};cases={};transfer={}
    volume_unit=2*np.pi*air["densityKgM3"]*air["specificHeatJkgK"]*(air["wallTemperatureK"]-air["ambientTemperatureK"])*d**3
    for name,r in studied["records"].items():
        f=r["flow"];b=r["heat"];n=f["step"];old=heat/f"{f['branch']}_half_1/thermal.npz"
        if n==1:old=ROOT/"outputs/v353_validation/coarse01/fields.npz"
        new=heat/name/"thermal.npz"
        assert sha(old)==b["initialThetaSHA256"] and sha(new)==b["fieldSHA256"]
        with np.load(old) as z:before=z["theta"][0] if n==1 else z["theta"]
        with np.load(new) as z:after=z["theta"]
        change=(after-before).mean(axis=0)[model["tri"]]
        # 半径和温差在每个三角形内均为线性，三点二次精度求积独立于时间算子的集总质量。
        integral=0.
        for corner in range(3):
            weights=np.full(3,1/6);weights[corner]=2/3
            integral+=np.sum(ar*(change@weights)*(pts[:,:,1]@weights))/3
        storage=integral*volume_unit/(f["dtDU"]*d/.05)
        assert abs(storage-b["storageW"])<1e-14
        table=np.loadtxt(heat/name/"wall.csv",delimiter=",",skiprows=1)
        assert sha(heat/name/"wall.csv")==b["wallSHA256"]
        np.testing.assert_allclose(table[:,1],surface,rtol=0,atol=1e-15)
        power=float(-surface@table[:,3])
        assert abs(power-b["farW"]-storage)<1e-12
        tables[name]=table
        cases[name]=dict(storageW=storage,farW=b["farW"],outwardW=power,flowSHA256=f["frameSHA256"],
            fieldSHA256=b["fieldSHA256"],wallSHA256=b["wallSHA256"],areaM2=float(surface.sum()),
            elapsedPhysicalS=f["elapsedPhysicalS"],timeStepS=f["dtDU"]*d/.05,flowResiduals=f["last20Maximum"],
            flowSeconds=f["seconds"],heatSeconds=r["heatSeconds"],heatPeakRSSBytes=r["heatPeakRSSBytes"],
            flowPeakRSSBytes=f["peakRSSBytes"],heatResidual=b["originalResidual"],massDivergence=b["normalizedDivergence"],
            projectionFraction=b["relativeProjection"],projectionFractionOfIncrement=b["projectionRelativeToIncrement"],
            nearWallMinK=float(table[:,4].min()),nearWallMaxK=float(table[:,4].max()),
            peakOutwardFace=int(table[np.argmin(table[:,3]),0]),peakOutwardWm2=float(-table[:,3].min()),regionW=b["regionW"])
    region=next(g for g in defn["regions"] if g["name"]=="base")
    base=slice(region["firstFace"]-1,region["lastFace"])
    deg=(np.arange(48)+.5)*7.5;phi=np.deg2rad(deg)
    responses={}
    for branch,step in (("full",1),("half",2)):
        suffix=f"_{branch}_{step}"
        ctrl=tables["control"+suffix][base];pert=tables["disturbed"+suffix][base]
        # 面表保存绝对压强（Pa）与入壁热流（W/m²）；绘图为相对对照的微帕、向外微瓦每平方米。
        delta=(pert[:,2:4]-ctrl[:,2:4])*[1e6,-1e6]
        amp=2*np.mean(delta*np.exp(-1j*phi[:,None]),axis=0)
        responses[branch]=dict(pressureDifferenceMicroPa=delta[:,0].tolist(),outwardFluxDifferenceMicroWm2=delta[:,1].tolist(),
            harmonicReal=amp.real.tolist(),harmonicImag=amp.imag.tolist(),harmonicMagnitude=abs(amp).tolist())
        # 质量投影对速度增量是线性的，直接量化投影对扰动与对照之差的修正。
        states=[np.loadtxt(flow/(k+suffix)/f"restart_{step:05d}.csv",skiprows=1,delimiter=",") for k in ("control","disturbed")]
        relative=states[1].copy()
        relative[:,5:8]=model["seedVelocity"]+(states[1][:,5:8]-states[0][:,5:8])
        transfer[branch]=project(model,relative)[-1]["projectionRelativeToIncrement"]
    h1=np.array(responses["full"]["harmonicReal"])+1j*np.array(responses["full"]["harmonicImag"])
    h2=np.array(responses["half"]["harmonicReal"])+1j*np.array(responses["half"]["harmonicImag"])
    record=dict(cases=cases,comparisons=studied["comparisons"],baseAngleDeg=deg.tolist(),baseResponses=responses,
                pairedProjectionRelativeToVelocityFluxSignal=transfer,relativeHarmonicTimeStepDifference=(abs(h1-h2)/abs(h2)).tolist())
    (DOC/"evidence.json").write_bytes((json.dumps(record,indent=2,allow_nan=False)+"\n").encode())
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"font.size":10,"axes.unicode_minus":False})
    fig,axes=plt.subplots(2,1,sharex=True,figsize=(7.2,5),layout="constrained")
    for ax,key,y in zip(axes,("pressureDifferenceMicroPa","outwardFluxDifferenceMicroWm2"),("壁压差 / µPa","向外热流差 / (µW/m²)")):
        ax.plot(deg,responses["full"][key],ls="--",color="#788797",label="一个整步")
        ax.plot(deg,responses["half"][key],color="#185c91",label="两个半步")
        ax.set_ylabel(y);ax.axhline(0,color="#777",lw=.6);ax.spines[["right","top"]].set_visible(False)
    axes[0].set_title("底面扇片 · 弱扰动减无扰动 · t = 0.600 s")
    axes[0].legend(frameon=False,ncol=2,loc="upper right")
    axes[1].set(xlabel="周向角 / °（+Y 为 0°，+Z 为 90°）",xlim=(0,360),xticks=np.arange(0,361,60))
    fig.savefig(DOC/"base_response.png",dpi=160);plt.close(fig)
    print("paired projection",transfer,"harmonic timestep difference",record["relativeHarmonicTimeStepDifference"])


if __name__=="__main__":
    audit()
