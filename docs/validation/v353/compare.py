"""独立从三档全周场、原 STL 和面表核算能量、网格误差及实际截面曲线。"""
import hashlib
from io import StringIO
import json
from pathlib import Path
import sys
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
CASE=ROOT/"cases/simple_launcher"
OUT=ROOT/"outputs/v353_validation"
sys.path.insert(0,str(ROOT/"scripts"))
from simpleLauncherTemperature import transfer
from simpleLauncherTunnel import geometry


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def section(tri,table,x):
    """与真实三角面求交，一个扇区内的两段按其在截面上的实际长度加权。"""
    yz=[]; values=[]; widths=[]
    for t,row in zip(tri,table):
        if t[:,0].min()<x<t[:,0].max():
            hits=[]
            for a,b in zip(t,np.roll(t,-1,axis=0)):
                if (a[0]-x)*(b[0]-x)<0:
                    hits.append(a[1:]+(b[1:]-a[1:])*(x-a[0])/(b[0]-a[0]))
            assert len(hits)==2
            yz.append(np.mean(hits,axis=0))
            values.append(-1000*row[3:])
            widths.append(np.linalg.norm(hits[1]-hits[0]))
    center=np.array(yz);widths=np.array(widths);values=np.array(values)
    sectors=np.floor(np.mod(np.arctan2(center[:,1],center[:,0]),2*np.pi)*48/(2*np.pi)).astype(int)
    assert np.all(np.bincount(sectors)==2)
    return np.array([np.bincount(sectors,weights=v*widths)/np.bincount(sectors,weights=widths) for v in values.T])


def save_separate_wall(table, evidence):
    """另存周向给温试验的完整原面；温界与热流同用原面面积作积分。"""
    _,tri,d,wall=geometry()
    with np.load(OUT/"outer01/fields.npz") as f:
        wet=f["xy"][f["wet"]]
        wall_theta=f["wallTheta"][1]
    length=lambda nodes:np.r_[0.,np.cumsum(np.linalg.norm(np.diff(nodes,axis=0),axis=1))]
    stations=np.interp(length(wet),length(wall),np.arange(len(wall)))
    same=abs(stations-np.rint(stations))<1e-9
    stations[same]=np.rint(stations[same])
    mapping=transfer(tri,wall,stations,d)
    centre=tri.mean(axis=1)
    sector=np.floor(np.mod(np.arctan2(centre[:,2],centre[:,1]),2*np.pi)*48/(2*np.pi)).astype(int)
    given=np.empty(len(tri))
    for j in range(48):
        which=sector==j
        given[which]=300.+.001*(mapping[which]@wall_theta[j])
    assert given.min()<300.00091 and given.max()>300.00109
    stream=StringIO()
    np.savetxt(stream,np.column_stack((table[:,:3],given,table[:,4])),delimiter=",",comments="",
        fmt=["%d"]+["%.17g"]*4,header="faceID,areaM2,pressurePa,givenWallK,heatFluxIntoWallWm2")
    raw=stream.getvalue().encode("utf8")
    meta=dict(velocityMps=[-.05,0,0],ambientTemperatureK=300.,ambientPressurePa=101325.,
        givenWall="300.001 K + 0.0001 K sin²(pi (x-x0)/(x1-x0)) cos(phi-0.37) on the cylindrical section; elsewhere 300.001 K",
        azimuth="48 sectors centred at (j+0.5)*2*pi/48 from +Y toward +Z; face-averaged prescribed wall temperature",
        pressureMethod="SU2 8.6.0 axisymmetric incompressible laminar flow",
        heatMethod="steady meridional P1 heat transport with periodic circumferential diffusion; positive into the wall",
        flowSHA256=evidence["sourceFlowSHA256"],fieldSHA256=evidence["fieldsSHA256"],
        wallSHA256=hashlib.sha256(raw).hexdigest(),areaM2=evidence["areaM2"],outwardW=evidence["pattern"]["outwardW"])
    with ZipFile(DOC/"circumferential_wall.zip","w") as z:
        for name,data in [("wall.csv",raw),("condition.json",json.dumps(meta,indent=2,ensure_ascii=False).encode("utf8"))]:
            z.writestr(ZipInfo(name),data,compress_type=ZIP_DEFLATED)


def compare():
    geom=json.loads((CASE/"geometry.json").read_bytes());stl=CASE/geom["file"]
    assert digest(stl)==geom["fileSHA256"]
    dtype=np.dtype([("n","<f4",3),("p","<f4",(3,3)),("a","<u2")])
    tri=np.fromfile(stl,offset=84,dtype=dtype)["p"].astype(float)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    evidence={};powers=[];tables=[]
    for g in ("coarse","fine","outer"):
        folder=OUT/(g+"01")
        summary=json.loads((folder/"results.json").read_bytes())
        rec=json.loads((folder/"debugger_process.json").read_bytes())
        assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
        assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
        assert rec["peakFamilyRSSBytes"]<4*1024**3 and summary["accepted"] and all(summary["gates"].values())
        assert digest(folder/"fields.npz")==summary["fieldsSHA256"] and digest(folder/"wall.csv")==summary["wallSHA256"]
        data=np.loadtxt(folder/"wall.csv",delimiter=",",skiprows=1)
        np.testing.assert_array_equal(data[:,0],np.arange(1,len(area)+1))
        np.testing.assert_allclose(data[:,1],area,rtol=0,atol=1e-15)
        power=-area@data[:,3:]
        for j,label in enumerate(("uniform","pattern","rotated")):
            assert abs(power[j]-summary[label]["outwardW"])<1e-15
            for region in geom["regions"]:
                faces=slice(region["firstFace"]-1,region["lastFace"])
                assert abs(-area[faces]@data[faces,3+j]-summary[label]["regionW"][region["name"]])<1e-15
        with np.load(folder/"fields.npz") as field:
            t=field["theta"];w=field["wallTheta"];wet=field["wet"];axis=field["axis"]
            assert t.shape[:2]==(3,48) and np.isfinite(t).all()
            np.testing.assert_allclose(t[:,:,wet],w,rtol=0,atol=1e-12)
            assert np.ptp(t[:,:,axis],axis=1).max()<1e-12
            assert t.min()>=-1e-12 and t.max()<=w.max()+1e-12
            assert np.max(abs(np.roll(t[1],12,axis=0)-t[2]))<1e-9
            assert 48*(t.shape[2]-len(axis))+len(axis)==summary["thermalDOF"]
        tables.append(data);powers.append(power)
        evidence[g]=dict(summary,observedSeconds=rec["elapsedSeconds"],peakRSSBytes=rec["peakFamilyRSSBytes"])
    grid=np.abs(powers[0]/powers[1]-1);far=np.abs(powers[1]/powers[2]-1)
    assert max(grid)<.05 and max(far)<.05
    save_separate_wall(tables[2],evidence["outer"])
    cut=.7
    q=section(tri,tables[2],cut)
    # 以原 STL 的实际柱段轴向端点重建给定温差，与解出的热流分别作图。
    group=next(r for r in geom["regions"] if r["name"]=="cylinder")
    xs=tri[group["firstFace"]-1:group["lastFace"],:,0]
    deg=(np.arange(48)+.5)*7.5
    wall=1+.1*np.cos(np.deg2rad(deg)-.37)*np.sin(np.pi*(cut-xs.min())/np.ptp(xs))**2
    np.testing.assert_allclose(np.roll(q[1],12),q[2],rtol=1e-7,atol=1e-9)
    evidence.update(gridPowerDifference=grid.tolist(),farPowerDifference=far.tolist(),sectionX=cut,
        sectorCenterDeg=deg.tolist(),wallMinusAirMilliK=wall.tolist(),uniformOutwardMilliWm2=q[0].tolist(),
        patternOutwardMilliWm2=q[1].tolist(),rotatedOutwardMilliWm2=q[2].tolist())
    (DOC/"evidence.json").write_bytes((json.dumps(evidence,indent=2,allow_nan=False)+"\n").encode())
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False,"font.size":10})
    fig,(top,bottom)=plt.subplots(2,1,figsize=(7.2,4.9),sharex=True,layout="constrained")
    top.plot(deg,wall,color="#196193",label="圆柱段周向给温");top.axhline(1,color="#667685",ls="--",label="等温壁")
    bottom.plot(deg,q[1],color="#196193");bottom.plot(deg,q[0],color="#667685",ls="--")
    top.set_title("低速风洞 · X = 0.700 m 圆柱截面");top.set_ylabel("给定壁温 − 来流温度 / mK")
    bottom.set(xlabel="周向角 / °（+Y 为 0°，+Z 为 90°）",ylabel="向外热流密度 / (mW/m²)",xlim=(0,360),xticks=np.arange(0,361,60))
    top.legend(frameon=False,loc="lower left")
    for ax in (top,bottom):
        ax.spines[["right","top"]].set_visible(False);ax.grid(axis="y",alpha=.18)
    fig.savefig(DOC/"circumference.png",dpi=160);plt.close(fig)
    print("section range",q[1].min(),q[1].max(),"grid differences",grid.tolist(),"domain",far.tolist())


if __name__=="__main__":
    compare()
