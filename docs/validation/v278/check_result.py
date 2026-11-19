"""独立积分、从三角形坐标识别适用面带，复核三维页和完整截面 CSV。"""
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

ROOT=Path(__file__).resolve().parents[3]


def band_mean(angle,width):
    # 直接按物理弧长积分，与 MATLAB 使用的闭式多项式矩独立。
    radius=.033
    center=radius*np.deg2rad(angle)
    def value(s):
        x=s/(2*radius)
        return np.sqrt(177000)*(.9449-.5100*x**2-.5956*x**4)
    return quad(value,center-width/2,center+width/2,epsabs=1e-13)[0]/width


def check(gui,out,modules=None):
    out.mkdir(parents=True,exist_ok=False)
    stations=np.genfromtxt(gui/"stations.csv",delimiter=",",names=True,encoding="utf-8-sig")
    measured=np.genfromtxt(ROOT/"tests/data/lewis_smooth_cylinder/stations.csv",delimiter=",",names=True)
    for key in measured.dtype.names:
        np.testing.assert_allclose(stations[key],measured[key],atol=1e-10,rtol=0)
    predictions=np.array([band_mean(angle,.0051) for angle in stations["angleDeg"]])
    np.testing.assert_allclose(stations["predictedNu"],predictions,atol=1e-8,rtol=0)
    relative=predictions/measured["figureNu"]-1
    np.testing.assert_allclose(relative,stations["relativeError"],atol=1e-12,rtol=0)
    workflow=json.loads((gui/"workflow.json").read_bytes());report=workflow["report"]
    assert abs(np.sqrt(np.mean(relative**2))-report["stationRmsRelative"])<1e-12
    page=gui/"workbench"
    data=json.loads((page/"case-data.js").read_text("utf8").strip().split("=",1)[1][:-1])
    fields={f["key"]:f for f in data["fields"]}
    for k in ("wallInput","wallSolved","pressure","heatFlux"):
        assert set(fields[k]["status"])=={0} and all(v is None for v in fields[k]["values"])
    assert data["wallCoverage"]["wholeSurfaceHeatPowerW"] is None and not data["researchIdentity"]["productionQualified"]
    nv=fields["nusselt"];vertices=np.asarray(data["vertices"]);tri=vertices[data["faces"]]
    assert data["defaultField"]=="nusselt" and abs(tri[:,:,0]).max()==.076
    vectors=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    centers=tri.mean(axis=1)
    assert np.all(np.einsum("ij,ij->i",vectors[:,1:],centers[:,1:])>0)
    area=.5*np.linalg.norm(vectors,axis=1)
    np.testing.assert_allclose(area,data["faceArea"],rtol=1e-12,atol=1e-18)
    angle=np.arctan2(tri[:,:,2],tri[:,:,1])
    middle=np.rad2deg(np.unwrap(angle,axis=1).min(axis=1))+.5
    middle=(middle+180)%360-180
    valid=(np.max(abs(tri[:,:,0]),axis=1)<=.033)&(abs(middle)<=54.50000000001)
    assert np.array_equal(valid,np.asarray(nv["status"])==1) and len(area)==2160 and sum(valid)==220
    np.testing.assert_allclose([nv["values"][i] for i in np.flatnonzero(valid)],
                               [band_mean(a,.033*np.deg2rad(1)) for a in middle[valid]],rtol=0,atol=1e-8)
    assert all(nv["values"][i] is None for i in np.flatnonzero(~valid))
    assert abs(area[valid].sum()-report["nuAreaM2"])<1e-14
    assert abs(area.sum()-report["sideAreaM2"])<1e-14
    # 中跨切面读取两个三角面组成的一度面带，合格与缺项都保留在导出中。
    js="""
const fs=require('fs'),p=require('path'),[folder,modules,target]=process.argv.slice(1);
global.window={};require(p.join(folder,'case-data.js'));const m=window.ROCKET_CASE;
const s=require(p.join(modules,'section.js')),v=require(p.join(modules,'viewer.js'));
const f=m.fields.find(field=>field.key==='nusselt'),cut=s.geometry(m.vertices,m.faces,0);
if(cut.components.length!==1||!cut.components[0].closed||cut.components[0].invalidAngleSegments)throw Error('截面轮廓不闭合');
fs.writeFileSync(target,s.csv(cut,f,'单工况',v.caseCaption(m,f,0),{model:m}),{flag:'wx'});
"""
    target=out/"section.csv"
    subprocess.run(["node","-e",js,str(page.resolve()),str((modules or page).resolve()),str(target.resolve())],check=True)
    with target.open(encoding="utf-8-sig",newline="") as stream:
        rows=list(csv.DictReader(stream))
    assert len(rows)==720 and len({r["面元编号"] for r in rows})==len(rows)
    for row in rows:
        i=int(row["面元编号"])-1;val=row[nv["label"]+" / "+nv["unit"]]
        assert min(tri[i,:,0])<0<max(tri[i,:,0])
        assert row["算法"]==nv["methodNames"][nv["methodIndex"][i]]
        assert "177000" in row["工况"]
        if valid[i]:
            assert float(val)==nv["values"][i]
        else:
            assert val==""
    plot(stations,out/"station_comparison.png")
    checks=dict(points=len(stations),maxQuadratureDifferenceNu=float(max(abs(predictions-stations["predictedNu"]))),
                maxReadDifferenceNu=float(max(abs(measured["secondReadNu"]-measured["figureNu"]))),
                wholeMeshFaces=len(area),numericalNuFaces=int(sum(valid)),sectionSegments=len(rows),
                sectionComputedSegments=sum(valid[int(r["面元编号"])-1] for r in rows).item(),
                pointVsBandMaxNu=float(max(abs(stations["pointNu"]-predictions))),
                rmsRelative=float(np.sqrt(np.mean(relative**2))),maxRelative=float(max(abs(relative))))
    (out/"checks.json").write_text(json.dumps(checks,indent=2),encoding="utf8")
    return checks


def plot(s,path):
    plt.rcParams.update({"font.family":"Microsoft YaHei","axes.unicode_minus":False,"font.size":10})
    fig,(top,bottom)=plt.subplots(2,1,sharex=True,figsize=(7.4,5.8),gridspec_kw={"height_ratios":[2,1]})
    x=s["angleDeg"]
    top.plot(x,s["predictedNu"],"o-",color="#224e76",lw=1.5,label="Frössling · 铜带平均")
    top.errorbar(x,s["figureNu"],yerr=[s["figureNu"]-s["readLowNu"],s["readHighNu"]-s["figureNu"]],
                 fmt="D",color="#99532b",ms=5,capsize=3,label="图 9 实验读数 · 竖线为读图范围")
    top.set_ylabel("局部努塞尔数 Nu_D");top.legend(frameon=False,loc="lower left")
    bottom.axhline(0,color="#707c80",lw=.8)
    bottom.plot(x,100*s["relativeError"],"o-",color="#224e76")
    bottom.set_ylabel("预测相对差 / %");bottom.set_xlabel("铜带中心距迎流驻点的角度 / °")
    bottom.set_xticks(x)
    for ax in (top,bottom):
        ax.grid(axis="y",color="#e3e7e7")
        ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    fig.suptitle("NASA Lewis 光滑圆柱   ·   低湍流空气   ·   Re_D = 177 000",fontsize=13)
    fig.tight_layout(rect=(.015,.015,.99,.95))
    fig.savefig(path,dpi=155);plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gui",type=Path)
    parser.add_argument("new_output",type=Path)
    parser.add_argument("--web-modules",type=Path)
    args=parser.parse_args()
    print(check(args.gui,args.new_output,args.web_modules))
