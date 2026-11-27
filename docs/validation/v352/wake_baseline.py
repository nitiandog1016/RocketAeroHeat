"""从三档已验收的低速保存场读取尾迹轴线，估算直接扩展热网的存储量。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DOC = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha


def audit(job=None):
    case=ROOT/"cases/simple_launcher"
    archive=case/"wind_tunnel.zip"
    assert sha(archive)==json.loads((case/"receipt.json").read_bytes())["sha256"]
    with ZipFile(archive) as z: release=json.loads(z.read("summary.json"))
    if job is None:
        saved=json.loads((ROOT/"docs/validation/v350/evidence.json").read_bytes())["fresh"]["runDirectory"]
        job=ROOT/"outputs"/Path(saved.replace("\\","/")).name
    job=job.resolve()
    plt.rcParams.update({"font.family":["Microsoft YaHei","DejaVu Sans"],"axes.unicode_minus":False,"font.size":10})
    fig,ax=plt.subplots(figsize=(7.2,3.9),layout="constrained")
    result={}
    for grid,label,color,style in (("coarse","粗网 14D","#778899","--"),("fine","细网 14D","#b86a1e","-."),
                                   ("outer","扩域细网 22D","#15599b","-")):
        folder=job/grid; frame=folder/"restart.csv"
        info=json.loads((folder/"condition.json").read_bytes())
        rec=json.loads((folder/"engine_process.json").read_bytes())
        assert rec["status"]=="terminal" and rec["exitCode"]==0 and all(rec["inputsStable"].values())
        assert sha(frame)==rec["outputSHA256"]["restart.csv"]==release["runs"][grid]["restartSHA256"]
        assert sha(folder/"condition.json")==rec["inputSHA256"][str((folder/"condition.json").resolve())]
        wall=np.array(info["wall"]); base_d=2*wall[wall[:,0]==0,1].max()
        state=np.loadtxt(frame,skiprows=1,delimiter=",")
        axis=state[(state[:,2]==0)&(state[:,1]<=0)]
        dist=-axis[:,1]/base_d; u=axis[:,4]
        assert len(dist)==info["layers"]+1 and dist[0]==0 and u[0]==0 and np.all(np.diff(dist)>0)
        intervals=np.flatnonzero((u[:-1]>0)&(u[1:]<=0))
        assert len(intervals)==1
        i=intervals.item()
        zero=dist[i]+(dist[i+1]-dist[i])*u[i]/(u[i]-u[i+1])
        mesh=job/("heat_"+grid)/"temperature.npz"
        with np.load(mesh) as heat:
            xy=heat["xy"];tri=heat["triangles"]
            assert len(xy)==release["runs"][grid]["heatNodes"]
            # 热网三角形扫过一格周向。接轴的顶点退化，每个三棱柱需要 1、2 或 3 个四面体。
            on_axis=np.isclose(xy[:,1],0.,atol=1e-12,rtol=0)
            cells=48*int((3-on_axis[tri].sum(axis=1)).sum())
            points=48*int((~on_axis).sum())+int(on_axis.sum())
        result[grid]=dict(flowSHA256=sha(frame),heatMeshSHA256=sha(mesh),reversalLengthOverBaseDiameter=zero,
                          bracketingNodeDistances=dist[i:i+2].tolist(),peakReverseSpeedOverU=float(u.max()),
                          sweptHeatPoints=points,linearTets=cells,cooTripletBytesOnly=cells*16*24,
                          meridianNodes=len(xy),meridianTriangles=len(tri))
        selection=dist<=3.5
        ax.plot(dist[selection],u[selection],style,color=color,lw=1.9,marker=".",ms=3.5,label=label)
        ax.plot(zero,0,marker="o",mfc="white",mec=color,ms=5)
    ax.axhline(0,color="#48515b",lw=.8,zorder=0)
    ax.set(xlim=(0,3),xlabel="距平底的下游距离 / 底面直径",ylabel="有符号轴向速度 / 来流速率",
           title="0.05 m/s 风洞工况 · 稳态轴对称尾迹中心线")
    ax.text(.03,.06,"正值表示空气流向平底；圆点标出回流末端。",transform=ax.transAxes,fontsize=9)
    ax.spines[["right","top"]].set_visible(False);ax.legend(frameon=False,loc="upper right")
    fig.savefig(DOC/"wake.png",dpi=160);plt.close(fig)
    (DOC/"baseline.json").write_bytes((json.dumps(result,indent=2,allow_nan=False)+"\n").encode())
    print(json.dumps({k:v["reversalLengthOverBaseDiameter"] for k,v in result.items()}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flow-run",type=Path,help="已验收的三档低速风洞流场目录")
    audit(parser.parse_args().flow_run)
