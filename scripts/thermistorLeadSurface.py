"""核验带引线热敏珠的原生终点和壁反力，按完整暴露面守恒分配并分别积分。"""
import argparse
from pathlib import Path
import os
import zipfile

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[variable]="1"
import numpy as np

from beadSurfaceResult import state_and_reaction
from thermistorBead import AREA,D,DT,P0,RHO,ROOT,T0,U,read,sha,write_json
from thermistorLeadMesh import original
from runThermistorLeadCase import configuration
from thermistorReaction import terminal


def collect(parent,out):
    record=terminal(parent)
    g=read(parent/"geometry.json")
    c=read(parent/"case.json")
    smooth=c["caseKind"]=="thermistor-B-paired-smooth"
    assert smooth or c["caseKind"]=="thermistor-B-axial-two-leads"
    assert c["wireDiameterM"]==(0. if smooth else .054e-3)
    assert c["wallK"]==T0+DT and c["speedMps"]==U and c["sectorDeg"]==30
    assert (parent/"case.cfg").read_text("ascii")==configuration()
    body,links,area=original(smooth)
    assert body==read(parent/"body.json") and sha(parent/"body.json")==g["bodySHA256"]
    for name in ("case.json","geometry.json","body.json","nodes.npz"):
        assert sha(parent/name)==record["inputSHA256"][str(parent/name)]
    assert sha(parent/"mesh.su2")==g["meshSHA256"] and sha(parent/"nodes.npz")==g["nodeMapSHA256"]
    with np.load(parent/"nodes.npz") as data:
        wall,near=data["wallMap"],data["firstMap"]
        xyz,ids=data["xyz"],data["ids"]
        assert np.array_equal(area,data["area"])
        assert np.array_equal(wall,ids[links[:,:,0],links[:,:,1]])
        assert np.array_equal(near,ids[links[:,:,0],links[:,:,1]+c["meridianIntervals"]+1])
    nodes=np.unique(wall)
    # 核对两个不同收敛端点与各自的独立壁反力读数，每个读数只使用自己的输入态。
    base_cfg=configuration()
    assert (parent/"plus20/case.cfg").read_text("ascii")==base_cfg.replace("ITER= 7000\n","ITER= 20\n")
    assert sha(parent/"plus20/initial_flow.csv")==sha(parent/"restart.csv")
    states=[]
    for source,probe in ((parent,parent/"reaction0"),(parent/"plus20",parent/"reaction20")):
        cfg=(source/"case.cfg").read_text("ascii")
        suffix="WRT_RESTART_COMPACT= NO\n"
        expected=cfg.replace("ITER= 7000\n","ITER= 1\n").replace("\nITER= 20\n","\nITER= 1\n")
        expected=expected.replace("(wall, VALUE, 1.0)","(wall, FLUX, 0.0)").replace(
            "VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)","VOLUME_OUTPUT= (COORDINATES, SOLUTION, RESIDUAL)")
        assert (probe/"case.cfg").read_text("ascii")==expected+suffix
        states.append(state_and_reaction(source,probe,nodes,xyz))
    state,reaction,res,rec=states[0]
    angle=np.pi/6
    normal=max(np.max(abs(state[ids[0],7])),np.max(abs(state[ids[-1],7]*np.cos(angle)-state[ids[-1],6]*np.sin(angle))))
    assert normal<1e-8
    release=-reaction[nodes].sum()
    change=abs(-states[1][1][nodes].sum()/release-1)
    e=read(parent/"ledger/energy.json")
    assert e["inputFrameSHA256"]==sha(parent/"restart.csv") and e["sampledThetaExact"]
    assert e["maxWallNodeReactionDifferenceW"]<1e-11 and abs(e["relativeDifference"])<.005 and change<.002
    assert abs(e["wallReleaseW"]-release)<1e-10
    assert sha(parent/"ledger/thermistor_observe.exe")==e["observerSHA256"]
    # 完整面上的三角顶点积分权重等于面积的 1/3，反射边与极点按真实节点累计。
    lump=np.bincount(wall.ravel(),weights=np.repeat(area/3,3),minlength=len(xyz))
    q=np.mean(reaction[wall]/lump[wall],axis=1)
    gauge=state[wall,4].mean(axis=1)*RHO*U**2
    temperature=T0+DT*state[near,9].mean(axis=1)
    watts=-area*q
    assert abs(watts.sum()/release-1)<1e-12 and np.isfinite(q).all()
    groups=np.array(body["regionID"])
    rows=[]
    for i,name in enumerate(body["regionNames"],1):
        f=np.flatnonzero(groups==i)
        peak=int(f[np.argmax(-q[f])])
        rows.append(dict(regionID=i,name=name,faces=len(f),areaM2=float(area[f].sum()),outwardW=float(watts[f].sum()),
                         pressurePa=[float(P0+gauge[f].min()),float(P0+gauge[f].max())],
                         intoWallWm2=[float(q[f].min()),float(q[f].max())],peakFaceID=peak+1))
    if smooth:
        bead=watts.sum()
        tip=wires=0.
    else:
        bead=watts[(groups>=3)&(groups<=5)].sum()
        tip=watts[(groups==1)|(groups==7)].sum()
        wires=watts[(groups==2)|(groups==6)].sum()
    assert abs(bead+wires+tip-release)<1e-12
    exp=read(ROOT/"docs/validation/v282/evidence/comparison.json")["experiment"]
    out.mkdir(exist_ok=False)
    table=out/"wall.csv"
    np.savetxt(table,np.column_stack((np.arange(len(area)),groups,area,P0+gauge,gauge,q,temperature)),
               fmt=["%d","%d"]+["%.17g"]*5,delimiter=",",comments="",
               header="faceID,regionID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,nearWallCellTemperatureK")
    summary=dict(areaM2=float(area.sum()),faces=len(q),regions=rows,beadAndJointsW=float(bead),
                 wireSidesW=float(wires),wireEndsW=float(tip),totalOutwardW=float(release),
                 beadMeanHWm2K=float(bead/(AREA*DT)),experimentMeanHWm2K=exp["hWm2K"],
                 jointExperimentVerifiedM2=0,peakFaceID=int(np.argmin(q))+1,relativeChange20=float(change),
                 residualWindow=res,plus20Window=states[1][2],symmetryNormalMaxOverU=float(normal),energy=e,
                 meshSHA256=g["meshSHA256"],bodySHA256=g["bodySHA256"],wallSHA256=sha(table),
                 caseSHA256=sha(parent/"case.json"),sourceFrameSHA256=sha(parent/"restart.csv"),
                 sourceReactionSHA256=sha(parent/"reaction0/restart.csv"),nativeSeconds=rec["elapsedSeconds"])
    write_json(out/"summary.json",summary)
    with zipfile.ZipFile(out/"surface.zip","x") as z:
        for p in (parent/"case.json",parent/"body.json",table,out/"summary.json"):
            item=zipfile.ZipInfo(p.name,date_time=(2026,1,1,0,0,0))
            z.writestr(item,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED)
    return summary


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("parent",type=Path)
    args=cli.parse_args()
    folder=args.parent.resolve();folder.relative_to(ROOT/"outputs")
    collect(folder,folder/"surface")
