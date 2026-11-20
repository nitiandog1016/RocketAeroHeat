"""核对独立强标量计算、同阶段远场收支及重建到每一张原面的压力、热流和温度。"""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import CP,D,DT,P0,RHO,T0,U,read,sha
from thermistorReaction import terminal
from beadSurfaceResult import projection,thermal_frame


def reaction(parent,probe):
    terminal(probe)
    assert sha(parent/"restart.csv")==sha(probe/"initial_flow.csv")
    assert sha(parent/"mesh.su2")==sha(probe/"mesh.su2")
    cfg=(parent/"case.cfg").read_text("ascii")
    cfg,n=re.subn(r"^MARKER_WALL_SPECIES= \(wall, VALUE, (?:0\.0|0\.5)\)$","MARKER_WALL_SPECIES= (wall, FLUX, 0.0)",cfg,flags=re.M)
    assert n==1
    cfg=re.sub(r"^ITER= \d+$","ITER= 1",cfg,flags=re.M).replace("COORDINATES, SOLUTION, PRIMITIVE","COORDINATES, SOLUTION, RESIDUAL")
    assert (probe/"case.cfg").read_text("ascii")==cfg+"WRT_RESTART_COMPACT= NO\n"
    with (probe/"restart.csv").open() as f:
        cols=next(csv.reader(f))
    raw=np.loadtxt(probe/"restart.csv",delimiter=",",skiprows=1,usecols=(0,cols.index("Residual_Species_0")))
    assert np.array_equal(raw[:,0],np.arange(len(raw))) and np.isfinite(raw).all()
    return raw[:,1]


def audit(folder,scale):
    reference=ROOT/"outputs/v282_validation/strong_outer01"
    strong=folder/"strong"
    stages,residuals=thermal_frame(strong)
    plus,res20=thermal_frame(folder/"plus20")
    assert sha(folder/"plus20/initial_flow.csv")==sha(strong/"restart.csv")
    assert sha(reference/"mesh.su2")==sha(strong/"mesh.su2")==sha(folder/"plus20/mesh.su2")
    cfg=(reference/"case.cfg").read_text("ascii")
    line="MARKER_WALL_SPECIES= (wall, VALUE, 1.0)"
    strong_cfg=cfg.replace(line,f"MARKER_WALL_SPECIES= (wall, VALUE, {scale})")
    assert (strong/"case.cfg").read_text("ascii")==strong_cfg
    assert (folder/"plus20/case.cfg").read_text("ascii")==re.sub(r"^ITER= \d+$","ITER= 20",strong_cfg,flags=re.M)
    initial=np.loadtxt(strong/"initial_flow.csv",delimiter=",",skiprows=1)
    original_initial=np.loadtxt(reference/"initial_flow.csv",delimiter=",",skiprows=1)
    assert np.array_equal(initial[:,:-1],original_initial[:,:-1]) and np.array_equal(initial[:,-1],scale*original_initial[:,-1])
    with np.load(reference/"nodes.npz") as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    body=read(reference/"body.json");v=np.asarray(body["V"]);f=np.asarray(body["F"])-1
    area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)*.5
    assign,wall=projection(v,f,ids[:,:65])
    loaded=np.loadtxt(strong/"restart.csv",delimiter=",",skiprows=1)
    theta=loaded[:,-1]
    np.testing.assert_allclose(loaded[:,1:4],xyz,rtol=0,atol=1e-13)
    assert np.array_equal(loaded[:,0],np.arange(len(xyz))) and np.isfinite(loaded).all()
    assert theta.min()>=0 and theta.max()<=scale and np.all(theta[wall]==scale)
    assert np.all(loaded[wall,5:8]==0) and np.min(P0+RHO*U**2*loaded[:,4])>0
    factor=12*D**2*RHO*U*CP*DT
    rhs=reaction(strong,folder/"reaction0")
    rhs20=reaction(folder/"plus20",folder/"reaction20")
    dual=assign.T@area
    q=assign@(rhs[wall]*factor/dual)
    power=area@q
    near_ids=np.zeros(len(xyz),int);near_ids[ids[:,:65]]=ids[:,65:130]
    near=T0+DT*(assign@theta[near_ids[wall]])
    p=P0+RHO*U**2*(assign@loaded[wall,4])
    original=ROOT/"outputs/v282_validation/result_outer01/wall.csv"
    assert sha(original)=="08fea8e578f50c45168b2c221e9e4d8150e64b8bbafa871d4aa80e8fe9640bfe"
    ref=np.genfromtxt(original,delimiter=",",names=True)
    assert np.array_equal(ref["faceID"],np.arange(len(f))) and np.array_equal(ref["areaM2"],area)
    expected_q=ref["heatFluxWm2"]*scale
    qerr=float(max(abs(q-expected_q)))
    terr=float(max(abs(near-(T0+scale*(ref["nearWallCellTemperatureK"]-T0)))))
    perr=float(max(abs(p-ref["pressurePa"])))
    assert terr<1e-3 and perr<1e-3
    # 未修改强边界的另一读法在同一输入态上记录壁面与远场项。按原生全局编号排序。
    frames=[]
    for name in ("control","observed"):
        step=folder/name;obs=read(step/"step_process.json")
        assert obs["exitCode"]==0 and not obs["livePids"] and obs["stopReason"] is None
        assert sha(step/"initial_flow.csv")==sha(strong/"restart.csv") and sha(step/"mesh.su2")==sha(strong/"mesh.su2")
        assert (step/"case.cfg").read_text("ascii")==re.sub(r"^ITER= \d+$","ITER= 1",strong_cfg,flags=re.M)
        frames.append(np.loadtxt(step/"restart.csv",delimiter=",",skiprows=1))
    assert sha(folder/"control/case.cfg")==sha(folder/"observed/case.cfg")
    assert np.max(abs(frames[0]-frames[1]))<2e-7
    sample=np.loadtxt(folder/"observed/original_scalar.csv",delimiter=",",skiprows=1)
    sample=sample[np.argsort(sample[:,0])]
    assert np.array_equal(sample[:,0],loaded[:,0]) and np.array_equal(sample[:,1],theta)
    wall_diff=float(np.max(abs(rhs[wall]+sample[wall,2]))*factor)
    assert wall_diff<1e-11
    with (folder/"observed/scalar_stages.csv").open() as stream:
        stage={k:float(x) for k,x in csv.reader(stream)}
    wall_W=(stage["beforeWall"]-stage["afterWall"])*factor
    far_W=(stage["afterFar"]-stage["beforeFar"])*factor
    assert abs(power+wall_W)<1e-11 and abs(rhs.sum()*factor+far_W)<1e-11
    power20=rhs20[wall].sum()*factor
    if scale==0:
        assert np.all(q==0) and np.all(near==T0) and np.all(rhs==0) and np.all(rhs20==0)
        assert power==wall_W==far_W==power20==0
        relq=balance=relative_power=None
    else:
        relq=float(max(abs((q-expected_q)/expected_q)))
        relative_power=float(abs(power/(area@expected_q)-1))
        balance=float(abs(1-far_W/wall_W))
        assert max(relq,relative_power)<.0005 and abs(power20/power-1)<.002 and balance<.005
    rows=np.column_stack((np.arange(len(f)),p,q,near))
    path=folder/"direct_wall.csv"
    if path.exists():
        np.testing.assert_array_equal(rows,np.loadtxt(path,delimiter=",",skiprows=1))
    else:
        with path.open("x",encoding="ascii") as stream:
            np.savetxt(stream,rows,delimiter=",",comments="",fmt=["%d","%.17g","%.17g","%.17g"],
                       header="faceID,pressurePa,heatFluxWm2,nearWallCellTemperatureK")
    return dict(wallTemperatureK=T0+DT*scale,faces=len(f),areaM2=area.sum(),wallPowerIntoW=power,
                continuedPowerIntoW=power20,wallReleaseW=wall_W,farOutW=far_W,relativeBoundaryDifference=balance,
                maxFluxErrorWm2=qerr,maxRelativeFluxError=relq,relativePowerError=relative_power,
                maxFirstNodeErrorK=terr,maxPressureErrorPa=perr,temperatureRangeK=[T0+DT*min(theta),T0+DT*max(theta)],
                maxObservedReactionErrorW=wall_diff,sixResiduals=residuals,continuedResiduals=res20,
                inputSHA256=sha(strong/"initial_flow.csv"),meshSHA256=sha(strong/"mesh.su2"),
                configSHA256=sha(strong/"case.cfg"),frameSHA256=sha(strong/"restart.csv"),
                tableSHA256=sha(folder/"direct_wall.csv"),strongSeconds=stages["elapsedSeconds"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();report={name:audit(args.directory.resolve()/name,scale) for name,scale in (("half",.5),("zero",0.))}
    with args.output.open("x",encoding="utf8") as stream:
        json.dump(report,stream,indent=2,allow_nan=False);stream.write("\n")
