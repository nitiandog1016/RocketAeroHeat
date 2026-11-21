"""从名义尺寸准备带引线热敏珠，使用项目内 SU2 求解同一轴向 3 m/s 冻结工况。"""
import argparse
import os
from pathlib import Path
import sys
import uuid

for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[variable]="1"
import numpy as np

from cm011Tunnel import laminar_config
from thermistorBead import AREA, CP, D, DT, ENGINE, ENGINE_SHA, K, MU, P0, PE, RE, RHO, ROOT, T0, U
from thermistorBead import read, sha, solve, write_json
from thermistorLeadMesh import prepare_mesh
from thermistorReaction import observe as read_reaction, terminal
from beadSurfaceResult import thermal_frame
from runSU2NativeDebugger import observe

PLAN=ROOT/"docs/validation/v297/PLAN.md"


def configuration():
    # 六项原残差与旧算例一致。原生弱焓保留作诊断，以强温差标量换算流体温度。
    cfg=laminar_config(RE,PE,7000)
    changes={"AXISYMMETRIC= YES":"AXISYMMETRIC= NO","MARKER_SYM= (axis)":"MARKER_SYM= (sym0, sym1)",
             "MUSCL_FLOW= YES":"MUSCL_FLOW= NO","RMS_VELOCITY-Y, RMS_ENTHALPY)":
             "RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY, RMS_SPECIES_0)",
             "RESTART_ASCII, SURFACE_PARAVIEW_ASCII, PARAVIEW_ASCII":"RESTART_ASCII, SURFACE_PARAVIEW_ASCII",
             "OUTPUT_WRT_FREQ= (1500, 500, 18000)":"OUTPUT_WRT_FREQ= (500, 500)"}
    for before,after in changes.items():
        assert cfg.count(before)==1,before
        cfg=cfg.replace(before,after)
    return cfg+f"""RESTART_SOL= YES
READ_BINARY_RESTART= NO
SOLUTION_FILENAME= initial_flow
RESTART_FILENAME= restart
KIND_SCALAR_MODEL= SPECIES_TRANSPORT
SPECIES_INIT= 0.0
DIFFUSIVITY_MODEL= CONSTANT_DIFFUSIVITY
DIFFUSIVITY_CONSTANT= {1/PE:.17g}
MARKER_WALL_SPECIES= (wall, VALUE, 1.0)
MUSCL_SPECIES= NO
CONV_NUM_METHOD_SPECIES= SCALAR_UPWIND
SPECIES_CLIPPING= NO
"""


def midbody_sources(smooth):
    plan="v299" if smooth else "v300"
    return [ROOT/f"docs/validation/{plan}/PLAN.md",ROOT/"scripts/thermistorMidbodyMesh.py"]


def build(out,first,smooth=False,midbody=False,third_grid=False):
    assert sha(ENGINE)==ENGINE_SHA
    stats=prepare_mesh(out,first,smooth,midbody,third_grid)
    with np.load(out/"nodes.npz") as data:
        xyz=data["xyz"]
        walls=np.unique(data["wallMap"])
    state=np.zeros((len(xyz),10))
    state[:,0]=np.arange(len(xyz));state[:,1:4]=xyz
    state[:,5]=-1.;state[walls,5]=0.
    state[:,8]=1-298.15;state[walls,8]=2-298.15;state[walls,9]=1.
    np.savetxt(out/"initial_flow.csv",state,fmt=["%d"]+["%.17g"]*9,delimiter=",",comments="",
               header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"')
    (out/"case.cfg").write_text(configuration(),encoding="ascii")
    write_json(out/"case.json",dict(caseKind="thermistor-B-paired-smooth" if smooth else "thermistor-B-axial-two-leads",
               meridianIntervals=stats["meridianIntervals"],sourceAreaM2=AREA,diameterM=D,
               wireDiameterM=0. if smooth else .054e-3,exposedStraightLengthM=0. if smooth else .300e-3,
               jointRadiusM=0. if smooth else .027e-3,speedMps=U,airPa=P0,airK=T0,wallK=T0+DT,rho=RHO,mu=MU,cp=CP,k=K,
               sectorDeg=30,angularCells=4,firstLayerM=first*D,outerRadiusOverD=20))
    sources=[PLAN,Path(__file__).resolve(),ROOT/"scripts/thermistorLeadMesh.py",ROOT/"scripts/thermistorBead.py",
             ROOT/"scripts/cm011Tunnel.py",out/"geometry.json",out/"nodes.npz",out/"body.json",out/"case.json"]
    if smooth:
        sources.extend((ROOT/"docs/validation/v298/PLAN.md",ROOT/"scripts/thermistorCapMesh.py"))
    if midbody or third_grid:
        sources.extend(midbody_sources(smooth))
    if third_grid:
        sources.extend((ROOT/"docs/validation/v304/PLAN.md", ROOT/"cases/thermistor_b/leadsMidbody3.json"))
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
               threads=1,maximumSeconds=720,jobToken=str(uuid.uuid4()),
               sourceFiles=[p.relative_to(ROOT).as_posix() for p in sources]))


def prepare(out,first,smooth=False,midbody=False,third_grid=False):
    out.mkdir(parents=True,exist_ok=False)
    script=Path(__file__).resolve()
    files=[script,ROOT/"scripts/thermistorLeadMesh.py",PLAN,ROOT/"third_party/gmsh/source_lock.json"]
    options=["--smooth-control"] if smooth else []
    if smooth:files.extend([ROOT/"scripts/thermistorCapMesh.py",ROOT/"docs/validation/v298/PLAN.md"])
    if midbody or third_grid:
        files.extend(midbody_sources(smooth))
        options.append("--third-midbody" if third_grid else "--refine-midbody")
    if third_grid:
        files.extend((ROOT/"docs/validation/v304/PLAN.md", ROOT/"cases/thermistor_b/leadsMidbody3.json"))
    rec=observe([sys.executable,str(script),"mesh-worker",str(out),"--first",str(first),*options],out,
                {str(f):sha(f) for f in files},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["supervisorError"],out/"debugger.log"


def probes(parent):
    thermal_frame(parent)
    read_reaction(parent,parent/"reaction0",release=True)
    read_reaction(parent,parent/"plus20",release=False)
    read_reaction(parent/"plus20",parent/"reaction20",release=True)
    # 同时以原边界装配及释放壁反力读出同一状态，加总远场、内域与壁面的各项。
    sys.path.insert(0,str(ROOT/"docs/validation/v282"))
    from observe_bead import audit
    audit(parent,parent,parent/"ledger")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("action",choices=("prepare","mesh-worker","solve","probes"))
    cli.add_argument("folder",type=Path)
    cli.add_argument("--first",type=float,choices=(.002,.001),default=.002)
    cli.add_argument("--smooth-control",action="store_true",help="以相同椭球替换两端引线，保留远区体网")
    grids=cli.add_mutually_exclusive_group()
    grids.add_argument("--refine-midbody",action="store_true",help="局部加密热敏珠中段的过渡区")
    grids.add_argument("--third-midbody",action="store_true",help="在带线热敏珠第二档中段网上细化 0.015D 至 0.1D 过渡区")
    args=cli.parse_args()
    out=args.folder.resolve();out.relative_to(ROOT/"outputs")
    if (args.smooth_control or args.refine_midbody or args.third_midbody) and args.action not in ("prepare","mesh-worker"):
        cli.error("几何选择仅用于准备输入，推进与审核读取已保存的工况")
    if args.third_midbody and args.smooth_control:
        cli.error("第三档对照用于带线热敏珠")
    if (args.refine_midbody or args.third_midbody) and args.first!=.001:
        cli.error("中段过渡区对照需要 --first .001")
    if args.action=="prepare":prepare(out,args.first,args.smooth_control,args.refine_midbody,args.third_midbody)
    elif args.action=="mesh-worker":build(out,args.first,args.smooth_control,args.refine_midbody,args.third_midbody)
    elif args.action=="probes":probes(out)
    else:
        g=read(out/"geometry.json")
        assert sha(out/"mesh.su2")==g["meshSHA256"] and sha(out/"body.json")==g["bodySHA256"]
        assert (out/"case.cfg").read_text("ascii")==configuration() and sha(ENGINE)==ENGINE_SHA
        solve(out)
        terminal(out)
