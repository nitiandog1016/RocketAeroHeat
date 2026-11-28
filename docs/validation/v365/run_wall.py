"""在固定 M0.30 民用试件上运行两组等温壁、轴对称层流近壁网格。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

DOC = Path(__file__).resolve().parent
ROOT = DOC.parents[2]
BASE = ROOT/"outputs/v365_validation"
sys.path[:0] = [str(ROOT/"scripts"), str(DOC.parent/"v359")]
from simpleLauncherTunnel import CASE, geometry
from cm011Tunnel import volume, write_volume_mesh
from plate_check import EXE, EXE_SHA, read, save, sha
from runSU2NativeDebugger import observe

GRIDS = {"normal4": (72,4e-6), "normal2": (88,2e-6)}
FLIGHT = DOC.parent/"v358/evidence/cases.json"


def build_grid(grid):
    """保留原壁节点，首层沿法向展开，再平滑接到原同心远界。"""
    definition, _, diameter, meridian = geometry()
    layers, spacing = GRIDS[grid]
    wall = meridian*diameter
    _, cells, markers = volume(meridian,layers,14,1)
    center = np.array([definition["lengthM"]/2,0.])
    radial = wall-center
    radii = np.linalg.norm(radial,axis=1)
    radial /= radii[:,None]
    differences = np.diff(wall,axis=0)
    unit = differences/np.linalg.norm(differences,axis=1)[:,None]
    tangent = np.vstack((unit[0],unit[1:]+unit[:-1],unit[-1]))
    normals = np.column_stack((-tangent[:,1],tangent[:,0]))
    normals /= np.linalg.norm(normals,axis=1)[:,None]
    normals[[0,-1]] = [[-1.,0.],[1.,0.]]
    points = np.empty((layers+1,len(wall),2))
    for i,(point,ray,normal,radius) in enumerate(zip(wall,radial,normals,radii)):
        length = 14*diameter-radius
        powers = np.arange(layers)
        ratio = brentq(lambda g:spacing*np.sum(g**powers)-length,1,2)
        distance = np.r_[0.,np.cumsum(spacing*ratio**powers)]
        distance[-1] = length
        blend = .004/(distance+.004)*(1-distance/length)
        points[:,i] = point+distance[:,None]*(ray+blend[:,None]*(normal-ray))
    assert (points[:,:,1]>=0).all() and (points[:,[0,-1],1]==0).all()
    np.testing.assert_array_equal(points[0],wall)
    np.testing.assert_allclose(np.linalg.norm(points[-1]-center,axis=1),14*diameter,rtol=1e-12)
    # 将求解第一轴定为顺流方向；反射坐标同时调整单元和边界绕序。
    points[:,:,0] = definition["lengthM"]-points[:,:,0]
    xy = points.reshape(-1,2)
    cells = cells[:,::-1]
    markers = {key:pairs[:,::-1] for key,pairs in markers.items()}
    corners = xy[cells]
    sides = np.roll(corners,-1,axis=1)-corners
    following = np.roll(sides,-1,axis=1)
    if np.min(sides[:,:,0]*following[:,:,1]-sides[:,:,1]*following[:,:,0])<=0:
        raise ValueError("单元翻转，拒绝启动")
    return definition,xy,cells,markers,points.shape[:2],wall


def prepare(grid,folder):
    """登记几何与来流，写出原生状态与归一化选项。"""
    model,xy,cells,markers,shape,wall = build_grid(grid)
    archive=read(FLIGHT)
    condition,=[case for case in archive["cases"] if case["mach"]==.30]
    assert model["fileSHA256"]==archive["surfaceSHA256"] and condition["wallGivenK"]==300 and condition["altitudeM"]==5000
    speed=condition["speedMps"]
    temperature=condition["freeTemperatureK"]
    pressure=condition["freePressurePa"]
    assert abs(speed-.3*np.sqrt(1.4*287.05287*temperature))<1e-7
    folder.mkdir(parents=True,exist_ok=False)
    write_volume_mesh(folder/"mesh.su2",xy,cells,markers)
    save(folder/"condition.json",dict(grid=grid,shape=list(shape),meridianM=wall.tolist(),lengthM=model["lengthM"],
        surfaceSHA256=model["fileSHA256"],speedMps=speed,freeTemperatureK=temperature,freePressurePa=pressure,
        densityKgM3=pressure/(287.05287*temperature),givenWallK=300,transitionModel="laminar throughout",
        initial="native uniform freestream, no restart",firstSpacingM=GRIDS[grid][1]))
    fields="RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_ENERGY"
    cfg=f"""SOLVER= NAVIER_STOKES
KIND_TURB_MODEL= NONE
AXISYMMETRIC= YES
MATH_PROBLEM= DIRECT
RESTART_SOL= NO
SYSTEM_MEASUREMENTS= SI
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.05287
INIT_OPTION= TD_CONDITIONS
FREESTREAM_OPTION= TEMPERATURE_FS
FREESTREAM_TEMPERATURE= {temperature:.17g}
FREESTREAM_PRESSURE= {pressure:.17g}
REF_DIMENSIONALIZATION= FREESTREAM_PRESS_EQ_ONE
MACH_NUMBER= 0.30
AOA= 0
SIDESLIP_ANGLE= 0
REF_LENGTH= {model['lengthM']:.17g}
REF_AREA= {np.pi*(model['maxDiameterM']/2)**2:.17g}
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716e-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
PRANDTL_LAM= 0.72
MARKER_ISOTHERMAL= (wall,300)
MARKER_FAR= (far)
MARKER_SYM= (axis)
MARKER_MONITORING= (wall)
MARKER_PLOTTING= (wall)
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
TIME_DISCRE_FLOW= EULER_IMPLICIT
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
RAMP_MUSCL= YES
KIND_MUSCL_RAMP= ITERATION
RAMP_MUSCL_COEFF= (0,1,200)
MGLEVEL= 0
CFL_NUMBER= 3
CFL_ADAPT= YES
CFL_ADAPT_PARAM= (0.5,1.3,1,100,0.01,50)
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 0.001
LINEAR_SOLVER_ITER= 20
ITER= 2500
CONV_FIELD= ({fields})
CONV_RESIDUAL_MINVAL= -10
CONV_STARTITER= 500
SCREEN_OUTPUT= (INNER_ITER, {fields}, TOTAL_HEATFLUX)
SCREEN_WRT_FREQ_INNER= 50
HISTORY_OUTPUT= (ITER,RMS_RES,HEAT)
TABULAR_FORMAT= CSV
MESH_FILENAME= mesh.su2
OUTPUT_FILES= (RESTART_ASCII,SURFACE_PARAVIEW_ASCII)
OUTPUT_WRT_FREQ= (250,250)
WRT_RESTART_COMPACT= NO
WRT_RESTART_OVERWRITE= NO
WRT_SURFACE_OVERWRITE= NO
VOLUME_OUTPUT= (COORDINATES,SOLUTION,PRIMITIVE)
"""
    (folder/"case.cfg").write_text(cfg,encoding="ascii")
    # 体初场另存准确的无量纲守恒量。原生迭代 0 帧是第一次更新后的场。
    u=speed/np.sqrt(287.05287*temperature)
    np.savez_compressed(folder/"initial.npz",xy=xy,conservative=np.tile([1.,u,0.,1/.4+.5*u*u],(len(xy),1)))


def execute(folder):
    assert sha(EXE)==EXE_SHA
    sources=[Path(__file__),DOC/"PLAN.md",FLIGHT,CASE/"geometry.json",CASE/"launcher-m.stl",
             ROOT/"scripts/simpleLauncherTunnel.py",ROOT/"scripts/cm011Tunnel.py",ROOT/"scripts/runSU2NativeDebugger.py",
             folder/"case.cfg",folder/"condition.json",folder/"initial.npz",folder/"mesh.su2",EXE]
    receipt=observe([str(EXE),"case.cfg","--threads","2"],folder,{str(p.resolve()):sha(p) for p in sources},300,160*1024**2)
    print(json.dumps({k:receipt[k] for k in ("elapsedSeconds","debuggerExitCode","stopReason","survivingProcessIDs")}),flush=True)
    if receipt["debuggerExitCode"] or receipt["stopReason"] or "Exit Success (SU2_CFD)" not in (folder/"debugger.log").read_text(errors="replace"):
        raise RuntimeError(f"检查原生终态及已保存帧 {folder}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid",choices=GRIDS)
    parser.add_argument("folder",type=Path)
    args=parser.parse_args()
    folder=args.folder.resolve()
    folder.relative_to(BASE)
    prepare(args.grid,folder)
    execute(folder)
