"""66 mm 加热圆柱的全风洞三维非定常入口试算；保留全部安装壁和原生终态。"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import psutil
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha

EXE=ROOT/"tools/runtime/su2/8.6.0/bin/SU2_CFD.exe"
EXE_SHA="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
BASE=ROOT/"outputs/v279_validation"
NPHI,NR,NX=72,40,24
D,W,H=.066,.152,.686
H1=8e-6


def save(path,value):
    with path.open("x",encoding="utf8") as stream:
        json.dump(value,stream,indent=2,allow_nan=False)


def stretch(count,length):
    """两端壁面向中间光滑放大，使用米制首间隔控制双曲正切映射。"""
    t=np.linspace(0,1,count+1)
    def positions(b):
        return .5*(1+np.tanh(b*(2*t-1))/np.tanh(b))
    b=brentq(lambda b:length*positions(b)[1]-H1,.01,12)
    return positions(b)


def mesh(folder):
    phi=2*np.pi*np.arange(NPHI)/NPHI
    direction=np.column_stack((np.cos(phi),np.sin(phi)))
    wall=D/2*direction
    square=direction/np.max(abs(direction),axis=1)[:,None]
    outer=square*np.column_stack((np.where(square[:,0]>=0,5*D,12*D),np.full(NPHI,H/2)))
    s=stretch(NR,(H-D)/2)[:,None,None]
    # 第一层沿圆柱法线；外部逐步变为矩形试验段，角点与四条外边界重合。
    yz=wall+s*((H-D)/2)*direction+(s*s*(3-2*s))*(outer-(H/2)*direction)
    yz[1]=wall+H1*direction
    x=W*(stretch(NX,W)-.5)
    shape=(NX+1,NR+1,NPHI)
    xyz=np.empty(shape+(3,));xyz[...,0]=x[:,None,None];xyz[...,1:]=yz
    xyz=xyz.reshape(-1,3)
    ids=np.arange(len(xyz)).reshape(shape)
    a=ids[:-1,:-1];b=ids[:-1,1:];e=ids[1:,:-1];f=ids[1:,1:]
    roll=lambda v:np.roll(v,-1,axis=-1)
    cells=np.stack((a,b,roll(b),roll(a),e,f,roll(f),roll(e)),axis=-1).reshape(-1,8)
    pts=xyz[cells];minimum=float("inf")
    # 三个参数方向的有向边在全部八角点检查；六面体连接由同一结构化编号建立。
    for ax in range(8):
        rx=(1,0,3,2,5,4,7,6)[ax]
        ry=(3,2,1,0,7,6,5,4)[ax]
        rz=(4,5,6,7,0,1,2,3)[ax]
        er=(pts[:,rx]-pts[:,ax])*(1 if ax in (0,3,4,7) else -1)
        et=(pts[:,ry]-pts[:,ax])*(1 if ax in (0,1,4,5) else -1)
        ex=(pts[:,rz]-pts[:,ax])*(1 if ax<4 else -1)
        minimum=min(minimum,float(np.einsum("ij,ij->i",np.cross(er,et),ex).min()))
    assert minimum>0
    block=cells.reshape(NX,NR,NPHI,8)
    skin=block[:,0,:][:,:,[0,4,7,3]].reshape(-1,4)
    ends=np.concatenate((block[0][:,:,[0,3,2,1]].reshape(-1,4),block[-1][:,:,[4,5,6,7]].reshape(-1,4)))
    outside=block[:,-1,:][:,:,[1,2,6,5]]
    mid=phi+np.pi/NPHI
    inlet=np.cos(mid)>abs(np.sin(mid));outlet=-np.cos(mid)>abs(np.sin(mid))
    markers={"cylinder":skin,"mounts":ends,"inlet":outside[:,inlet].reshape(-1,4),
             "outlet":outside[:,outlet].reshape(-1,4),"ceiling_floor":outside[:,~(inlet|outlet)].reshape(-1,4)}
    # 独立按每个单元六面统计邻接，要求唯一外边界与具名标记完全相等。
    subfaces=((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    unique,counts=np.unique(np.sort(np.concatenate([cells[:,f] for f in subfaces]),axis=1),axis=0,return_counts=True)
    all_outer=np.sort(np.concatenate(list(markers.values())),axis=1)
    assert max(counts)==2 and len(np.unique(all_outer,axis=0))==len(all_outer)
    assert np.array_equal(unique[counts==1],np.unique(all_outer,axis=0))
    with (folder/"mesh.su2").open("x",encoding="ascii") as stream:
        stream.write("NDIME= 3\nNELEM= %d\n"%len(cells))
        np.savetxt(stream,cells,fmt="12 "+" ".join(["%d"]*8))
        stream.write("NPOIN= %d\n"%len(xyz));np.savetxt(stream,xyz,fmt="%.17g")
        stream.write("NMARK= %d\n"%len(markers))
        for name,faces in markers.items():
            stream.write("MARKER_TAG= %s\nMARKER_ELEMS= %d\n"%(name,len(faces)))
            np.savetxt(stream,faces,fmt="9 %d %d %d %d")
    np.savez_compressed(folder/"wall.npz",xyz=xyz[ids[:,0]],nodeIDs=ids[:,0],x=x,phi=phi)
    return xyz,dict(nodes=len(xyz),cells=len(cells),minimumCornerJacobianM3=minimum,
                    markers={k:len(v) for k,v in markers.items()},spanFirstM=x[1]-x[0],
                    radialFirstM=float(np.linalg.norm(yz[1]-yz[0],axis=1).min()))


def initial(folder,xyz,c):
    """保存完整有量纲冷启动场。用壁面阻尼的圆柱势流作初猜，尾迹加弱跨流扰动。"""
    u=c["velocityMps"];rho=c["densityKgM3"];mu=c["dynamicViscosityPaS"]
    x,y,z=xyz.T;r2=y*y+z*z;r=np.sqrt(r2);gap=abs(r-D/2)
    side=W/2-abs(x);vertical=H/2-abs(z)
    damp=(1-np.exp(-gap/.00035))*(1-np.exp(-abs(side)/.002))*(1-np.exp(-abs(vertical)/.002))
    uy=-u*(1-(D/2)**2*(y*y-z*z)/r2**2)*damp
    uz=u*(2*(D/2)**2*y*z/r2**2+.001*np.exp(-((y+2*D)/D)**2-(z/D)**2))*damp
    solid=(gap<1e-12)|(abs(side)<1e-12)|(abs(vertical)<1e-12)
    uy[solid]=0;uz[solid]=0
    temp=c["airTemperatureK"]+(c["wallTemperatureK"]-c["airTemperatureK"])*np.exp(-gap/.00035)
    zero=np.zeros(len(xyz));one=np.ones(len(xyz))
    k=1.5*(.005*u)**2;omega=rho*k/(5*mu)
    distance=np.minimum(gap,np.minimum(abs(side),abs(vertical)))
    w=omega+60*mu/(rho*.075*(distance+H1)**2)
    turb=k*damp;turb[solid]=0
    intensity=.5;retheta=1173.51-589.428*intensity+.2196/intensity**2
    fields=np.column_stack((np.arange(len(xyz)),xyz,zero,zero,uy,uz,
        c["cpJkgK"]*(temp-298.15),turb,w,one,retheta*one,zero,one))
    names="PointID x y z Pressure Velocity_x Velocity_y Velocity_z Enthalpy Turb_Kin_Energy Omega LM_gamma LM_Re_t LM_gamma_sep LM_gamma_eff".split()
    np.savetxt(folder/"seed_00000.csv",fields,delimiter=",",fmt=["%d"]+["%.17g"]*(len(names)-1),
               comments="",header=",".join('"'+n+'"' for n in names))


def configuration(c):
    u=c["velocityMps"];temp=c["airTemperatureK"];tw=c["wallTemperatureK"]
    return f"""SOLVER= INC_RANS
KIND_TURB_MODEL= SST
KIND_TRANS_MODEL= LM
INC_ENERGY_EQUATION= YES
INC_DENSITY_MODEL= CONSTANT
FLUID_MODEL= CONSTANT_DENSITY
INC_NONDIM= DIMENSIONAL
INC_DENSITY_INIT= {c['densityKgM3']:.17g}
INC_VELOCITY_INIT= (0, {-u:.17g}, 0)
INC_TEMPERATURE_INIT= {temp:.17g}
SPECIFIC_HEAT_CP= {c['cpJkgK']:.17g}
VISCOSITY_MODEL= CONSTANT_VISCOSITY
MU_CONSTANT= {c['dynamicViscosityPaS']:.17g}
CONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY
THERMAL_CONDUCTIVITY_CONSTANT= {c['conductivityWmK']:.17g}
PRANDTL_TURB= 0.90
FREESTREAM_TURBULENCEINTENSITY= 0.005
FREESTREAM_TURB2LAMVISCRATIO= 5
MARKER_INLET_TURBULENT= (inlet, 0.005, 5)
INC_INLET_TYPE= VELOCITY_INLET
INC_OUTLET_TYPE= PRESSURE_OUTLET
MARKER_INLET= (inlet, {temp:.17g}, {u:.17g}, 0, -1, 0)
MARKER_OUTLET= (outlet, 0)
MARKER_ISOTHERMAL= (cylinder, {tw:.17g})
MARKER_HEATFLUX= (mounts, 0, ceiling_floor, 0)
MARKER_PLOTTING= (cylinder)
MARKER_MONITORING= (cylinder)
TIME_DOMAIN= YES
TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER
TIME_STEP= {.05*D/u:.17g}
MAX_TIME= 1.0
RESTART_SOL= YES
READ_BINARY_RESTART= NO
RESTART_ITER= 1
SOLUTION_FILENAME= seed
TIME_ITER= 2
INNER_ITER= 120
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CONV_NUM_METHOD_FLOW= FDS
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
CFL_NUMBER= 3
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 0.001
LINEAR_SOLVER_ITER= 16
CONV_FIELD= (RMS_PRESSURE, RMS_VELOCITY-X, RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY, RMS_TKE, RMS_DISSIPATION, RMS_INTERMITTENCY, RMS_RE_THETA_T)
CONV_RESIDUAL_MINVAL= -8
MESH_FILENAME= mesh.su2
MESH_FORMAT= SU2
VOLUME_OUTPUT= (COORDINATES, SOLUTION, INTERMITTENCY_SEP, INTERMITTENCY_EFF, TEMPERATURE, HEAT_FLUX, Y_PLUS)
HISTORY_OUTPUT= (ITER, RMS_RES, HEAT)
SCREEN_OUTPUT= (TIME_ITER, INNER_ITER, RMS_PRESSURE, RMS_ENTHALPY, RMS_TKE, RMS_INTERMITTENCY)
HISTORY_WRT_FREQ_INNER= 1
SCREEN_WRT_FREQ_INNER= 5
WRT_RESTART_COMPACT= NO
OUTPUT_FILES= (RESTART_ASCII, SURFACE_PARAVIEW_ASCII)
OUTPUT_WRT_FREQ= (1, 1)
RESTART_FILENAME= restart
SURFACE_FILENAME= cylinder_surface
"""


def prepare(folder):
    folder.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(BASE/"condition01/condition.json",folder/"condition.json")
    c=json.loads((folder/"condition.json").read_bytes())
    xyz,info=mesh(folder);initial(folder,xyz,c)
    (folder/"case.cfg").write_text(configuration(c),encoding="ascii")
    save(folder/"mesh.json",info)
    paths=[EXE,Path(__file__).resolve(),ROOT/"lewisCylinderCondition.m",folder/"condition.json",
           folder/"wall.npz",folder/"case.cfg",folder/"mesh.su2",folder/"seed_00000.csv"]
    assert sha(EXE)==EXE_SHA
    save(folder/"inputs.json",{str(p):sha(p) for p in paths})


def solve(folder):
    inputs=json.loads((folder/"inputs.json").read_bytes())
    assert all(sha(Path(p))==h for p,h in inputs.items()) and sha(EXE)==EXE_SHA
    command=[str(EXE),"case.cfg","--threads","2"]
    start=time.monotonic();peak=0;disk=0;reason=None
    with (folder/"native.log").open("xb") as log:
        proc=subprocess.Popen(command,cwd=folder,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW,
                              env=dict(os.environ,OMP_NUM_THREADS="2"))
        try:
            save(folder/"started.json",dict(pid=proc.pid,command=command,inputsSHA256=sha(folder/"inputs.json")))
            handle=psutil.Process(proc.pid)
            while proc.poll() is None:
                try: peak=max(peak,handle.memory_info().rss)
                except psutil.NoSuchProcess: continue
                disk=max(disk,sum(p.stat().st_size for p in folder.iterdir() if p.is_file()))
                if peak>2*1024**3 or disk>240*1024**2 or time.monotonic()-start>=420:
                    reason="resource budget";break
                time.sleep(.2)
        finally:
            if proc.poll() is None:proc.kill()
            proc.wait()
    disk=max(disk,sum(p.stat().st_size for p in folder.iterdir() if p.is_file()))
    if disk>240*1024**2 and reason is None:reason="final output exceeds disk budget"
    receipt=dict(exitCode=proc.returncode,stopReason=reason,elapsedSeconds=time.monotonic()-start,
                 peakSampledRssBytes=peak,peakSampledFolderBytes=disk,pid=proc.pid,command=command,
                 inputsStable=all(sha(Path(p))==h for p,h in inputs.items()),
                 outputs={p.name:sha(p) for p in folder.iterdir() if p.is_file() and str(p) not in inputs})
    save(folder/"terminal.json",receipt)
    if proc.returncode or reason or not receipt["inputsStable"]:raise RuntimeError("原生试算终止，请查看本次回执")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("prepare","solve"))
    parser.add_argument("folder",type=Path)
    args=parser.parse_args();folder=args.folder.resolve();folder.relative_to(BASE)
    (prepare if args.action=="prepare" else solve)(folder)
