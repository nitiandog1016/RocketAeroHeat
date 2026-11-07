"""生成具名Run14层流热控制；固定物理输入和官方模型，不修改真实箭体默认值。"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import time

COMMIT='f53f4609a0331d48efee69a4e4f3c3598378cc03'


def prepare(out):
    start=time.monotonic();out=Path(out).resolve();assert not out.exists();out.mkdir(parents=True)
    # 模块一：物理尺寸由原论文图1核对，220mm末端保持，不拉长迎合越域压力测站。
    grid='''-- Run14物理外壁与官方两块外域；坐标为本验证算例局部轴向+X。
config.dimensions=2
config.axisymmetric=true
config.grid_format="gziptext"
local L1,L2,R=0.1017,0.220,0.0325
local rise=math.tan(math.rad(30))*(L2-L1)
local a0=Vector3:new{x=0,y=R};local a1=Vector3:new{x=0,y=R+0.005}
local b0=Vector3:new{x=L1,y=R};local b1=Vector3:new{x=L1-0.005,y=R+0.020}
local c0=Vector3:new{x=L2,y=R+rise};local c1=Vector3:new{x=L2,y=R+rise+0.025}
local wallCluster=RobertsFunction:new{end0=true,end1=false,beta=1.01}
local noseCluster=RobertsFunction:new{end0=true,end1=false,beta=1.2}
registerFluidGrid{grid=StructuredGrid:new{
  psurface=CoonsPatch:new{p00=a0,p10=b0,p11=b1,p01=a1},niv=161,njv=97,
  cfList={north=noseCluster,south=noseCluster,east=wallCluster,west=wallCluster}},
  fsTag="inflow",bcTags={north="inflow",west="inflow",south="wall"}}
registerFluidGrid{grid=StructuredGrid:new{
  psurface=CoonsPatch:new{p00=b0,p10=c0,p11=c1,p01=b1},niv=161,njv=97,
  cfList={east=wallCluster,west=wallCluster}},
  fsTag="initial",bcTags={north="inflow",east="outflow",south="wall"}}
identifyGridConnections()
'''
    # 模块二：官方N2 Sutherland黏度+导热，消除120.4K读取200K端点的问题。
    gas='''-- 固定官方N2默认Sutherland输运；不传Prandtl以免覆盖导热模型。
model="IdealGas"
species={"N2"}
'''
    job='''-- 具名理想N2层流控制，忽略振动是明示模型假设，不替换真实箭体变cp模型。
config.dimensions=2
config.axisymmetric=true
config.solver_mode="steady"
config.field_format="rawbinary"
config.grid_format="gziptext"
config.viscous=true
config.turbulence_model="none"
config.reacting=false
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
config.flux_calculator="adaptive_hanel_ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
local _,_,gm=setGasModel("nitrogen.gas")
local gs=GasState:new{gm}
local f=assert(io.open("native-gas.csv","w"))
f:write("T,p,rho,a,e,cp,mu,k\\n")
for _,T in ipairs{100,120.4,150,200,295.2,500,1000,2000} do
 gs.p=31.88;gs.T=T;gs.massf={N2=1}
 gm:updateThermoFromPT(gs);gm:updateSoundSpeed(gs);gm:updateTransCoeffs(gs)
 f:write(string.format("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\\n",T,gs.p,gs.rho,gs.a,gs.u,gm:Cp(gs),gs.mu,gs.k))
end
f:close()
local inflow=FlowState:new{p=31.88,T=120.4,velx=2304,massf={N2=1}}
local initial=FlowState:new{p=31.88/5,T=120.4,massf={N2=1}}
makeFluidBlocks({inflow=InFlowBC_Supersonic:new{flowState=inflow},
 wall=WallBC_NoSlip_FixedT:new{Twall=295.2,group="wall"},
 outflow=OutFlowBC_FixedP:new{p_outside=31.88/5}}, {inflow=inflow,initial=initial})
config.boundary_groups_for_loads="wall"
-- 原生两阶段精度推进：50步一阶启动后转二阶；自动CFL按残差反馈，不扫描。
NewtonKrylovGlobalConfig{number_of_phases=2,max_steps_in_initial_phases={50},
 max_newton_steps=2000,stop_on_relative_residual=1e-8,
 min_cfl=0.01,max_cfl=1000,use_line_search=true,use_physicality_check=true,
 max_consecutive_bad_steps=2,use_preconditioner=true,preconditioner="ilu",ilu_fill=0,
 max_linear_solver_iterations=50,max_linear_solver_restarts=0,
 write_loads=true,steps_between_loads_update=50,steps_between_snapshots=250,total_snapshots=3,
 write_snapshot_on_last_step=true,write_loads_on_last_step=true,
 steps_between_status=10,steps_between_diagnostics=1}
NewtonKrylovPhase:new{use_auto_cfl=true,start_cfl=0.1,max_cfl=1000,
 threshold_relative_residual_for_cfl_growth=0.9,auto_cfl_exponent=0.75,
 residual_interpolation_order=1,jacobian_interpolation_order=1,linear_solve_tolerance=0.01}
NewtonKrylovPhase:new{start_cfl=1,residual_interpolation_order=2,jacobian_interpolation_order=1}
'''
    files={'grid.lua':grid,'gas.inp':gas,'job.lua':job}
    for name,text in files.items():(out/name).write_text(text,encoding='utf8',newline='\n')
    m=dict(schema='RAH.EilmerRun14Control.v1',sourceCommit=COMMIT,
        gas='calorically-perfect N2 / official Sutherland viscosity and conductivity',
        turbulence='laminar',vibrationIncluded=False,flowMode='explicitVelocity',
        pressurePa=31.88,temperatureK=120.4,velocityMps=2304.,wallTemperatureK=295.2,
        reportedDensityKgM3=.000881,reportedVibrationalTemperatureK=2467,
        cylinderLengthM=.1017,wallEndXM=.22,cylinderRadiusM=.0325,flareAngleDeg=30,
        originalMeasurementNormalizationM=.1017016,wallFaces=320,cells=30720,
        maxSteps=2000,initialOrderSteps=50,fieldSnapshotsEvery=250,loadsEvery=50,
        physicalQualification=False,newQualifiedAreaM2=0,
        excludedFromWallComparison=['P30','P32'],excludeReason='x exceeds physical 220mm outer-wall end; retained as unresolved',
        measurementUncertainty=None,finiteGaugeAreaKnown=False,
        experimentSource='https://chimeracfd.com/professional/public/aiaa2004_0529.pdf',
        filesSHA256={k:hashlib.sha256((out/k).read_bytes()).hexdigest() for k in files},
        elapsedSeconds=time.monotonic()-start)
    (out/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf8')
    shutil.copyfile(__file__,out/'executed_preparation.py');return m


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True)
    m=prepare(p.parse_args().output);print(json.dumps({'prepared':True,'cells':m['cells'],'physicsQualified':False}))
