# -*- coding: utf-8 -*-
"""生成具名压缩—膨胀Euler输入；本模块不运行原生程序、不准入长算。"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

COMMIT = 'f53f4609a0331d48efee69a4e4f3c3598378cc03'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(folder, grid_intervals=(60,2,20)):
    """模块一：原X约定下挤出10度压缩/10度回转控制，所有端面具名。"""
    folder = Path(folder).resolve()
    if folder.exists():
        raise ValueError('目标已存在，禁止覆盖或当作自动续算')
    if tuple(grid_intervals) not in ((60,2,20),(120,2,40)):
        raise ValueError('当前只登记原控制或流向/法向各加密2倍的控制')
    nx, ny, nz = grid_intervals
    ss = np.linspace(0, 1, nx+1)
    yy = np.linspace(-.02, .02, ny+1)
    ee = np.linspace(0, 1, nz+1)
    wall = np.clip(ss-.2, 0, .4)*math.tan(math.radians(10))
    points = np.array([[-s, y, w+e*(.6-w)] for s, w in zip(ss, wall) for y in yy for e in ee])
    def node(i,j,k): return (i*(ny+1)+j)*(nz+1)+k
    cells=[]; markers={name:[] for name in ('wall','outer','inflow','outflow','span')}
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                a=[node(i+1,j,k),node(i,j,k),node(i,j+1,k),node(i+1,j+1,k),
                   node(i+1,j,k+1),node(i,j,k+1),node(i,j+1,k+1),node(i+1,j+1,k+1)]
                cells.append(a)
                # 下列均为标准六面体外向面，不用最近邻推定边界身份。
                for active, name, order in (
                    (k==0,'wall',(0,3,2,1)),(k==nz-1,'outer',(4,5,6,7)),
                    (i==0,'inflow',(1,2,6,5)),(i==nx-1,'outflow',(3,0,4,7)),
                    (j==0,'span',(0,1,5,4)),(j==ny-1,'span',(2,3,7,6))):
                    if active: markers[name].append([a[q] for q in order])
    folder.mkdir(parents=True)
    with (folder/'mesh.su2').open('w',encoding='ascii',newline='\n') as stream:
        stream.write('NDIME= 3\nNELEM= %d\n'%len(cells))
        for a in cells: stream.write('12 '+' '.join(map(str,a))+'\n')
        stream.write('NPOIN= %d\n'%len(points))
        for i,p in enumerate(points): stream.write(' '.join(format(v,'.17g') for v in p)+' '+str(i)+'\n')
        stream.write('NMARK= %d\n'%len(markers))
        for name,faces in markers.items():
            stream.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(name,len(faces)))
            for a in faces: stream.write('9 '+' '.join(map(str,a))+'\n')
    # 模块二：固定成熟Euler算子。无湍流/热壁项，不复用旧Wilcox准入。
    files={
        'gas.inp':'model = "ThermallyPerfectGas"\nspecies = {"air"}\n',
        'grid.lua':'''-- 新挤出控制网格，非真实20度箭体。
config.dimensions = 3
config.grid_format = "rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="mesh.su2",fmt="su2text",
 scale=1.0,expect_gmsh_order_for_wedges=false},fsTag="initial"}
''',
        'job.lua':'''-- 仅生成，未原生验证。局部伪时间不代表飞行时间。
config.dimensions = 3
config.solver_mode = "time-marching"
config.field_format = "rawbinary"
config.grid_format = "rawbinary"
config.viscous = false
config.reacting = false
config.turbulence_model = "none"
config.adjust_invalid_cell_data = false
config.max_invalid_cells = 0
local nsp,nmodes,gm=setGasModel("rah-air.gas")
local gs = GasState:new{gm}
gs.p=5000.0; gs.T=300.0; gs.massf={air=1.0}
gm:updateThermoFromPT(gs); gm:updateSoundSpeed(gs); gm:updateTransCoeffs(gs)
local inflow=FlowState:new{p=gs.p,T=gs.T,massf=gs.massf,velx=-3.0*gs.a,vely=0.0,velz=0.0}
local receipt=assert(io.open("native-input-state.json","w"))
receipt:write(string.format('{"p":%.17g,"T":%.17g,"rho":%.17g,"e":%.17g,"a":%.17g,"cp":%.17g,"velx":%.17g}',
 gs.p,gs.T,gs.rho,gs.u,gs.a,gm:Cp(gs),-3.0*gs.a)); receipt:close()
makeFluidBlocks({wall=WallBC_WithSlip:new{group="wall"},span=WallBC_WithSlip:new{},
 outer=InOutFlowBC_Ambient:new{flowState=inflow},
 inflow=InFlowBC_Supersonic:new{flowCondition=inflow},outflow=OutFlowBC_Simple:new{}},
 {initial=inflow})
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
config.gasdynamic_update_scheme="tvd_rk3"
config.with_local_time_stepping=true
config.local_time_stepping_limit_factor=1000
config.with_super_time_stepping=false
config.residual_smoothing=false
config.stringent_cfl=true
config.cfl_value=0.3
config.cfl_count=1
config.fixed_time_step=false
config.dt_init=1.0e-10
config.dt_max=1.0e-3
config.max_time=1.0
config.max_step=400
config.max_attempts_for_step=1
config.print_count=1
config.write_transient_residuals=true
config.dt_plot_value=1.0
config.write_flow_solution_at_step=200
'''}
    for name, text in files.items(): (folder/name).write_bytes(text.encode('utf-8'))
    triangles=[]
    for a,b,c,d in markers['wall']: triangles += [[a,b,c],[a,c,d]]
    surface=dict(vertices=points.tolist(),facesOneBased=(np.array(triangles)+1).tolist(),
                 faceIDs=list(range(1,len(triangles)+1)),identity='new extruded-control wall; not original rocket IDs')
    (folder/'surface.json').write_text(json.dumps(surface,separators=(',',':'))+'\n')
    manifest=dict(schema='RAH.HybridEulerPreparation.v1',candidate='hybrid-compression-expansion-control-01',
        sourceCommit=COMMIT,preparedOnly=True,nativePrepared=False,executionAuthorized=False,
        model='eilmer-air-cea-v1',energyReference='fixed official air CEA polynomial reference',
        condition=dict(Mach=3.,pressurePa=5000.,temperatureK=300.,velocityDirection=[-1.,0.,0.]),
        geometry=dict(cells=len(cells),vertices=len(points),wallQuads=len(markers['wall']),
                      wallTriangles=len(triangles),gridIntervals=[nx,ny,nz],full3DTarget=False,
                      analyticalVolumeM3=.04*(.6-.24*math.tan(math.radians(10)))),
        nativeConfig=dict(viscous=False,turbulence_model='none',gasdynamic_update_scheme='tvd_rk3',
                          max_step=400,write_flow_solution_at_step=200,cfl_value=.3),
        pending=['native prep/API and actual gas confirmation','same-EOS shock/expansion reference',
                 'ghost and complete memory accounting','finite native budget and admission',
                 'boundary/refinement sensitivity and convergence'],
        newPhysicalQualification=False,newSolvedHeatAreaM2=0,wallHeatFluxAvailable=False,
        filesSHA256={n:sha(folder/n) for n in ['mesh.su2','surface.json']+list(files)})
    (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    result=prepare(parser.parse_args().out)
    print(json.dumps(dict(preparedOnly=True,geometry=result['geometry'],executionAuthorized=False)))
