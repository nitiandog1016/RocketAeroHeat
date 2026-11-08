"""从既有真实棱柱层提取有母单元映射的小控制；不改原网格坐标或层数。"""
from pathlib import Path
import json,time,shutil,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();case=BASE/'control01';case.mkdir(exist_ok=True)
    source=ROOT/'outputs/v208_validation/layer_mesh02/exact_mesh.npz'
    assert sha(source)=='ae2f2d7224e44e5cf8f11984f4032248d9c2a12d90860fb92fcf2b5d9ff00ec7'
    z=np.load(str(source));P=z['points'];C=z['cells'];nw=len(z['originalWallFaceIDs']);nv=len(P)//33
    wall=C[:nw,:3];center=P[wall].mean(axis=1)
    pick=int(np.argmin(np.linalg.norm(center-[-.46,0,.12],axis=1)))
    cols=np.flatnonzero(np.isin(wall,wall[pick]).any(axis=1))
    selected=np.concatenate([cols+k*nw for k in range(32)])
    assert np.all(C[selected,:3]%nv==np.tile(wall[cols],(32,1)))
    selected_mask=np.zeros(len(C),bool);selected_mask[selected]=True
    nodes=np.unique(C[selected]);local=np.full(len(P),-1,int);local[nodes]=np.arange(len(nodes))
    owner=z['owner'];other=np.full(len(owner),-1,int);other[:len(z['neighbour'])]=z['neighbour']
    om=selected_mask[owner];nm=np.zeros(len(owner),bool);nm[other>=0]=selected_mask[other[other>=0]]
    boundaries=np.flatnonzero(om^nm);tags=dict(wall=[],outer=[],exchange=[]);donors=[]
    for face_id in boundaries:
        raw=z['faces'][face_id,:z['faceSizes'][face_id]];face=local[raw];assert np.all(face>=0)
        if other[face_id]>=0:
            donor=other[face_id] if om[face_id] else owner[face_id]
            tags['exchange'].append(face.tolist())
            donors.append(dict(faceVertices=sorted(face.tolist()),motherCellID=int(z['originalWedgeRows'][donor]),
                vertices=P[C[donor]].tolist()))
        elif np.all(raw<nv):tags['wall'].append(face.tolist())
        elif np.all(raw>=32*nv):tags['outer'].append(face.tolist())
        else:raise ValueError('控制触及原区域侧界，没有完整母donor')
    assert len(tags['wall'])==len(cols)==len(tags['outer'])
    # 模块一：只重编号控制文件；所有顶点坐标/棱柱连接直接取原资产。
    with (case/'mesh.su2').open('x') as stream:
        stream.write('NDIME= 3\nNELEM= %d\n'%len(selected))
        for i,row in enumerate(local[C[selected]]):stream.write('13 '+' '.join(map(str,row))+' '+str(i)+'\n')
        stream.write('NPOIN= %d\n'%len(nodes))
        for i,row in enumerate(P[nodes]):stream.write(' '.join(format(x,'.17g') for x in row)+' '+str(i)+'\n')
        stream.write('NMARK= 3\n')
        for tag,rows in tags.items():
            stream.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(tag,len(rows)))
            for row in rows:stream.write(('5 ' if len(row)==3 else '9 ')+' '.join(map(str,row))+'\n')
    (case/'donors.json').write_text(json.dumps(dict(schema='RAH.NativeFrozenControl.v1',
        origin=P[wall[pick]].mean(axis=0).tolist(),donors=donors),indent=2)+'\n')
    (case/'grid.lua').write_text('''-- 模块二：原生只读取已有棱柱控制，保持楔体节点顺序。
config.dimensions=3
config.grid_format="rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="mesh.su2",fmt="su2text",scale=1.0,expect_gmsh_order_for_wedges=false},fsTag="initial"}
''',encoding='utf8')
    (case/'phase-air.lua').write_text('model="rah-phase-aware-dry-air-eilmer-v1"\n')
    (case/'job.lua').write_text('''-- 模块三：仅原生算子控制，禁止牛顿/时间发展；人工层顶不是物理终端资格。
config.dimensions=3
config.solver_mode="steady"
config.field_format="rawbinary"
config.grid_format="rawbinary"
config.viscous=true
config.reacting=false
config.turbulence_model="spalart_allmaras_edwards"
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
setGasModel("phase-air.lua")
local initial=FlowState:new{p=2500.0,T=300.0,velx=-1000.0,vely=0.0,velz=200.0,nuhat=1.0e-5}
makeFluidBlocks({wall=WallBC_NoSlip_FixedT:new{Twall=300.0,group="wall"},
 outer=InFlowBC_Supersonic:new{flowCondition=initial},
 exchange=InFlowBC_Supersonic:new{flowCondition=initial,label="exchange"}}, {initial=initial})
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.thermo_interpolator="rhou"
config.extrema_clipping=false
config.unstructured_limiter="hvenkat"
config.smooth_limiter_coeff=1.0
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
config.viscous_least_squares_type="weighted_qr"
NewtonKrylovGlobalConfig{number_of_phases=1,max_newton_steps=0,use_preconditioner=true,
 preconditioner="ilu",ilu_fill=0,max_linear_solver_iterations=50,max_linear_solver_restarts=0,
 use_physicality_check=true,use_line_search=true,stop_on_relative_residual=1.0e-8}
NewtonKrylovPhase:new{residual_interpolation_order=2,jacobian_interpolation_order=1,
 use_residual_smoothing=false,start_cfl=0.1}
''',encoding='utf8')
    np.savez(case/'mapping.npz',motherCellRows=z['originalWedgeRows'][selected],motherNodeIDs=z['originalNodeIDs'][nodes],
        originalWallFaceIDs=z['originalWallFaceIDs'][cols],points=P[nodes],cells=local[C[selected]])
    report=dict(passed=True,selectedColumns=len(cols),cells=len(selected),vertices=len(nodes),
        boundaries={k:len(v) for k,v in tags.items()},motherMeshSHA256=sha(source),layers=32,
        firstHeightM=2e-6,CFDRuns=0,physicalQualification=False,seconds=time.monotonic()-start)
    (case/'design_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':raise ValueError('v213控制提取已结束；首节点约定失败保留，复用control02，不重新选域。')
