# -*- coding: utf-8 -*-
"""Eilmer Euler生产场读取与壁邻映射；连接关系决定身份，数值有限不授外缘资格。"""
from pathlib import Path
import json
import numpy as np
from readEilmerNative import read_grid, compare_source, read_field, sha, need
from eilmerAirThermo import bind_model, evaluate
from hybridExternalState import exchange, fingerprint


def audit_inputs(case, preparation):
    """模块一：逐文件准备来源、实际算子和边界合同；与旧RANS入口隔离。"""
    case=Path(case);receipt=json.loads(Path(preparation).read_text())
    need(receipt.get('schema')=='RAH.HybridEulerNativePreparation.v1' and receipt.get('passed'),
         '缺少已通过Euler准备回执')
    # receipt可在Windows复核WSL路径；不放行别的case身份。
    def local(p):
        if str(case).startswith('/'): return Path(p)
        return Path(p.replace('/mnt/e/','E:/',1))
    need(local(receipt['case']).resolve()==case.resolve(),'准备回执算例身份错误')
    for p,digest in receipt['nativeInputSHA256'].items():
        need(sha(local(p))==digest,'准备输入变化: '+p)
    cfg=json.loads((case/'lmrsim/config').read_text());ctl=json.loads((case/'lmrsim/control').read_text())
    expected=dict(solver_mode='time-marching',grid_format='rawbinary',field_format='rawbinary',
        viscous=False,reacting=False,turbulence_model='none',dimensions=3,axisymmetric=False,
        gasdynamic_update_scheme='tvd_rk3',flux_calculator='ausmdv',interpolation_order=2,
        spatial_deriv_calc='least_squares',spatial_deriv_locn='cells',with_local_time_stepping=True,
        with_super_time_stepping=False,residual_smoothing=False,adjust_invalid_cell_data=False,max_invalid_cells=0)
    for key,val in expected.items():need(cfg.get(key)==val,'实际Euler配置不符: '+key)
    need(cfg['cfl_schedule_values']==[.3] and ctl['max_step']==400 and
         ctl['write_flow_solution_at_step']==200,'推进设计变化')
    model=bind_model(case/'rah-air.gas',sha(case/'rah-air.gas'))
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,error=compare_source(grid,case/'mesh.su2')
    need((len(grid['cells']),len(grid['points']),len(grid['faces']))==(2400,3843,8560),'控制规模改变')
    bc=dict(wall='wall_with_slip',span='wall_with_slip',inflow='inflow_supersonic',
            outflow='outflow_simple_flux',outer='inoutflow_ambient')
    need(set(grid['boundaries'])==set(bc),'边界集合改变')
    for name,kind in bc.items():
        b=cfg['block_0']['boundary_%d'%grid['boundaries'][name]['index']]
        need(b['type']==kind and not b['is_wall_with_viscous_effects'],'Euler边界错误: '+name)
    field=read_field(case,0,len(grid['cells']))
    volume=float(field['vol'].sum());need(np.all(field['vol']>0),'非正原生体积')
    analytical=.04*(.6-.24*np.tan(np.radians(10)))
    need(abs(volume-analytical)<2e-13,'原生体积与解析不符')
    ng=sum(len(b['ids']) for b in grid['boundaries'].values())
    # 固定ufluidblock按每外边界一ghost分配；不把结构网格的层数乘进来。
    return dict(passed=True,cells=2400,vertices=3843,faces=8560,ghostCells=ng,
                cellObjectsIncludingGhost=2400+ng,volumeM3=volume,sourceCoordinateError=error,
                boundaryFaces={k:len(b['ids']) for k,b in grid['boundaries'].items()},
                actualGasModel=model,actualNumerics=expected,
                initialTemperatureK=[float(field['T'].min()),float(field['T'].max())],
                physicalQualification=False),grid,mapping


def read_wall_frame(case, snapshot, preparation):
    """模块二：原壁三角面→原生四边形→唯一相邻单元；不是面重构/边界层外缘。"""
    case=Path(case);audit,grid,mapping=audit_inputs(case,preparation)
    fields=read_field(case,snapshot,len(grid['cells']))
    thermo=evaluate(fields['p'],fields['T']);errors={}
    for key,native in [('densityKgm3','rho'),('internalEnergyJkg','e'),('soundSpeedMps','a')]:
        error=float(np.max(abs(thermo[key]-fields[native])/np.maximum(1,abs(thermo[key]))))
        need(error<2e-10,'原生EOS不相容: '+native);errors[native]=error
    need(np.all(fields['massf-air']==1),'非纯伪空气')
    surface=json.loads((case/'surface.json').read_text());xyz=np.array(surface['vertices'])
    need(np.array_equal(xyz,grid['points']),'表面顶点与原生不一致')
    source_quads=mapping['wall']['originalNodes'];triangles=[]
    for a,b,c,d in source_quads: triangles.extend([[a,b,c],[a,c,d]])
    need(np.array_equal(np.array(triangles)+1,surface['facesOneBased']),'表面三角剖分与原壁不一致')
    rows=mapping['wall']['nativeRows'];boundary=grid['boundaries']['wall']
    owners=np.repeat(boundary['owners'][rows],2);fid=np.repeat(boundary['ids'][rows],2)
    vel=np.column_stack([fields['vel.'+k] for k in 'xyz']);p=fields['p'][owners];T=fields['T'][owners]
    g=evaluate(p,T);out={key:g[key].tolist() for key in ['pressurePa','temperatureK','densityKgm3',
                          'internalEnergyJkg','enthalpyJkg','entropyJkgK']}
    out['velocityMps']=vel[owners].tolist()
    out['totalEnthalpyJkg']=(g['enthalpyJkg']+np.sum(vel[owners]**2,axis=1)/2).tolist()
    v=xyz[np.array(triangles)];areas=np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)/2
    condition=json.loads((case/'manifest.json').read_text())['condition'];model=audit['actualGasModel']
    source=[case/'lmrsim/snapshots'/('%04d'%snapshot)/'fluid-0000',case/'lmrsim/snapshots/fluid.metadata',
            case/'lmrsim/grid/grid-0000',case/'lmrsim/config',case/'lmrsim/control',case/'rah-air.gas',
            case/'manifest.json',case/'surface.json',Path(preparation),Path(__file__),
            Path(__file__).with_name('eilmerAirThermo.py'),Path(__file__).with_name('hybridExternalState.py'),
            Path(__file__).with_name('readEilmerNative.py')]
    sources={str(p.resolve()):sha(p) for p in source}
    frame=dict(schema='RAH.HybridEulerFrame.v1',solver='Euler',turbulenceModel='none',
               surfaceSHA256=fingerprint(surface),conditionSHA256=fingerprint(condition),modelSHA256=fingerprint(model),
               sourceFilesSHA256=sources,snapshot=snapshot,faceIDs=surface['faceIDs'],
               originalAreaM2=float(areas.sum()),faceCentersM=v.mean(axis=1).tolist(),fields=out)
    result=exchange(frame,surface,condition,model,sources,evaluate)
    result.update(sampling='adjacent-cell centre; no wall reconstruction or edge qualification',
                  nativeCellIDsZeroBased=owners.tolist(),nativeFaceIDsZeroBased=fid.tolist(),
                  cellCentersM=np.column_stack([fields['pos.'+k][owners] for k in 'xyz']).tolist(),
                  faceCentersM=v.mean(axis=1).tolist(),faceAreasM2=areas.tolist(),
                  fullVolumeEOSErrors=errors)
    return result,audit
