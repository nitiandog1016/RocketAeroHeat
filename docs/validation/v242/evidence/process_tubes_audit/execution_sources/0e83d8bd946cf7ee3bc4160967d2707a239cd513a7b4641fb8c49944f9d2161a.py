"""模块职责：批量读取具名SU2研究场，分清当前守恒态、原生诊断和给定壁温。"""
from pathlib import Path
import csv
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirContract import MODEL,SCHEMA
from su2EnergyReference import contract,SHIFT_JKG
from runSU2RocketPhaseV222 import sha

FLOW=['Density','Momentum_x','Momentum_y','Momentum_z','Energy']
RESIDUAL=['Residual_'+name for name in FLOW]
COLUMNS=['PointID','x','y','z']+FLOW+['Nu_Tilde','Pressure','Temperature','Mach','Pressure_Coefficient',
    'Velocity_x','Velocity_y','Velocity_z','Laminar_Viscosity','Skin_Friction_Coefficient_x',
    'Skin_Friction_Coefficient_y','Skin_Friction_Coefficient_z','Heat_Flux','Y_Plus','Eddy_Viscosity']+RESIDUAL+['Residual_Nu_Tilde']

def inverse_energy(air,e):
    """模块一：Cv分段线性表的稳定二次根；只定位合法分段，不裁能量或温度。"""
    e=np.asarray(e,dtype=float);knots=air.H-air.R*air.T
    if not e.size or not np.isfinite(e).all() or np.any((e<knots[0])|(e>knots[-1])):raise ValueError('FRAME_ENERGY_DOMAIN')
    i=np.minimum(np.searchsorted(knots,e,side='right')-1,len(knots)-2)
    delta=e-knots[i];cv=air.Cp[i]-air.R;slope=air.slopes[i]
    t=air.T[i]+2*delta/(cv+np.sqrt(cv*cv+2*slope*delta))
    # 精确端点从原表取值；这不是对非法输入的裁剪。
    t=np.where(e==knots[0],air.T[0],np.where(e==knots[-1],air.T[-1],t))
    air.thermal(t)
    return t

def read_fields(path,metadata,expected,coordinates):
    """模块二：先查来源和全局编号，再批量解析；旧无标记场不可混入。"""
    if not isinstance(metadata,dict) or set(metadata)!=set(expected):raise ValueError('FRAME_SCHEMA')
    for name,value in expected.items():
        if metadata[name]!=value:raise ValueError('FRAME_IDENTITY_'+name)
    if metadata['schema']!=SCHEMA or metadata['model']!=MODEL or metadata['energy']!=contract():raise ValueError('FRAME_MODEL')
    if metadata['payloadSHA256']!=sha(path):raise ValueError('FRAME_BYTES')
    with Path(path).open(encoding='utf8') as stream:columns=next(csv.reader(stream,skipinitialspace=True))
    if columns!=COLUMNS:raise ValueError('FRAME_COLUMNS_'+str(columns))
    data=np.loadtxt(str(path),delimiter=',',skiprows=1)
    if data.shape!=(len(coordinates),len(columns)) or not np.isfinite(data).all():raise ValueError('FRAME_SHAPE_OR_FINITE')
    ids=data[:,0].astype(np.int64)
    if not np.array_equal(ids,data[:,0]) or not np.array_equal(np.sort(ids),np.arange(len(coordinates))):raise ValueError('FRAME_IDS')
    if not np.array_equal(data[:,1:4],coordinates[ids]):raise ValueError('FRAME_COORDINATES')
    data=data[np.argsort(ids)];v={name:data[:,i] for i,name in enumerate(columns)}
    rho=v['Density']
    if np.any(rho<=0):raise ValueError('FRAME_DENSITY')
    velocity=np.column_stack([v['Momentum_'+a]/rho for a in ('x','y','z')])
    air=PhaseAwareViscousAir();e=v['Energy']/rho-.5*np.sum(velocity**2,axis=1)-SHIFT_JKG
    current_t=inverse_energy(air,e);current_p=rho*air.R*current_t;current=air.evaluate(current_p,current_t)
    primitive=air.evaluate(v['Pressure'],v['Temperature'])
    residual=np.column_stack([v[k] for k in RESIDUAL]);physical_residual=residual.copy()
    physical_residual[:,4]-=SHIFT_JKG*residual[:,0]
    return dict(values=v,current=current,currentVelocity=velocity,primitive=primitive,
        nativeResidual=residual,physicalResidual=physical_residual,physicalQualification=False,
        solidTemperature=None,wallHeatFlux=None)
