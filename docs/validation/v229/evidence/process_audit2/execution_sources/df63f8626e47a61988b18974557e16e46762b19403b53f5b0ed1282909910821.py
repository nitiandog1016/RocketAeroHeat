"""模块职责：独立Euler五守恒量的文件身份/能量/全局编号合同；不改旧NS/SA合同。"""
from pathlib import Path
import csv
import numpy as np
from runSU2RocketPhaseV222 import sha,ROOT
from su2EnergyReference import contract as energy_contract,flow_packet,convert_flow,SOLVER_REFERENCE,PHYSICAL_REFERENCE
from phaseAwareViscousAir import PhaseAwareViscousAir

MODEL='rah-su2-phase-air-euler-eplus300k-v1'
SCHEMA='rah-su2-phase-euler-frame-v1'
ROLES=('manufactured-control','initialized-state','research-state')

def bind(path,mesh,config,engine,producer,role):
    """绑定来源不授物理资格，气体键/无湍流必须与独立原生构建一致。"""
    if role not in ROLES:raise ValueError('EULER_FRAME_ROLE')
    cfg={}
    for line in Path(config).read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    if cfg.get('FLUID_MODEL')!='RAH_PHASE_AIR_EULER' or cfg.get('SOLVER')!='EULER' or cfg.get('KIND_TURB_MODEL','NONE')!='NONE':raise ValueError('EULER_CONFIG_IDENTITY')
    return dict(schema=SCHEMA,model=MODEL,energy=energy_contract(),meshSHA256=sha(mesh),sourceConfigSHA256=sha(config),
      engineSHA256=sha(engine),producerSHA256=sha(producer),payloadSHA256=sha(path),role=role,physicalQualification=False,
      turbulence='none',conservativeVariables=5,wallHeatFluxRole='not-solved',solidTemperatureRole='not-solved',
      sourceThermoSHA256=sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json'))

def read_restart(path,metadata,expected,mesh):
    """逐项身份后核对全局ID和实际热状态；原生诊断零值不读取为热结果。"""
    if not isinstance(metadata,dict) or set(metadata)!=set(expected):raise ValueError('EULER_FRAME_SCHEMA')
    for k,v in expected.items():
        if metadata[k]!=v:raise ValueError('EULER_FRAME_IDENTITY_'+k)
    fixed=dict(schema=SCHEMA,model=MODEL,energy=energy_contract(),physicalQualification=False,turbulence='none',
      conservativeVariables=5,wallHeatFluxRole='not-solved',solidTemperatureRole='not-solved')
    if any(metadata.get(k)!=v for k,v in fixed.items()) or metadata['role'] not in ROLES:raise ValueError('EULER_FRAME_MODEL')
    if metadata['payloadSHA256']!=sha(path) or metadata['meshSHA256']!=sha(mesh):raise ValueError('EULER_FRAME_BYTES')
    with Path(path).open(encoding='utf8') as f:rows=list(csv.DictReader(f,skipinitialspace=True))
    try:
        ids=np.array([int(r['PointID']) for r in rows]);xyz=np.array([[float(r[k]) for k in ('x','y','z')] for r in rows])
        U=np.array([[float(r[k]) for k in ('Density','Momentum_x','Momentum_y','Momentum_z','Energy')] for r in rows])
    except (KeyError,TypeError,ValueError):raise ValueError('EULER_FRAME_COLUMNS')
    if not len(ids) or len(set(ids.tolist()))!=len(ids) or np.any(ids<0) or not np.isfinite(xyz).all() or not np.isfinite(U).all():raise ValueError('EULER_FRAME_IDS')
    if np.any(U[:,0]<=0):raise ValueError('EULER_FRAME_DENSITY')
    lines=Path(mesh).read_text(encoding='utf8').splitlines();starts=[i for i,s in enumerate(lines) if s.strip().startswith('NPOIN=')]
    if len(starts)!=1:raise ValueError('EULER_MESH_FORMAT')
    start=starts[0];count=int(lines[start].split('=',1)[1].split()[0]);coordinates=np.empty((count,3));seen=set()
    for line in lines[start+1:start+1+count]:
        fields=line.split()
        if len(fields)!=4:raise ValueError('EULER_MESH_DIMENSION')
        ident=int(fields[3])
        if ident in seen or ident<0 or ident>=count:raise ValueError('EULER_MESH_IDS')
        seen.add(ident);coordinates[ident]=[float(v) for v in fields[:3]]
    if len(seen)!=count or len(ids)!=count or set(ids.tolist())!=seen:raise ValueError('EULER_FRAME_COVERAGE')
    if not np.array_equal(xyz,coordinates[ids]):raise ValueError('EULER_FRAME_COORDINATES')
    physical=convert_flow(flow_packet(U,SOLVER_REFERENCE,'state'),PHYSICAL_REFERENCE)['values']
    air=PhaseAwareViscousAir();rho=U[:,0];v=U[:,1:4]/rho[:,None];e=U[:,4]/rho-.5*np.sum(v*v,axis=1)-300000
    from scipy.optimize import brentq
    lo,hi=[float(air.thermal(t)[1]-air.R*t) for t in (60.,2000.)]
    if np.any((e<lo)|(e>hi)):raise ValueError('EULER_FRAME_ENERGY_DOMAIN')
    T=np.array([brentq(lambda t:float(air.thermal(t)[1]-air.R*t)-x,60.,2000.,xtol=1e-10) for x in e])
    p=rho*air.R*T;air.evaluate(p,T)
    return dict(ids=ids,xyz=xyz,nativeConservative=U,physicalConservative=physical,pressurePa=p,fluidTemperatureK=T,
      velocityMps=v,role=metadata['role'],physicalQualification=False,wallHeatFlux=None,solidTemperature=None)
