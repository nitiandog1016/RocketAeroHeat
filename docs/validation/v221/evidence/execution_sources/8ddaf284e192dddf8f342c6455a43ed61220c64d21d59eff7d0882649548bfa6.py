"""模块职责：绑定新SU2能量身份与实际保存文件；未知诊断字段不变成已解结果。"""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from su2EnergyReference import contract as energy_contract,flow_packet,convert_flow,SOLVER_REFERENCE,PHYSICAL_REFERENCE
from phaseAwareViscousAir import PhaseAwareViscousAir
MODEL='rah-su2-phase-air-saneg-eplus300k-v1'
SCHEMA='rah-su2-phase-native-frame-v1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def bind(path,mesh,config,engine,producer,role):
    """只绑定已存在文件；调用者须先核对实际运行源，绑定本身不授数值资格。"""
    if role not in ('manufactured-control','initialized-state','research-state'):raise ValueError('FRAME_ROLE')
    return dict(schema=SCHEMA,model=MODEL,energy=energy_contract(),meshSHA256=sha(mesh),
        sourceConfigSHA256=sha(config),engineSHA256=sha(engine),producerSHA256=sha(producer),
        payloadSHA256=sha(path),role=role,physicalQualification=False,
        turbulence='SA-neg/WITHFT2',saVariable='nonconservative-nu-tilde',wallTemperatureRole='prescribed',
        sourceThermoSHA256=sha(Path(__file__).resolve().parents[1]/'data/thermo/phase_aware_dry_air_v1.json'))

def read_restart(path,metadata,expected,mesh):
    """核对全部身份再读取守恒态；拒绝旧包、无标记、篡改文件或其他网格。"""
    if not isinstance(metadata,dict) or set(metadata)!=set(expected):raise ValueError('FRAME_SCHEMA')
    for name,value in expected.items():
        if metadata[name]!=value:raise ValueError('FRAME_IDENTITY_'+name)
    if metadata['schema']!=SCHEMA or metadata['model']!=MODEL or metadata['energy']!=energy_contract():raise ValueError('FRAME_MODEL')
    if metadata['payloadSHA256']!=sha(path) or metadata['meshSHA256']!=sha(mesh):raise ValueError('FRAME_BYTES')
    with Path(path).open(encoding='utf8') as stream:rows=list(csv.DictReader(stream,skipinitialspace=True))
    names=['Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    try:
        ids=np.array([int(r['PointID']) for r in rows]);xyz=np.array([[float(r[k]) for k in ('x','y','z')] for r in rows])
        native=np.array([[float(r[k]) for k in names] for r in rows]);sa=np.array([float(r['Nu_Tilde']) for r in rows])
    except (KeyError,TypeError,ValueError):raise ValueError('FRAME_COLUMNS')
    if len(ids)==0 or len(set(ids.tolist()))!=len(ids) or np.any(ids<0) or not np.isfinite(xyz).all() or not np.isfinite(sa).all():raise ValueError('FRAME_IDS')
    packet=convert_flow(flow_packet(native,SOLVER_REFERENCE,'state'),PHYSICAL_REFERENCE)
    # 比焓零点已转换；热状态仍须同物性反算，不能仅有限值就接受。
    air=PhaseAwareViscousAir();rho=native[:,0];vel=native[:,1:4]/rho[:,None]
    e=native[:,4]/rho-.5*np.sum(vel*vel,axis=1)-300000
    from scipy.optimize import brentq
    lo,hi=[float(air.thermal(t)[1]-air.R*t) for t in (60.,2000.)]
    if np.any((e<lo)|(e>hi)):raise ValueError('FRAME_ENERGY_DOMAIN')
    t=np.array([brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-v,60.,2000.,xtol=1e-10) for v in e])
    p=rho*air.R*t;air.evaluate(p,t)
    # 原生原始输出可能含未计算的Heat_Flux/Y_Plus=0。此接口故意不读取这些诊断列。
    return dict(ids=ids,xyz=xyz,nativeConservative=native,physicalConservative=packet['values'],
        nuTilde=sa,pressurePa=p,fluidTemperatureK=t,role=metadata['role'],physicalQualification=False,
        wallHeatFlux=None,solidTemperature=None)
