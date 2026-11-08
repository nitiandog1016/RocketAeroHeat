"""真实全壁原生场的严格研究读取：原面映射、温度身份与能量零点分离。"""
from pathlib import Path
import json,re
import numpy as np
from prepareOpenFoamControl import ROOT,sha
from openFoamAirThermo import evaluate,HF,MODEL,bind_engineering_source
from openFoamFields import read_state,audit_state
from openFoamConservation import geometry,boundary_values,native_mass_flux

PROFILE='of13-original-rocket-sync-development-v1'
BASE=ROOT/'outputs/v201_validation'
FIELDS=('p','T','rho','U','e','phi')
END_MARK=b'// ************************************************************************* //'


def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))


def local_path(value):
    s=str(value).replace('\\','/')
    for prefix in ['E:/RocketAeroHeat/','/mnt/e/RocketAeroHeat/']:
        if s.startswith(prefix):return ROOT/s[len(prefix):]
    return Path(s)


def mapped_geometry(case):
    """模块一：验证原面编号、连接、体积质心和法向距离，不采用空间近邻。"""
    g=geometry(case);m=read_json(Path(case)/'wall_mapping.json')
    surface=read_json(ROOT/'outputs/v194_validation/rocket_design03/surface.json')
    ids=np.array(m['nativeFaceIDsZeroBased'],int);owner=np.array(m['nativeOwnerCellsZeroBased'],int)
    vertices=np.array(surface['vertices']);connections=np.array(surface['facesOneBased'],int)-1
    if m['originalFaceIDs']!=surface['faceIDs'] or surface['faceIDs']!=list(range(1,19753)):raise ValueError('原面编号不完整')
    if len(ids)!=19752 or len(np.unique(ids))!=19752 or ids.min()<len(g['neighbor']) or ids.max()>=len(g['owner']):raise ValueError('原面双射越界/重复')
    if not np.array_equal(g['points'][:len(vertices)],vertices):raise ValueError('原顶点改变')
    if not np.array_equal(np.sort(g['faces'][ids],axis=1),np.sort(connections,axis=1)):raise ValueError('原连接映射错误')
    if not np.array_equal(g['owner'][ids],owner):raise ValueError('原生owner映射错误')
    if m['originalConnections']!=surface['facesOneBased']:raise ValueError('来源连接错误')
    centers=vertices[connections].mean(axis=1);positions=g['centers'][owner]
    if not np.allclose(centers,m['wallCentersM'],rtol=0,atol=1e-13):raise ValueError('原面中心错误')
    if not np.allclose(positions,m['samplePositionsM'],rtol=0,atol=1e-13):raise ValueError('实际体积质心错误')
    inward=-g['sf'][ids]/g['area'][ids,None]
    distance=np.einsum('ij,ij->i',positions-centers,inward)
    if np.any(distance<=0) or not np.allclose(distance,m['normalDistanceM'],rtol=0,atol=1e-13):raise ValueError('原面法向采样距离错误')
    if abs(g['area'][ids].sum()-1.9646292565715553)>1e-13:raise ValueError('原面积分母错误')
    if not np.allclose(g['area'][ids],m['areaM2'],rtol=0,atol=1e-13):raise ValueError('逐原面面积错误')
    return g,m,ids,owner,surface


def verify_identity(case,run,admission):
    """模块二：源任务、气体、工况及准确输入绑定；已修改脚本以执行副本核对。"""
    case=Path(case)
    if run.get('profile')!=PROFILE or admission.get('profile')!=PROFILE:raise ValueError('错误运行身份')
    if admission.get('case')!=case.relative_to(ROOT).as_posix():raise ValueError('源目录错误')
    process=run['process']
    if process.get('liveGroupAfterCleanup') or not process.get('inputStable') or not all(process['inputStable'].values()):raise ValueError('输入变化或任务未清点')
    if not process.get('logExportPassed'):raise ValueError('日志导出未验证')
    if sha(BASE/'rocket_run01/admission.json')!=run['admissionSHA256']:raise ValueError('准入回执SHA改变')
    required=['exact_mesh.npz','wall_mapping.json','patches.json','thermo_contract.json',
        'constant/physicalProperties','system/controlDict','system/fvSchemes','system/fvSolution',
        'constant/polyMesh/points','constant/polyMesh/faces','constant/polyMesh/owner',
        'constant/polyMesh/neighbour','constant/polyMesh/boundary','0/p','0/T','0/U']
    hashes={str(local_path(p).resolve()):h for p,h in admission['inputSHA256'].items()}
    for name in required:
        p=case/name
        if sha(p)!=hashes.get(str(p.resolve())):raise ValueError('实际输入SHA改变 '+name)
    model=read_json(case/'thermo_contract.json');source=bind_engineering_source()
    if model['model']!=MODEL or model['source']['sha256']!=source['sha256'] or model['hfJkg']!=HF:raise ValueError('气体/参考焓改变')
    condition=read_json(ROOT/'outputs/v194_validation/rocket_design03/manifest.json')['condition']
    if admission['condition']!=condition:raise ValueError('错误来流工况')
    u=np.array(condition['velocityDirection'])*condition['Mach']*float(evaluate(condition['pressurePa'],condition['temperatureK'])['soundSpeedMps'])
    if not np.allclose(u,admission['velocityMps'],rtol=0,atol=1e-10):raise ValueError('混用旧气体来流速度')
    return model


def read_frame(name,expected_time,diagnostic=False,run_override=None,admission_override=None):
    """模块三：异常默认拒绝；显式诊断只读取当时完整有效场，不转授外缘/热资格。"""
    if not re.fullmatch(r'[0-9.eE+-]+',name) or float(name)<=0 or float(name)!=expected_time:raise ValueError('源时刻不符/初场不能冒充发展场')
    case=BASE/'rocket01';run=run_override or read_json(BASE/'rocket_run01/report.json')
    admission=admission_override or read_json(BASE/'rocket_run01/admission.json')
    if not run.get('passed') and not diagnostic:raise ValueError('异常/未完成任务默认拒绝，仅允许显式只读诊断')
    model=verify_identity(case,run,admission)
    if name not in run['saved']:raise ValueError('没有监督绑定的完整保存场')
    saved=run['saved'][name]
    if saved['timeSeconds']!=expected_time or set(saved['fieldsSHA256'])!=set(FIELDS):raise ValueError('保存字段/时刻错误')
    folder=case/name
    for n in FIELDS:
        path=folder/n;data=path.read_bytes()
        if not data.rstrip().endswith(END_MARK) or sha(path)!=saved['fieldsSHA256'][n]:raise ValueError('字段未完成或SHA改变 '+n)
        head=data[:2400].decode('utf-8')
        location=re.search(r'\blocation\s+"([^"]+)";',head)
        if not location or float(location.group(1))!=expected_time:raise ValueError('字段location时刻错误 '+n)
        if not re.search(r'\bobject\s+'+re.escape(n)+r'\s*;',head):raise ValueError('字段名称错误')
    g,m,ids,owners,surface=mapped_geometry(case);state=read_state(folder,len(g['cells']));state_audit=audit_state(state)
    boundary={n:boundary_values(folder,n,state[n],g) for n in state};audit_state(boundary)
    native_mass_flux(folder,g) # 同时核对全部实际phi保存数量/边界，不将缺通量默认为零。
    rows=ids-len(g['neighbor'])
    for n in ['p','T']:
        if not np.allclose(boundary[n][rows],state[n][owners],rtol=2e-12,atol=1e-10):raise ValueError('壁面零梯度状态不符 '+n)
    p=state['p'][owners];T=state['T'][owners];rho=state['rho'][owners];u=state['U'][owners];es=state['e'][owners]
    gas=evaluate(p,T);engineering_e=es+HF;h=engineering_e+p/rho;h0=h+np.einsum('ij,ij->i',u,u)/2
    h_error=float(np.max(abs(h-gas['enthalpyJkg'])/np.maximum(1,abs(gas['enthalpyJkg']))))
    if h_error>2e-8:raise ValueError('工程焓/原生内能零点不相容')
    clock=folder/'uniform/time';snapshot=None;native_time=None
    if clock.is_file():
        text=clock.read_text(encoding='utf-8');snapshot=int(re.search(r'\bindex\s+(\d+);',text).group(1))
        native_time=float(re.search(r'\bvalue\s+([\deE.+-]+);',text).group(1))
        # timePrecision=12控制目录名，字典value保留17位；按实际格式核对，不能要求两个浮点逐位相同。
        if float(format(native_time,'.12g'))!=expected_time:raise ValueError('原生时间字典错误')
    fields=dict(pressurePa=p,fluidTemperatureK=T,densityKgm3=rho,velocityMps=u,nativeSensibleInternalEnergyJkg=es,
        engineeringInternalEnergyJkg=engineering_e,enthalpyJkg=h,totalEnthalpyJkg=h0,entropyJkgK=gas['entropyJkgK'])
    metadata=dict(schema='RAH.OpenFOAMOriginalWallResearch.v1',profile=PROFILE,timeSeconds=expected_time,
        nativeTimeSeconds=native_time,timeIndex=snapshot,nativeTimeDictionaryPresent=clock.is_file(),
        timeDictionarySHA256=sha(clock) if clock.is_file() else None,restartReady=False,
        status='interrupted-run-diagnostic' if not run['passed'] else 'unqualified-development',
        faceCount=19752,originalAreaM2=float(g['area'][ids].sum()),temperatureIdentity='wall-adjacent fluid; not Te/Tr/solid Tw',
        samplePositionKind='exact native-owner tetrahedron volume centroid',method='OpenFOAM 13 shockFluid / NASA7 / synchronized Euler / upwind / Kurganov',
        energyReference='engineering e = native es + hf; h = e + p/rho; h0 = h + |U|^2/2',hfJkg=HF,
        fieldsSHA256=saved['fieldsSHA256'],runReportSHA256=sha(BASE/'rocket_run01/report.json'),admissionSHA256=run['admissionSHA256'],
        stateAudit=state_audit,engineeringEnthalpyMaxScaled=h_error,attachedEdgeQualified=False,
        engineeringHeatReady=False,physicalQualification=False,wallHeatFluxWm2=None,solidWallTemperatureK=None,newQualifiedHeatAreaM2=0)
    return metadata,fields,g,m,ids,owners
