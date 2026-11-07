"""新NIST原全壁研究场读取：明确模型/能量参考，禁止隐式套旧NASA7。"""
from pathlib import Path
import json,re
import numpy as np
from prepareOpenFoamControl import ROOT,sha
from nistIntegratedAir import evaluate,contract,MODEL
from openFoamFields import read_state,audit_state
from openFoamConservation import boundary_values,native_mass_flux
from readOpenFoamRocketState import mapped_geometry,local_path,FIELDS,END_MARK

BASE=ROOT/'outputs/v202_validation';PROFILE='of13-nist-original-rocket-sync-v1'


def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))


def verify_identity(case,run,admission):
    """模块一：绑定实际任务、原面输入、具名气体和新声速来流。"""
    case=Path(case)
    if run.get('profile')!=PROFILE or admission.get('profile')!=PROFILE:raise ValueError('错误运行身份')
    if run.get('model')!=MODEL or admission.get('model')!=MODEL:raise ValueError('错误气体模型')
    if admission.get('case')!=case.relative_to(ROOT).as_posix():raise ValueError('源目录错误')
    process=run['process']
    if process.get('liveGroupAfterCleanup') or not process.get('inputStable') or not all(process['inputStable'].values()):raise ValueError('输入变化或进程未清点')
    if not process.get('logExportPassed') or not run.get('actualSolverDispatchVerified'):raise ValueError('日志/实际分派未验证')
    if sha(BASE/'rocket_run01/admission.json')!=run['admissionSHA256']:raise ValueError('准入回执SHA改变')
    hashes={str(local_path(p).resolve()):h for p,h in admission['inputSHA256'].items()}
    required=['exact_mesh.npz','wall_mapping.json','patches.json','thermo_contract.json',
        'constant/physicalProperties','system/controlDict','system/fvSchemes','system/fvSolution',
        'constant/polyMesh/points','constant/polyMesh/faces','constant/polyMesh/owner',
        'constant/polyMesh/neighbour','constant/polyMesh/boundary','0/p','0/T','0/U']
    for name in required:
        p=case/name
        if sha(p)!=hashes.get(str(p.resolve())):raise ValueError('实际输入SHA改变 '+name)
    model=read_json(case/'thermo_contract.json')
    if model!=contract():raise ValueError('气体数据/表示/参考焓改变')
    library=BASE/'native01/library/librahNistAir.so'
    if sha(library)!=admission['registrationLibrarySHA256'] or sha(library)!=hashes.get(str(library.resolve())):raise ValueError('实际注册库改变')
    condition=read_json(ROOT/'outputs/v194_validation/rocket_design03/manifest.json')['condition']
    if admission['condition']!=condition:raise ValueError('错误来流工况')
    u=np.array(condition['velocityDirection'])*condition['Mach']*float(evaluate(condition['pressurePa'],condition['temperatureK'])['soundSpeedMps'])
    if not np.allclose(u,admission['velocityMps'],rtol=0,atol=1e-10):raise ValueError('混用旧气体速度')
    return model


def read_frame(name,expected_time,diagnostic=False,run_override=None,admission_override=None):
    """模块二：完整原生场、边界与原面严格映射；异常仅显式诊断读取。"""
    if not re.fullmatch(r'[0-9.eE+-]+',name) or float(name)<=0 or float(name)!=expected_time:raise ValueError('初场/源时刻错误')
    case=BASE/'rocket01';run=run_override or read_json(BASE/'rocket_run01/report.json')
    admission=admission_override or read_json(BASE/'rocket_run01/admission.json')
    if not run.get('passed') and not diagnostic:raise ValueError('异常/未完成任务默认拒绝')
    verify_identity(case,run,admission)
    if name not in run['saved']:raise ValueError('没有监督绑定的完整保存场')
    saved=run['saved'][name]
    if saved['timeSeconds']!=expected_time or set(saved['fieldsSHA256'])!=set(FIELDS):raise ValueError('保存字段/时刻错误')
    folder=case/name
    for n in FIELDS:
        p=folder/n;data=p.read_bytes()
        if not data.rstrip().endswith(END_MARK) or sha(p)!=saved['fieldsSHA256'][n]:raise ValueError('字段未完成或SHA改变 '+n)
        head=data[:2400].decode('utf-8');location=re.search(r'\blocation\s+"([^"]+)";',head)
        if not location or float(location.group(1))!=expected_time:raise ValueError('字段location时刻错误')
        if not re.search(r'\bobject\s+'+re.escape(n)+r'\s*;',head):raise ValueError('字段名称错误')
    g,m,ids,owners,_=mapped_geometry(case);state=read_state(folder,len(g['cells']));state_audit=audit_state(state,evaluate)
    boundary={n:boundary_values(folder,n,state[n],g) for n in state};audit_state(boundary,evaluate);native_mass_flux(folder,g)
    rows=ids-len(g['neighbor'])
    for n in ['p','T']:
        if not np.allclose(boundary[n][rows],state[n][owners],rtol=2e-12,atol=1e-10):raise ValueError('壁面零梯度状态不符')
    p=state['p'][owners];T=state['T'][owners];rho=state['rho'][owners];u=state['U'][owners];es=state['e'][owners]
    gas=evaluate(p,T);h=es+p/rho;h0=h+np.einsum('ij,ij->i',u,u)/2
    error=float(np.max(abs(h-gas['enthalpyJkg'])/np.maximum(1,abs(gas['enthalpyJkg']))))
    if error>2e-8:raise ValueError('原生内能/新模型焓零点不相容')
    clock=folder/'uniform/time';text=clock.read_text(encoding='utf-8')
    index=int(re.search(r'\bindex\s+(\d+);',text).group(1));native_time=float(re.search(r'\bvalue\s+([\deE.+-]+);',text).group(1))
    if index<=0 or float(format(native_time,'.12g'))!=expected_time:raise ValueError('原生时间字典错误')
    fields=dict(pressurePa=p,fluidTemperatureK=T,densityKgm3=rho,velocityMps=u,nativeSensibleInternalEnergyJkg=es,
        enthalpyJkg=h,totalEnthalpyJkg=h0,entropyJkgK=gas['entropyJkgK'])
    meta=dict(schema='RAH.NISTOriginalWallResearch.v1',profile=PROFILE,model=MODEL,timeSeconds=expected_time,nativeTimeSeconds=native_time,timeIndex=index,
        timeDictionarySHA256=sha(clock),status='unqualified-development' if run['passed'] else 'interrupted-run-diagnostic',
        faceCount=19752,originalAreaM2=float(g['area'][ids].sum()),temperatureIdentity='wall-adjacent fluid; not Te/Tr/solid Tw',
        samplePositionKind='exact native-owner tetrahedron volume centroid',
        method='OpenFOAM 13 shockFluid / NIST integrated Cp / synchronized Euler / upwind / Kurganov',
        energyReference='hf=0; hs(298.15K)=0; es=hs-RT; h=es+p/rho; h0=h+|U|^2/2',
        fieldsSHA256=saved['fieldsSHA256'],runReportSHA256=sha(BASE/'rocket_run01/report.json'),admissionSHA256=run['admissionSHA256'],
        stateAudit=state_audit,enthalpyMaxScaled=error,attachedEdgeQualified=False,engineeringHeatReady=False,
        physicalQualification=False,wallHeatFluxWm2=None,solidWallTemperatureK=None,newQualifiedHeatAreaM2=0)
    return meta,fields,g,m,ids,owners
