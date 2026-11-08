"""只读核验真实能量拒绝：不将反算输入、完成日志或初场冒充已解流场。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from phaseAwareDryAir import PhaseAwareDryAir,MODEL,DOMAIN,pressure_ceiling
from readOpenFoamRocketState import mapped_geometry
from runOpenFoamRocketDevelopment import local_path
from auditOpenFoamDevelopment import log_progress

BASE=ROOT/'outputs/v206_validation'

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def parse_energy_failure(text):
    """模块一：保留原生17位输入与错误分支；不向表外反解温度。"""
    number=r'([\deE.+-]+)'
    pattern=r'RAH_PHASE_ENERGY_OUTSIDE_TABLE e='+number+r' elo='+number+r' ehi='+number+r' pOld='+number+r' seed='+number
    matches=re.findall(pattern,text)
    if len(matches)!=1:raise ValueError('缺少或重复能量越界记录')
    e,lo,hi,p,seed=map(float,matches[0])
    if not np.isfinite([e,lo,hi,p,seed]).all() or not e<lo<hi or p<=0:raise ValueError('错误分支与数值不符')
    return dict(internalEnergyJkg=e,tableLowerEnergyJkg=lo,tableUpperEnergyJkg=hi,oldPressurePa=p,seedTemperatureK=seed)

def verify_identity(case,run,admission):
    """模块二：实际类型、工况、准入、库/气体/输入逐项绑定；无完整场即无正向导出。"""
    profile='of13-phase-aware-air-original-rocket-sync-v1'
    if run.get('profile')!=profile or admission.get('profile')!=profile:raise ValueError('错误任务身份')
    if run.get('model')!=MODEL or admission.get('model')!=MODEL:raise ValueError('错误气体身份')
    if run.get('domain')!=DOMAIN or admission.get('domain')!=DOMAIN:raise ValueError('错误联合域')
    if admission.get('case')!=case.relative_to(ROOT).as_posix():raise ValueError('错误源目录')
    if not run.get('actualSolverDispatchVerified') or admission.get('energyType')!='rahPhaseAwareInternalEnergy':raise ValueError('实际能量类型未验证')
    process=run['process']
    if process.get('liveGroupAfterCleanup') or not process.get('inputStable') or not all(process['inputStable'].values()):raise ValueError('输入/进程清点未过')
    if not process.get('logExportPassed'):raise ValueError('日志未导出')
    if sha(BASE/'rocket_run01/admission.json')!=run['admissionSHA256']:raise ValueError('准入SHA改变')
    gas=PhaseAwareDryAir()
    if admission['contract']!=gas.contract() or read_json(case/'thermo_contract.json')!=gas.contract():raise ValueError('物性合同改变')
    for name,digest in admission['inputSHA256'].items():
        p=local_path(name)
        # CLI事后关闭；执行前副本必须与原准入SHA逐字一致。
        if p==ROOT/'scripts/runPhaseAwareAirRocket.py':p=BASE/'rocket_run01/runPhaseAwareAirRocket.py'
        if sha(p)!=digest:raise ValueError('准确输入或执行源改变 '+name)
    library=ROOT/'outputs/v206_validation/native01/library/librahPhaseAwareAir.so'
    if sha(library)!=admission['registrationLibrarySHA256']:raise ValueError('注册库改变')
    condition=read_json(ROOT/'outputs/v194_validation/rocket_design03/manifest.json')['condition']
    if admission['condition']!=condition:raise ValueError('实际工况改变')
    velocity=np.array(condition['velocityDirection'])*condition['Mach']*float(gas.evaluate(condition['pressurePa'],condition['temperatureK'])['soundSpeedMps'])
    if not np.allclose(velocity,admission['velocityMps'],rtol=0,atol=1e-10):raise ValueError('混用旧气体速度')
    return gas

def require_developed_frame(run,name):
    """模块三：异常任务/初态/没有保存证据均拒绝作为工程外缘。"""
    if not run.get('passed'):raise ValueError('未正常完成的真实候选不得导出工程外缘')
    if name not in run.get('saved',{}) or float(name)<=0:raise ValueError('没有完整发展帧')
    if not run.get('attachedEdgeQualified'):raise ValueError('尚未验证为边界层外缘')

def main():
    start=time.monotonic();out=BASE/'rocket_audit01';out.mkdir(exist_ok=False);case=BASE/'rocket01'
    for n in ['auditPhaseAwareAirRocket.py','phaseAwareDryAir.py','readOpenFoamRocketState.py','auditOpenFoamDevelopment.py']:
        shutil.copyfile(ROOT/'scripts'/n,out/n)
    run=read_json(BASE/'rocket_run01/report.json');a=read_json(BASE/'rocket_run01/admission.json')
    gas=verify_identity(case,run,a);process=run['process']
    if process['nativeExitCode']!=1 or process['stopReason']!='scientific-guard-after-exit':raise ValueError('不符已观测终态')
    text=(BASE/'rocket_run01/native/stdout.txt').read_text(encoding='utf-8')
    rejected=parse_energy_failure(text);progress=log_progress(text)
    p=rejected['oldPressurePa'];endpoint=gas.thermal(np.array([60.,2000.]))[1]-gas.R*np.array([60.,2000.])
    endpoint_error=float(np.max(abs(endpoint-[rejected['tableLowerEnergyJkg'],rejected['tableUpperEnergyJkg']])/np.maximum(1,abs(endpoint))))
    if endpoint_error>2e-10:raise ValueError('原生能量端点与新模型不一致')
    minimum_cv=float(gas.Cp.min()-gas.R)
    if minimum_cv<=0:raise ValueError('无法证明表内e严格单调')
    # 分段线性Cp的最小值在节点，Cv>0严格证明无表内根；不计算表外温度。
    deficit=rejected['tableLowerEnergyJkg']-rejected['internalEnergyJkg']
    g,m,ids,owners,surface=mapped_geometry(case)
    directories=[]
    for d in case.iterdir():
        if not d.is_dir():continue
        try:t=float(d.name)
        except ValueError:continue
        directories.append(dict(name=d.name,timeSeconds=t,files=sorted(x.name for x in d.rglob('*') if x.is_file())))
    if run['saved'] or any(d['timeSeconds']>0 for d in directories):raise ValueError('存在需另审核的保存场')
    co_error=abs(run['initialMaxCo']/a['initialMaxCoExpected']-1)
    if co_error>5e-6 or progress['announcedSteps']!=progress['completedLogSteps']+1:raise ValueError('首Co或完整步数不符')
    protected={k:sha(ROOT/k.replace('\\','/'))==v for k,v in read_json(ROOT/'outputs/v176_validation/protected_initial.json').items()}
    if not all(protected.values()):raise ValueError('并发保护文件改变')
    result=dict(passed=True,candidatePassed=False,candidateClosed=True,profile=run['profile'],model=MODEL,domain=DOMAIN,
        actualSolverDispatchVerified=True,nativeFatalPrecededSupervisorObservation=True,supervisorTerminatedRunningSolver=False,
        progress=progress,rejectedInput=rejected,energyDeficitJkg=deficit,energyEndpointMaxScaledError=endpoint_error,
        minimumTableCvJkgK=minimum_cv,noRootWithinDeclaredTemperatureDomain=True,
        oldPressureWithinDomainAt60K=bool(p<=pressure_ceiling(60)),pressureCeilingAt60KPa=float(pressure_ceiling(60)),
        actualInitialGuardVerified=run['actualInitialStateGuardVerified'],actualPerStepGuardVerified=run['actualPerStepGuardVerified'],
        finalPressureUnknown=True,seedIsNotAcceptedGlobalMinimum=True,
        temperatureRootComputedOutsideTable=False,rejectedInputIsNotAcceptedPhysicalState=True,
        firstIllegalNativeCellUnknown=True,uniqueFailureRootCauseUnknown=True,fullFieldMinimumTemperatureK=None,
        initialCoRelativeError=co_error,initialDeltaTSeconds=float(re.search(r'deltaT = ([\deE.+-]+)',text).group(1)),
        originalWallFaces=len(ids),originalAreaM2=float(g['area'][ids].sum()),originalWallExact=True,
        savedDirectories=directories,completeSavedDevelopmentFields=0,positiveDevelopedReaderTested=False,
        processSeconds=process['elapsedSeconds'],entrySeconds=run['elapsedSeconds'],parentTimeoutSeconds=960,parentElapsedSeconds=None,
        peakRSSBytes=process['peakObservedGroupRSSBytes'],stableInputCount=len(process['inputStable']),
        protected=protected,physicalQualification=False,numericallyConverged=False,attachedEdgeQualified=False,
        wallHeatFluxWm2=None,solidWallTemperatureK=None,newQualifiedHeatAreaM2=0,
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [BASE/'rocket_run01/report.json',BASE/'rocket_run01/admission.json',
            BASE/'rocket_run01/native/stdout.txt',Path(__file__).resolve(),ROOT/'scripts/phaseAwareDryAir.py']},elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    print(json.dumps({k:result[k] for k in ['passed','progress','rejectedInput','energyDeficitJkg','energyEndpointMaxScaledError','initialCoRelativeError','elapsedSeconds']}))

if __name__=='__main__':main()
