"""模块职责：登记、准备并监督唯一完整箭体候选；保留失败，不自动续算。"""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,sys,time
import psutil

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v222_validation'
CASE=BASE/'rocket01';DOC=ROOT/'docs/validation/v222'
ENGINE=ROOT/'outputs/v220_validation/build02/SU2_CFD.exe'
MESH=ROOT/'outputs/v167_validation/real01/mesh.su2'
EXECUTION_CLOSED=False

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def save(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8');os.replace(str(tmp),str(p))
def size_bytes(folder):return sum(p.stat().st_size for p in folder.rglob('*') if p.is_file())

def supervise(stage,argv,seconds,output_mib,locks,cwd,kind='offline',scientific=False):
    """模块一：采样整个子进程组；准确源/计划先快照，停止后核验输入和清点。"""
    folder=BASE/('process_'+stage);folder.mkdir(parents=True,exist_ok=False)
    if kind in ('offline','audit'):
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']==kind)
        seconds=min(seconds,{'offline':120,'audit':240}[kind]-used)
        assert seconds>0
    before={str(p):sha(p) for p in locks}
    snapshots=folder/'execution_sources';snapshots.mkdir()
    for p in locks:
        if p.suffix in ('.py','.md') and (p.parent==ROOT/'scripts' or p==DOC/'PLAN.md'):
            shutil.copyfile(p,snapshots/(before[str(p)]+p.suffix))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf8')
    save(folder/'admission.json',dict(argv=argv,limitSeconds=seconds,memoryBytes=4*1024**3,outputMiB=output_mib,
        inputSHA256=before,cwd=str(cwd),budgetKind=kind,cpuLimit=1,hardAddressSpaceLimit=False))
    start=time.monotonic();peak=0;virtual=0;family={};stop=None;error=None;seen='';last_scan=0;outsize=0
    with (folder/'stdout.txt').open('wb') as log:
        process=subprocess.Popen(argv,cwd=str(cwd),env=env,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        parent=psutil.Process(process.pid);parent.cpu_affinity(parent.cpu_affinity()[:1]);family[parent.pid]=parent
        save(folder/'live.json',dict(pid=parent.pid,startedUTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        try:
            while process.poll() is None:
                try:
                    for child in parent.children(recursive=True):family[child.pid]=child
                except psutil.NoSuchProcess:pass
                rss=0;vms=0
                for child in family.values():
                    try:mi=child.memory_info();rss+=mi.rss;vms+=mi.vms
                    except psutil.NoSuchProcess:pass
                peak=max(peak,rss);virtual=max(virtual,vms);elapsed=time.monotonic()-start
                if elapsed-last_scan>=1:
                    outsize=size_bytes(CASE if kind in ('cfd','parent') else folder)
                    if kind=='audit':outsize+=sum(size_bytes(p) for p in BASE.glob('audit*') if p.is_dir())
                    if scientific:
                        with (folder/'stdout.txt').open('rb') as tail:
                            tail.seek(max(0,tail.seek(0,2)-8192));seen=tail.read().decode('utf8',errors='replace')
                        error=bool(re.search(r'what\(\):|Error in|terminate called|Floating point exception|SIGSEGV|Divergence detected',seen,re.I))
                    last_scan=elapsed
                if elapsed>seconds:stop='time';break
                if max(rss,vms)>4*1024**3:stop='memory';break
                if outsize>output_mib*1024**2:stop='output';break
                if error:stop='scientific-log';break
                time.sleep(.1)
        finally:
            # 只终止本次登记的进程组，不影响用户或其他工作区进程。
            if process.poll() is None:
                for child in reversed(list(family.values())):
                    try:child.kill()
                    except psutil.NoSuchProcess:pass
            code=process.wait(timeout=5);_,alive=psutil.wait_procs(list(family.values()),timeout=5)
    stable=all(sha(Path(n))==h for n,h in before.items())
    report=dict(passed=code==0 and stop is None and stable and not alive,status='terminal',budgetKind=kind,argv=argv,
        exitCode=code,elapsedSeconds=time.monotonic()-start,limitSeconds=seconds,stopReason=stop,
        peakObservedGroupRSSBytes=peak,peakObservedGroupVirtualBytes=virtual,inputSHA256=before,inputsStable=stable,
        liveGroupAfterCleanup=[p.pid for p in alive],cpuLimit=1,outputMiB=output_mib,hardAddressSpaceLimit=False)
    save(folder/'report.json',report);print(json.dumps({k:report[k] for k in ('passed','exitCode','elapsedSeconds','stopReason','peakObservedGroupRSSBytes')}))
    return report['passed']

def prepare():
    """模块二：只复制原输入，显式冻结成熟组合及准确构建依赖；不运行引擎。"""
    assert sha(MESH)=='00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181'
    assert sha(ENGINE)=='badd860ac7a1399c06f0d411fcb13bd85c014da69c61a4b5288f28ff42b0f828'
    assert read(ROOT/'outputs/v221_validation/final01/report.json')['passed']
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    members=read(ROOT/'outputs/v220_validation/build02/members.json')
    for rows in members.values():
        for row in rows:assert sha(row['path'])==row['sha256']
    CASE.mkdir(exist_ok=False);shutil.copyfile(MESH,CASE/'mesh.su2')
    cfg={}
    for line in (MESH.parent/'case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    for key in ('SST_OPTIONS','FREESTREAM_TURBULENCEINTENSITY','FREESTREAM_TURB2LAMVISCRATIO','MU_REF','MU_T_REF','SUTHERLAND_CONSTANT','PRANDTL_LAM'):cfg.pop(key,None)
    cfg.update(KIND_TURB_MODEL='SA',SA_OPTIONS='(NEGATIVE, WITHFT2)',FREESTREAM_NU_FACTOR='3',RESTART_SOL='NO',
        FLUID_MODEL='RAH_PHASE_AIR',VISCOSITY_MODEL='RAH_PHASE_LINEAR',CONDUCTIVITY_MODEL='RAH_PHASE_LINEAR',
        KIND_TRANS_MODEL='NONE',USE_ACCURATE_TURB_JACOBIANS='NO',NUM_METHOD_GRAD_RECON='WEIGHTED_LEAST_SQUARES',
        ROE_KAPPA='0.5',ENTROPY_FIX_COEFF='0.001',LOW_MACH_PREC='NO',LOW_MACH_CORR='NO',
        ITER='200',CONV_FIELD='(RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_MOMENTUM-Z, RMS_ENERGY, RMS_NU_TILDE)',
        OUTPUT_WRT_FREQ='50',SCREEN_WRT_FREQ_INNER='10',VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)',
        SOLUTION_FILENAME='restart_flow')
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    # 只绑定输出阶段所需准确源；不重新构建或重复物性验证。
    sources=['Common/src/CConfig.cpp','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
        'SU2_CFD/src/output/CFlowCompOutput.cpp','SU2_CFD/src/iteration/CFluidIteration.cpp',
        'SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/solvers/CNSSolver.cpp']
    origin=ROOT/'outputs/v220_validation/source';bound={}
    for name in sources:
        p=origin/name;target=BASE/'source01'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target);bound[name]=sha(p)
    inputs=[MESH,ENGINE,CASE/'mesh.su2',CASE/'case.cfg',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
        ROOT/'outputs/v220_validation/build02/dependencies.json',ROOT/'outputs/v220_validation/build02/members.json',
        ROOT/'outputs/v221_validation/final01/report.json']
    save(CASE/'case_manifest.json',dict(model='rah-su2-phase-air-saneg-eplus300k-v1',inputSHA256={str(p):sha(p) for p in inputs},
        sourceSHA256=bound,dependencies=len(deps),cells=705490,points=335726,wallFaces=19752,wallAreaM2=1.9646292565715553,
        updates=200,expectedSavedIterations=[0,50,100,150,199],initialization='new-compatible-uniform-state',
        physicalQualification=False,newQualifiedAreaM2=0,automaticRetries=0))
    save(BASE/'prepare_report.json',dict(passed=True,dependencies=len(deps),meshSHA256=sha(MESH),configSHA256=sha(CASE/'case.cfg')))

def main(stage):
    if EXECUTION_CLOSED:raise SystemExit('v222 candidate closed; use NEXT_DECISION, no automatic continuation')
    BASE.mkdir(exist_ok=True)
    locks=[Path(__file__).resolve(),DOC/'PLAN.md']
    if stage=='prepare':return supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,256,locks,BASE)
    assert read(BASE/'prepare_report.json')['passed']
    manifest=read(CASE/'case_manifest.json')
    for n,h in manifest['inputSHA256'].items():assert sha(n)==h,n
    locks += [Path(n) for n in manifest['inputSHA256']]+[CASE/'case_manifest.json']
    if stage=='run':
        assert not (BASE/'process_engine').exists(),'Unique candidate already consumed'
        return supervise('parent',[sys.executable,str(Path(__file__).resolve()),'execute'],1020,2048,locks,CASE,'parent')
    if stage=='execute':return supervise('engine',[str(ENGINE),'case.cfg'],900,2048,locks,CASE,'cfd',True)
    raise ValueError('STAGE_NOT_ADMITTED')

if __name__=='__main__':
    if EXECUTION_CLOSED:raise SystemExit('v222 candidate closed; use NEXT_DECISION, no automatic continuation')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
