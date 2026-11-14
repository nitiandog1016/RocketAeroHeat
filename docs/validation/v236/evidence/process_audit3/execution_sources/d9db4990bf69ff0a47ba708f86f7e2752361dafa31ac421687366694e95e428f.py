"""模块职责：固定一次真实状态重启动的身份、资源与收益合同，监督全过程，绝不自动重试CFD。"""
from pathlib import Path
import os,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v236_validation';DOC=ROOT/'docs/validation/v236';CASE=BASE/'rocket01';BUILD=BASE/'build01'
OLD=ROOT/'outputs/v235_validation';ORIGIN=ROOT/'outputs/v222_validation/rocket01'
CPP=ROOT/'native/su2PhaseAwareAir/restartBranchV236.cpp';CLOSED=False
CHECKPOINT=ORIGIN/'restart_flow_000199.csv';MESH=ORIGIN/'mesh.su2'
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def configure_guard():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original=guard.size_bytes
    # 整版只计一次；旧监督器audit分支还会补加audit目录，此处这些补项返回零。
    guard.size_bytes=lambda p:original(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0


def prepare():
    """模块一：复用旧网格/场，不复制大资产；产生原面索引、机器合同及薄驱动依赖。"""
    assert read(OLD/'audit01/report.json')['passed']
    expected=read(ROOT/'outputs/v222_validation/audit01/report.json')['outputSHA256']['restart_flow_000199.csv']
    assert sha(CHECKPOINT)==expected
    gm=read(ROOT/'outputs/v222_validation/geometry01/report.json')
    assert sha(MESH)==gm['meshSHA256'] and sha(GEOMETRY)==gm['cacheSHA256']
    CASE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    cfg={k.strip():v.strip() for k,v in (line.split('=',1) for line in (ORIGIN/'case.cfg').read_text(encoding='utf8').splitlines() if '=' in line)}
    oldcfg=dict(cfg)
    cfg.update(RESTART_SOL='YES',MESH_FILENAME=os.path.relpath(str(MESH),str(CASE)).replace('\\','/'),
        SOLUTION_FILENAME=os.path.relpath(str(CHECKPOINT.with_suffix('')),str(CASE)).replace('\\','/'),
        ITER='400',OUTPUT_FILES='(RESTART_ASCII)',OUTPUT_WRT_FREQ='400',WRT_RESTART_OVERWRITE='YES')
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    g=np.load(str(GEOMETRY));a=np.column_stack([np.arange(1,19753),g['triangles'],g['faceArea']])
    np.savetxt(str(CASE/'original_faces.txt'),a,fmt=['%d']*4+['%.17g'])
    commands=read(OLD/'stage01/commands.json');cmd=commands['helper']
    cmd=[str(CPP) if x.endswith('singleDomainRunV235.cpp') else str(BUILD/'branch.d') if x.endswith('control.d') else
         str(BUILD/'branch.exe') if x.endswith('control.exe') else x for x in cmd]
    save(BUILD/'command.json',cmd)
    filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    cwd=ROOT/'tools/runtime/su2-thermal-build/build-v117'
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(BUILD/'inputs.d')],cwd=str(cwd),timeout=20)
    known=read(OLD/'stage01/dependencies.json');bound={}
    for name in (BUILD/'inputs.d').read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else cwd/p).resolve();h=sha(p)
        if p!=CPP:assert h==known[str(p)],str(p)
        bound[str(p)]=h
    for rows in read(OLD/'build01/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256'];bound[row['path']]=row['sha256']
    for p in (CHECKPOINT,MESH,GEOMETRY,CASE/'case.cfg',CASE/'original_faces.txt',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
              OLD/'build01/libSU2core.a',OLD/'build01/libSU2Common.a',OLD/'source01/bindings.json'):
        bound[str(p)]=sha(p)
    save(BUILD/'inputs.json',bound)
    contract=dict(branch='v222-state199-native-reset-v236',maximumUpdates=400,cfdSeconds=2400,buildSeconds=60,
        cpuLimit=1,sampledGroupMemoryBytes=4*1024**3,wholeVersionOutputBytes=1024**3,
        sourceCheckpointSHA256=expected,meshSHA256=sha(MESH),geometrySHA256=sha(GEOMETRY),
        originalFaces=19752,originalAreaM2=float(sum(g['faceArea'])),points=335726,cells=705490,
        exactHistoricalContinuation=False,initialSixStatesAndCoordinatesBitCheck=True,nativeCFLHistoryReset=True,
        initialRHSUnknown=True,physicalEnergyOffsetJkg=300000,physicalQualification=False,
        fixedNumericsPreserved=True,configChanges={k:dict(old=oldcfg[k],new=v) for k,v in cfg.items() if v!=oldcfg[k]},
        adjustmentUpdates=100,windowLength=50,checkUpdates=[200,300],significantWorseningRatio=1.25,
        geometricMeanImprovementRatio=.95,pressureChangeScreen=.01,heatChangeScreen=.05,
        consecutiveNoBenefitWindowsToStop=2,qualificationGatesUnchanged=True,automaticRetries=0,
        initialAndCompleteSaveUpdates=[0,1,100,200,300,400],wallSaveInterval=50,
        heatDiagnostic='current primitive k*(TnormalNeighbor-Tw)/normal projected distance; native cache separate',
        sourceObjectsRecompiled=0)
    save(BASE/'contract.json',contract)
    save(BASE/'prepare_report.json',dict(passed=True,boundInputs=len(bound),nativeCalls=0,realUpdates=0,contractSHA256=sha(BASE/'contract.json')))


def main(stage):
    if CLOSED:raise SystemExit('v236 branch closed; no automatic continuation')
    BASE.mkdir(exist_ok=True);configure_guard()
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',CPP,ROOT/'scripts/runSU2RocketPhaseV222.py']
    if stage=='prepare':
        locks += [ORIGIN/'case.cfg',OLD/'stage01/commands.json',OLD/'stage01/dependencies.json',OLD/'build01/members.json',CHECKPOINT,MESH,GEOMETRY]
        return guard.supervise('prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,1024,locks,BASE)
    assert read(BASE/'prepare_report.json')['passed'];bound=read(BUILD/'inputs.json')
    for n,h in bound.items():assert sha(n)==h,n
    locks += [Path(n) for n in bound]+[BUILD/'inputs.json',BUILD/'command.json',BASE/'contract.json']
    locks=list(dict.fromkeys(locks))
    if stage=='build':return guard.supervise('build',read(BUILD/'command.json'),60,1024,locks,ROOT/'tools/runtime/su2-thermal-build/build-v117','build')
    if stage=='run':
        assert read(BASE/'process_build/report.json')['passed']
        assert not (BASE/'process_rocket').exists(),'Unique actual branch already consumed'
        locks.append(BUILD/'branch.exe')
        return guard.supervise('rocket',[str(BUILD/'branch.exe'),'case.cfg','.',os.path.relpath(str(CHECKPOINT),str(CASE)),
             'original_faces.txt'],2400,1024,locks,CASE,'cfd',True)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v236 branch closed; no automatic continuation')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
