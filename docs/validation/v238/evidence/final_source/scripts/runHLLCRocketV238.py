"""模块职责：登记唯一相容均匀初场HLLC整箭候选，复用成熟对象，监督固定收益与资源预算。"""
from pathlib import Path
import csv,difflib,json,os,shutil,subprocess,sys,time
import numpy as np
import runHLLCCompatibilityV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=run.BASE;DOC=run.DOC;STAGE=BASE/'rocket_prepare01';BUILD=BASE/'build04';CASE=BASE/'rocket01'
CPP=ROOT/'native/su2PhaseAwareAir/uniformHLLCBranchV238.cpp';HOOK=ROOT/'native/su2PhaseAwareAir/PhaseStepObservationV238.hpp'
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz';MESH=ROOT/'outputs/v222_validation/rocket01/mesh.su2'
CLOSED=True


def prepare():
    assert read(BASE/'native_audit01/report.json')['passed']
    for p in (STAGE,BUILD,CASE,BASE/'failure_control01'):p.mkdir(exist_ok=False)
    original=ROOT/'native/su2PhaseAwareAir/restartBranchV236.cpp'
    assert sha(original)==read(ROOT/'outputs/v236_validation/build01/inputs.json')[str(original)]
    common=original.read_text(encoding='utf8').split('class BranchDriverV236 final:')[0]
    (STAGE/'commonBranchV238.hpp').write_text(common,encoding='utf8')
    old=ROOT/'outputs/v237_validation/stage02/CEulerSolver.cpp';assert sha(old)==read(ROOT/'outputs/v237_validation/stage02/inputs.json')[str(old)]
    before=old.read_text(encoding='utf8');assert before.count('rah_phase_step_v237::native_limit(')==1
    text='#include "'+HOOK.as_posix()+'"\n'+before.replace('rah_phase_step_v237::native_limit(','rah_phase_step_v238::native_limit(')
    (STAGE/'CEulerSolver.cpp').write_text(text,encoding='utf8')
    (STAGE/'failure_observer.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=old.name,tofile='failure-observer/'+old.name)),encoding='utf8')
    cmd=read(ROOT/'outputs/v237_validation/stage02/commands.json')['euler'];cmd=list(cmd)
    cmd[cmd.index('-c')+1]=str(STAGE/'CEulerSolver.cpp');cmd[cmd.index('-o')+1]=str(BUILD/'CEulerSolver.obj');cmd[cmd.index('-MF')+1]=str(BUILD/'CEulerSolver.d')
    helper=read(BASE/'native_prepare01/commands.json')['helper']
    helper=[str(CPP) if x.endswith('nativeHLLCV238.cpp') else x.replace(str(BASE/'build03'),str(BUILD)) if x.endswith(('.a','.exe','.d')) else x for x in helper]+['-I'+str(STAGE)]
    commands=dict(euler=cmd,helper=helper);save(STAGE/'commands.json',commands)
    known=read(BASE/'native_prepare01/inputs.json');known.update(read(ROOT/'outputs/v237_validation/stage02/inputs.json'));tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');deps={}
    for name,cmd in commands.items():
        filtered=[];i=0
        while i<len(cmd):
            x=cmd[i]
            if x in ('-o','-MF'):i+=2;continue
            if x in ('-c','-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            filtered.append(x);i+=1
        dep=STAGE/(name+'.d');subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(run.OLDBUILD),timeout=20)
        for raw in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(raw);p=(p if p.is_absolute() else run.OLDBUILD/p).resolve();h=sha(p)
            if p not in (CPP,HOOK) and p.parent!=STAGE:
                expected=known.get(str(p),tree.get(str(p)))
                if expected is None:expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/p.relative_to(run.GAS/'source'))]
                assert h==expected,str(p)
            deps[str(p)]=h
    for rows in read(BASE/'build03/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'];deps[r['path']]=r['sha256']
    oldCfg=ROOT/'outputs/v237_validation/rocket01/case.cfg'
    config={k.strip():v.strip() for k,v in (l.split('=',1) for l in oldCfg.read_text(encoding='utf8').splitlines() if '=' in l)};prior=dict(config)
    config.update(RESTART_SOL='NO',CONV_NUM_METHOD_FLOW='HLLC',MESH_FILENAME=os.path.relpath(str(MESH),str(CASE)).replace('\\','/'))
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in config.items())+'\n',encoding='utf8')
    shutil.copyfile(ROOT/'outputs/v236_validation/rocket01/original_faces.txt',CASE/'original_faces.txt')
    frame=read_frame(BASE/'native01/state_000000.bin',24,True)['values'];six=frame[0,5:11].copy()
    assert np.array_equal(frame[:,5:11],np.tile(six,(24,1)))
    air=PhaseAwareViscousAir();p=float(config['FREESTREAM_PRESSURE']);t=float(config['FREESTREAM_TEMPERATURE']);q=air.evaluate(p,t)
    speed=float(config['MACH_NUMBER'])*float(q['soundSpeedMps']);angle=np.deg2rad(float(config['AOA']));v=np.array([speed*np.cos(angle),0.,speed*np.sin(angle)])
    rho=float(q['densityKgm3']);mu=float(air.interpolated_transport(p,t)['dynamicViscosityPas'])
    expected=np.r_[rho,rho*v,rho*(float(q['internalEnergyJkg'])+300000+.5*np.sum(v*v)),3*mu/rho]
    error=float(np.max(abs(six-expected)/np.maximum(1,abs(expected))));assert error<3e-12,error
    np.savetxt(str(CASE/'uniform_six_state.txt'),six[None,:],fmt='%.17g')
    contract=read(ROOT/'outputs/v236_validation/contract.json');contract.pop('sourceCheckpointSHA256')
    contract.update(branch='phase-air-general-hllc-uniform-v238',uniformFreestreamInitial=True,sourceObjectsRecompiledThisVersion=3,
        fixedNumericsPreservedExceptNamedHLLC=True,sourceSmallUniformFrameSHA256=sha(BASE/'native01/state_000000.bin'),
        uniformSixStateSHA256=sha(CASE/'uniform_six_state.txt'),uniformIndependentMaxScaledError=error,
        sourceCheckpoint=None,exactHistoricalContinuation=False,nativeCFLHistoryReset=False,
        originalCFLInitial=.1,configChanges={k:dict(old=prior[k],new=v) for k,v in config.items() if v!=prior[k]},
        actualLimitedStatesRecorded=True,nativeFailureStateRecorded=True,forcedFinalNativeRestartRequired=True)
    assert sha(MESH)==contract['meshSHA256'] and sha(GEOMETRY)==contract['geometrySHA256']
    save(BASE/'rocket_contract.json',contract)
    for p in (MESH,GEOMETRY,CASE/'case.cfg',CASE/'original_faces.txt',CASE/'uniform_six_state.txt',BASE/'rocket_contract.json',DOC/'ROCKET_PLAN.md',BASE/'native_audit01/report.json'):
        deps[str(p)]=sha(p)
    save(STAGE/'inputs.json',deps);save(STAGE/'report.json',dict(passed=True,uniformSix=six.tolist(),independentInitialError=error,realCFDAdmitted=False))


def build():
    commands=read(STAGE/'commands.json');t=time.monotonic();subprocess.check_call(commands['euler'],cwd=str(run.OLDBUILD));save(BUILD/'compile.json',dict(eulerSeconds=time.monotonic()-t))
    ar=str(Path(commands['euler'][0]).parent/'ar.exe');members={};changed=0
    for name,rows in read(BASE/'build03/members.json').items():
        result=[]
        for r in rows:
            p=Path(r['path']);assert sha(p)==r['sha256']
            if p.name=='CEulerSolver.obj':p=BUILD/'CEulerSolver.obj';changed+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result]);members[name]=result
    assert changed==1;save(BUILD/'members.json',members);subprocess.check_call(commands['helper'],cwd=str(run.OLDBUILD))
    save(BUILD/'report.json',dict(passed=True,driverSHA256=sha(BUILD/'driver.exe'),observerWrapperOnly=True))


def admit():
    assert read(BASE/'process_failure_control/report.json')['passed']
    with (BASE/'failure_control01/failed_step.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==1 and rows[0]['reason']=='PHASE_STEP_NO_ADMISSIBLE_PROGRESS' and int(rows[0]['localID'])==7 and int(rows[0]['update'])==3 and int(rows[0]['globalID'])==-1
    for prefix in ('old','delta'):assert np.isfinite([float(rows[0][prefix+str(j)]) for j in range(5)]).all()
    bound=read(STAGE/'inputs.json')
    for p,h in bound.items():assert sha(p)==h,p
    assert sha(BUILD/'driver.exe')==read(BUILD/'report.json')['driverSHA256']
    for name,rows in read(BUILD/'members.json').items():
        bound[str(BUILD/name)]=sha(BUILD/name)
        for r in rows:assert sha(r['path'])==r['sha256'];bound[r['path']]=r['sha256']
    for p in (BUILD/'driver.exe',Path(__file__).resolve()):bound[str(p)]=sha(p)
    save(BASE/'rocket_inputs.json',bound);save(BASE/'rocket_admission.json',dict(passed=True,unique=True,contractSHA256=sha(BASE/'rocket_contract.json'),
        inputsSHA256=sha(BASE/'rocket_inputs.json'),failureObserverControl=True,allNecessaryNativeControlsPassed=True,physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v238 rocket execution closed; follow NEXT_DECISION')
    run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'ROCKET_PLAN.md',CPP,HOOK]
    if stage=='prepare':return run.guard.supervise('rocket_prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,1024,locks,BASE)
    for p,h in read(STAGE/'inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
    locks+=[STAGE/'inputs.json',STAGE/'commands.json']
    if stage=='build':
        remaining=180-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        return run.guard.supervise('rocket_build',[sys.executable,str(Path(__file__).resolve()),'task_build'],remaining,1024,locks,run.OLDBUILD,'build')
    if stage=='failure_control':
        assert read(BUILD/'report.json')['passed'];locks.append(BUILD/'driver.exe')
        return run.guard.supervise('failure_control',[str(BUILD/'driver.exe'),'--failure-control','.'],10,1024,locks,BASE/'failure_control01','query')
    if stage=='admit':return run.guard.supervise('rocket_admit',[sys.executable,str(Path(__file__).resolve()),'task_admit'],30,1024,locks,BASE)
    if stage=='run':
        assert read(BASE/'rocket_admission.json')['passed'] and not (BASE/'process_rocket').exists()
        for p,h in read(BASE/'rocket_inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
        locks+=[BASE/'rocket_admission.json',BASE/'rocket_inputs.json']
        return run.guard.supervise('rocket',[str(BUILD/'driver.exe'),'case.cfg','.','uniform_six_state.txt','original_faces.txt'],2400,1024,list(dict.fromkeys(locks)),CASE,'cfd',True)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 rocket execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    elif sys.argv[1]=='task_admit':admit()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
