"""模块职责：复用已编译记录器，唯一有界捕获原故障；不修改成熟数值设置。"""
from pathlib import Path
import json,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v231_validation';DOC=ROOT/'docs/validation/v231'
PREV=ROOT/'outputs/v230_validation';CASE=BASE/'capture01';ENGINE=PREV/'build02/rocket.exe'
CLOSED=False

def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False);CASE.mkdir(exist_ok=False)
    index=read(ROOT/'docs/validation/v230/evidence/local_asset_index.json')['files']
    assert sha(ENGINE)==index['build02/rocket.exe']['sha256']
    assert read(PREV/'process_build_observer/report.json')['passed']
    deps=read(PREV/'build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for p in ('rocketEulerRun.cpp','rocketEulerFailure.hpp'):
        assert sha(ROOT/'native/su2PhaseAwareAir'/p)==sha(PREV/'observer01'/p)
    src=(PREV/'rocket01/case.cfg').read_text(encoding='utf8')
    assert src.count('\nITER= 800\n')==1
    dst=src.replace('\nITER= 800\n','\nITER= 31\n')
    (CASE/'case.cfg').write_text(dst,encoding='utf8')
    shutil.copyfile(PREV/'rocket01/mesh.su2',CASE/'mesh.su2')
    source=PREV/'source01/SU2_CFD/src/iteration/CFluidIteration.cpp'
    assert sha(source)==read(PREV/'prepare02/sources.json')['SU2_CFD/src/iteration/CFluidIteration.cpp']['sha256']
    text=source.read_text(encoding='utf8');a=text.index('void CFluidIteration::Solve(');b=text.index('void CFluidIteration::',a+10)
    (out/'actual_solve.cpp.txt').write_text(text[a:b],encoding='utf8')
    assert 'Inner_Iter < nInner_Iter' in text[a:b]
    manifest={str(p):sha(p) for p in (ENGINE,CASE/'case.cfg',CASE/'mesh.su2',source,PREV/'build02/dependencies.json')}
    save(out/'report.json',dict(passed=True,changesOnly={'ITER':['800','31']},expectedMaxUpdates=31,
        actualSourceSHA256=sha(source),inputSHA256=manifest,dependencies=len(deps),CFDExecuted=False))
    save(CASE/'contract.json',dict(role='bounded-failure-capture',expectedOriginalCompletedUpdates=30,
      updatesLimit=31,wallSecondsLimit=60,expectedFailureExitCode=86,automaticRetries=0,
      rawStateRole='possibly partially updated; not accepted; no continuation',physicalQualification=False))

def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);run=read(BASE/'process_capture/report.json')
    assert run['inputsStable'] and not run['liveGroupAfterCleanup'] and run['stopReason'] is None
    assert run['exitCode']==86 and not run['passed'],run['exitCode']
    comparison={}
    for n in ('history.csv','physical_history.csv','state_000000.csv','state_000001.csv'):
        original=PREV/'rocket01'/n;current=CASE/n
        assert current.read_bytes()==original.read_bytes(),n
        comparison[n]=dict(sha256=sha(current),byteIdentical=True)
    termination={}
    for line in (CASE/'termination.txt').read_text(encoding='utf8').splitlines():
        k,v=line.split('=',1);termination[k]=v
    assert int(termination['last_observed_complete_update'])==30
    assert int(termination['native_inner_iteration'])==30
    raw=np.loadtxt(str(CASE/'termination_raw_state.csv'),delimiter=',',skiprows=1)
    assert raw.shape==(18027,6) and np.isfinite(raw).all()
    assert np.array_equal(np.sort(raw[:,0]),np.arange(18027))
    save(out/'report.json',dict(passed=True,originalFirst30TrajectoryByteIdentical=True,comparison=comparison,
      termination=termination,rawStateSHA256=sha(CASE/'termination_raw_state.csv'),rawStateAccepted=False,
      rawNodes=18027,normalCompletion=False,scientificAcceptance=False,newHeatAreaM2=0,newQualifiedAreaM2=0))
    print(json.dumps(dict(passed=True,trajectoryIdentical=True,termination=termination)))

def main(stage):
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',ENGINE,
        PREV/'build02/dependencies.json',PREV/'observer01/rocketEulerRun.cpp',PREV/'observer01/rocketEulerFailure.hpp',
        PREV/'rocket01/mesh.su2',PREV/'rocket01/case.cfg',PREV/'prepare02/sources.json']
    locks += [Path(n) for n in read(PREV/'build02/dependencies.json')]
    if stage=='prepare':cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];kind='offline';seconds=120;cwd=BASE
    elif stage=='capture':
        assert not (BASE/'process_capture').exists(),'Unique failure capture already consumed'
        report=read(BASE/'prepare01/report.json');assert report['passed']
        for n,h in report['inputSHA256'].items():assert sha(n)==h,n
        locks += [CASE/'case.cfg',CASE/'mesh.su2',CASE/'contract.json']
        cmd=[str(ENGINE),'case.cfg','.'];kind='cfd';seconds=60;cwd=CASE
    elif stage=='audit':
        locks += [p for p in CASE.iterdir() if p.is_file()]+[PREV/'rocket01/history.csv',PREV/'rocket01/physical_history.csv',
            PREV/'rocket01/state_000000.csv',PREV/'rocket01/state_000001.csv']
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_audit'];kind='offline';seconds=120;cwd=BASE
    else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise(stage,cmd,seconds,128,list(dict.fromkeys(locks)),cwd,kind,False)

if __name__=='__main__':
    if CLOSED:raise SystemExit('v231 failure capture closed; no retry or continuation')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_audit':audit()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
