"""模块职责：对新相态波系/原面研究连接执行有限监督；不启动CFD。"""
from pathlib import Path
import sys,shutil,time,json
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v225_validation';DOC=ROOT/'docs/validation/v225'
CLOSED=True
SOURCES=['phaseAwareAirCaloric.m','solvePhaseAwareWave.m','marchPhaseAwareWavePath.m','phaseAwareAirState.m',
 'phaseAwareAirReferenceState.m','phaseAwareMeshIdentity.m',
 'tests/testPhaseAwareWaveV225.m','data/thermo/phase_aware_dry_air_v1.json']

def inspect():
    from scipy.io import loadmat
    start=time.monotonic()
    d=loadmat(str(ROOT/'outputs/v104_validation/baseline02/current_case.mat'),squeeze_me=True,struct_as_record=False)
    m=d['current'];sc=m.surfaceCoord;fs=d['fs']
    report={'meshFields':m._fieldnames,'sectionFields':sc._fieldnames,'fs':{n:str(getattr(fs,n))[:150] for n in fs._fieldnames},
      'meshFaces':m.F.shape[0], 'sectionX':sc.sectionX.tolist(),'sectionRadius':sc.sectionR.tolist() if hasattr(sc,'sectionR') else None,
      'seconds':time.monotonic()-start}
    save(BASE/'inspect01.json',report);print(json.dumps(report))

def main(name):
    if CLOSED:raise SystemExit('v225 validation entry closed')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda unused:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if name=='inspect01':
        locks.append(ROOT/'outputs/v104_validation/baseline02/current_case.mat')
        return guard.supervise(name,[sys.executable,str(Path(__file__).resolve()),'task_inspect'],30,128,locks,BASE)
    if name=='audit01':
        return guard.supervise(name,[sys.executable,str(ROOT/'scripts/auditPhaseWaveV225.py')],120,128,
          locks+[ROOT/'scripts/auditPhaseWaveV225.py'],BASE,'audit')
    if name=='geometry01':
        out=BASE/name;out.mkdir(exist_ok=False)
        names=['inspectPhaseAwareRocketInlet.m','phaseAwareMeshIdentity.m','tests/completePhaseWaveGeometryV225.m']
        for rel in names:
            q=out/'executed_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,q)
        locks += [ROOT/n for n in names]+[BASE/'matlab02/waves.json',BASE/'audit01/rocket_inlet.json',
          ROOT/'outputs/v104_validation/baseline02/current_case.mat']
        expression="restoredefaultpath;maxNumCompThreads(1);addpath(pwd);addpath(fullfile(pwd,'tests'));completePhaseWaveGeometryV225('"+out.as_posix()+"');"
        return guard.supervise(name,['D:/matlab/bin/matlab.exe','-batch',expression],60,128,locks,ROOT,'matlab')
    assert name in ('matlab01','matlab02')
    out=BASE/name;out.mkdir(exist_ok=False)
    for rel in SOURCES:
        p=ROOT/rel;q=out/'executed_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    locks += [ROOT/n for n in SOURCES]+[ROOT/'outputs/v104_validation/baseline02/current_case.mat']
    seconds=240-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
    assert seconds>0
    if name=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
    expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseAwareWaveV225('"+out.as_posix()+"');"
    return guard.supervise(name,['D:/matlab/bin/matlab.exe','-batch',expression],min(120,seconds),128,locks,ROOT,'matlab')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v225 validation entry closed')
    if sys.argv[1]=='task_inspect':inspect()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
