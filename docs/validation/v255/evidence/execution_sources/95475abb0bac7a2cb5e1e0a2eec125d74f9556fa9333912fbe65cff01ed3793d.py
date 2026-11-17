"""模块职责：固定本轮预算，全部源类型预先快照，检查真实子进程与终态。"""
import sys,shutil
from pathlib import Path
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v255_validation';DOC=ROOT/'docs/validation/v255';CLOSED=False
OWN=['solveRegularizedEntropyIntegral.m','phaseAwareAirCaloricDrop.m','tests/testRegularizedEntropyV255.m']
LIB=['phaseAwareAirState.m','phaseAwareAirCaloric.m','data/thermo/phase_aware_dry_air_v1.json']
INPUT=['outputs/v254_validation/prepare01/source.json','outputs/v252_validation/matlab01/surface_base.json',
    'outputs/v253_validation/inlet02/inlet_profiles.npz','outputs/v226_validation/matlab02/rocket_metadata.json']


def main(stage):
    if CLOSED:raise SystemExit('v255 fixed regularized tasks closed; follow NEXT_DECISION')
    assert stage in ('prepare','matlab01','matlab02','audit','export')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original=guard.size_bytes;guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py']+[ROOT/n for n in OWN+LIB+INPUT]
    if stage=='prepare':
        script=ROOT/'scripts/prepareRegularizedEntropyV255.py';locks+=[script];argv=[sys.executable,str(script)];seconds=60;kind='offline'
    elif stage.startswith('matlab'):
        locks += [BASE/'prepare01/source.json',BASE/'prepare01/geometry.json']
        if stage=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
        out=BASE/stage;out.mkdir(exist_ok=False)
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'));seconds=min(180,360-used);kind='matlab'
        expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testRegularizedEntropyV255('"+out.as_posix()+"');"
        argv=['D:/matlab/bin/matlab.exe','-batch',expression]
    else:
        script=ROOT/'scripts'/('auditRegularizedEntropyV255.py' if stage=='audit' else 'exportRegularizedEntropyV255.py')
        locks += [script,BASE/'prepare01/source.json',BASE/'prepare01/geometry.json',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py']
        locks += list((BASE/'matlab01').glob('*.json'))+list((BASE/'matlab02').glob('*.json'))
        if stage=='export':
            locks += [BASE/'audit01/report.json',ROOT/'web/regularized-entropy-v255.js',ROOT/'web/regularized-entropy-v255.html',
                ROOT/'tests/testRegularizedEntropyViewerV255.js',ROOT/'outputs/v252_validation/workbench02/case-data.js',
                ROOT/'outputs/v252_validation/workbench02/index.html',ROOT/'outputs/v254_validation/workbench04/coupled.json',
                ROOT/'web/section.js',ROOT/'web/viewer.js',ROOT/'web/export.js']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='audit');seconds=120-used;kind='audit';argv=[sys.executable,str(script)]
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ('.py','.m','.md','.js','.html','.css','.json') and p.stat().st_size<2*1024**2:
            target=snap/(sha(p)+p.suffix)
            if not target.exists():shutil.copyfile(p,target)
    return supervise(guard,BASE,stage,argv,seconds,128,locks,ROOT,kind)


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
