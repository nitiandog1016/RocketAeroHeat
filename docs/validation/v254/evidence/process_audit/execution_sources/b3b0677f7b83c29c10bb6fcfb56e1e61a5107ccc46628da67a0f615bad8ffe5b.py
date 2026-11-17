"""模块职责：固定联合工程计算与独立审核预算，准确保存执行源码和终态。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v254_validation';DOC=ROOT/'docs/validation/v254';CLOSED=False
OWN=['solvePhaseAwareEntropyIntegral.m','tests/testCoupledEntropyV254.m']
LIB=['phaseAwareAirState.m','phaseAwareAirCaloric.m','data/thermo/phase_aware_dry_air_v1.json']
INPUT=[ROOT/'outputs/v253_validation/inlet01/shock_mass_source.json',ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz',
       ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json']

def main(stage):
    if CLOSED:raise SystemExit('v254 fixed coupled integral tasks closed; follow NEXT_DECISION')
    assert stage in ('prepare','matlab01','matlab02','audit','export')
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original=guard.size_bytes;guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py']+INPUT+[ROOT/n for n in OWN+LIB]
    if stage=='prepare':
        script=ROOT/'scripts/prepareCoupledEntropyV254.py';locks += [script,ROOT/'scripts/phaseAwareEntropyInlet.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py']
        return supervise(guard,BASE,stage,[sys.executable,str(script)],120,128,locks,ROOT,'offline')
    locks += [BASE/'prepare01/source.json']
    if stage.startswith('matlab'):
        if stage=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
        out=BASE/stage;out.mkdir(exist_ok=False)
        for name in OWN+LIB:
            p=out/'executed_source'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,p)
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
        expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testCoupledEntropyV254('"+out.as_posix()+"');"
        return supervise(guard,BASE,stage,['D:/matlab/bin/matlab.exe','-batch',expression],min(300,600-used),128,locks,ROOT,'matlab')
    script=ROOT/'scripts'/('auditCoupledEntropyV254.py' if stage=='audit' else 'exportCoupledEntropyV254.py')
    used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='audit')
    locks += [script,ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py']+list((BASE/'matlab01').glob('*.json'))+list((BASE/'matlab02').glob('*.json'))
    if stage=='export':
        locks += [BASE/'audit01/report.json',ROOT/'web/coupled-entropy-v254.html',ROOT/'web/coupled-entropy-v254.js',
            ROOT/'tests/testCoupledEntropyViewerV254.js',ROOT/'outputs/v252_validation/workbench02/case-data.js',
            ROOT/'outputs/v252_validation/workbench02/index.html',ROOT/'outputs/v252_validation/matlab01/surface_base.json',
            ROOT/'outputs/v252_validation/matlab01/base.json',ROOT/'web/viewer.js',ROOT/'web/section.js',ROOT/'web/export.js']
    return supervise(guard,BASE,stage,[sys.executable,str(script)],min(120,180-used),128,locks,ROOT,'audit')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
