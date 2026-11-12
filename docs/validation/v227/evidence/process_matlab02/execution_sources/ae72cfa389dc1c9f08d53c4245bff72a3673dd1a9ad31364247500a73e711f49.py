"""模块职责：绑定唯一原面热剖面任务及独立审核，预算覆盖整版输出。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v227_validation';DOC=ROOT/'docs/validation/v227';CLOSED=False
SOURCES=['solvePhaseAwareLaminarProfile.m','solvePhaseAwareNoseHeating.m','runPhaseAwareNoseHeatingCase.m',
 'buildPhaseAwareNoseHeatingViewerData.m','runAeroHeatSingleCase.m','buildAerothermalViewerData.m',
 'buildPhaseAwareNoseViewerData.m','phaseAwarePayloadSHA.m','solvePhaseAwareNoseEdge.m',
 'runPhaseAwareNoseCase.m','phaseAwareAirState.m','phaseAwareAirCaloric.m','phaseAwareMeshIdentity.m',
 'solvePhaseAwareWave.m','estimateNoseRadius.m','surfaceWallTemperature.m','calcFreestream.m',
 'tests/testPhaseNoseHeatingV227.m','data/thermo/phase_aware_dry_air_v1.json']
def main(stage):
    if CLOSED:raise SystemExit('v227 execution closed')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda unused:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    old=ROOT/'outputs/v226_validation/matlab02'
    locks += [old/'rocket_faces.csv',old/'rocket_metadata.json']
    if stage in ('audit01','audit02','render01'):
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='audit')
        assert used<120
        return guard.supervise(stage,[sys.executable,str(ROOT/'scripts/auditPhaseNoseHeatingV227.py'),stage],120-used,64,
          locks+[ROOT/'scripts/auditPhaseNoseHeatingV227.py',ROOT/'tests/testPhaseNoseHeatingViewerV227.js',
            ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'web/viewer.js',
            BASE/'matlab02/heating.json',BASE/'matlab02/rocket_workbench.json'],BASE,'audit')
    assert stage in ('matlab01','matlab02');out=BASE/stage;out.mkdir(exist_ok=False)
    if stage=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
    for rel in SOURCES:
        q=out/'executed_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,q)
    locks += [ROOT/n for n in SOURCES]+[ROOT/'outputs/v104_validation/baseline02/current_case.mat']
    used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
    expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseNoseHeatingV227('"+out.as_posix()+"');"
    return guard.supervise(stage,['D:/matlab/bin/matlab.exe','-batch',expression],min(180,360-used),128,locks,ROOT,'matlab')
if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
