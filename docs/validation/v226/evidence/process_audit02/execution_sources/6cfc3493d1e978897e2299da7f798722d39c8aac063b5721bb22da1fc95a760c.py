"""模块职责：事前绑定相态鼻部生产器/原面/工作台代码并监督有限验证。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v226_validation';DOC=ROOT/'docs/validation/v226';CLOSED=False
SOURCES=['phaseAwarePayloadSHA.m','solvePhaseAwareNoseEdge.m','buildPhaseAwareNoseViewerData.m',
 'runPhaseAwareNoseCase.m','runAeroHeatSingleCase.m','buildAerothermalViewerData.m',
 'phaseAwareAirState.m','phaseAwareAirCaloric.m','phaseAwareMeshIdentity.m','solvePhaseAwareWave.m',
 'estimateNoseRadius.m','surfaceWallTemperature.m','calcFreestream.m','tests/testPhaseAwareNoseV226.m',
 'buildCanonicalBody.m','calcFaceGeometry.m','fixNormalDirection.m','buildSurfaceCoordinates.m','analyzeLocalGeometry.m',
 'data/thermo/phase_aware_dry_air_v1.json']
def main(stage):
    if CLOSED:raise SystemExit('v226 execution closed')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda unused:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage in ('audit01','audit02'):
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_audit*/report.json'))
        return guard.supervise(stage,[sys.executable,str(ROOT/'scripts/auditPhaseNoseV226.py'),stage],120-used,128,
          locks+[ROOT/'scripts/auditPhaseNoseV226.py',ROOT/'tests/testPhaseNoseViewerV226.js',
            ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'web/viewer.js',
            BASE/'matlab02/rocket_faces.csv',BASE/'matlab02/rocket_metadata.json',BASE/'matlab02/rocket_workbench.json'],BASE,'audit')
    assert stage in ('matlab01','matlab02');out=BASE/stage;out.mkdir(exist_ok=False)
    for rel in SOURCES:
        q=out/'executed_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,q)
    locks += [ROOT/n for n in SOURCES]+[ROOT/'outputs/v104_validation/baseline02/current_case.mat']
    if stage=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
    remaining=240-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
    expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseAwareNoseV226('"+out.as_posix()+"');"
    return guard.supervise(stage,['D:/matlab/bin/matlab.exe','-batch',expression],min(120,remaining),128,locks,ROOT,'matlab')
if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
