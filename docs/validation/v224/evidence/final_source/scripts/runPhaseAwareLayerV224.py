"""模块职责：监督新工程层的有限MATLAB验证，绑定源码/数据并保留所有失败。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v224_validation';DOC=ROOT/'docs/validation/v224'
SOURCES=['phaseAwareAirState.m','phaseAwareAirReferenceState.m','phaseAwareMeshIdentity.m',
    'calcPhaseAwareLaminarHeating.m','calcSurfaceHeating.m','buildPhaseAwareLaminarViewerData.m',
    'buildAerothermalViewerData.m','runAeroHeatSingleCase.m','surfaceWallTemperature.m',
    'calcBoundaryLayerState.m','airReferenceTemperature.m','airRecoveryTemperature.m',
    'airThermalProperties.m','selectAirThermalModel.m','tests/testPhaseAwareLaminarHeating.m',
    'tests/testBoundaryLayerState.m','data/thermo/phase_aware_dry_air_v1.json']
CLOSED=True

def main(name):
    if CLOSED:raise SystemExit('v224 validation entry closed')
    assert name in ('matlab01','matlab02')
    BASE.mkdir(exist_ok=True);out=BASE/name;out.mkdir(exist_ok=False)
    for rel in SOURCES:
        p=ROOT/rel;q=out/'executed_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    locks=[ROOT/n for n in SOURCES]+[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    seconds=180-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
    assert seconds>0
    if name=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
    expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseAwareLaminarHeating('"+out.as_posix()+"');"
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda unused:size(BASE)
    return guard.supervise(name,['D:/matlab/bin/matlab.exe','-batch',expression],min(90,seconds),128,locks,ROOT,'matlab')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
