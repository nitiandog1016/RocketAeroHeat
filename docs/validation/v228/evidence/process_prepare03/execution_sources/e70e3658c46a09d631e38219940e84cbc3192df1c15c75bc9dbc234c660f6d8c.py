"""模块职责：只执行登记的壁列回归/离线准备，明确禁止构建和CFD。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,sha
BASE=ROOT/'outputs/v228_validation';DOC=ROOT/'docs/validation/v228';CLOSED=False
SOURCES=['solvePhaseAwareLaminarProfile.m','buildPhaseAwareNoseHeatingViewerData.m','buildPhaseAwareNoseViewerData.m',
 'buildAerothermalViewerData.m','phaseAwareAirState.m','phaseAwareAirCaloric.m','phaseAwarePayloadSHA.m',
 'phaseAwareMeshIdentity.m','tests/testPhaseNoseRecoveryV228.m','data/thermo/phase_aware_dry_air_v1.json']
def main(stage):
    if CLOSED:raise SystemExit('v228 execution closed')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda f:0 if f.parent==BASE and f.name.startswith('audit') else size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage in ('matlab01','matlab02'):
        out=BASE/stage;out.mkdir(exist_ok=False)
        if stage=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
        for rel in SOURCES:
            p=out/'executed_source'/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,p)
        locks += [ROOT/n for n in SOURCES]+[ROOT/'outputs/v104_validation/baseline02/current_case.mat',
          ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json',ROOT/'outputs/v226_validation/matlab02/rocket_faces.csv']
        locks += [ROOT/'docs/validation/v227/evidence/recovered'/n for n in ('heating.json','report.json','rocket_workbench.json')]
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
        expr="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseNoseRecoveryV228('"+out.as_posix()+"');"
        return guard.supervise(stage,['D:/matlab/bin/matlab.exe','-batch',expr],min(90 if stage=='matlab01' else 60,150-used),128,locks,ROOT,'matlab')
    assert stage in ('prepare01','prepare02','prepare03','audit01')
    script=ROOT/'scripts/preparePhaseEulerV228.py'
    locks += [script,ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',
      ROOT/'outputs/v194_validation/rocket01/mesh.su2',ROOT/'outputs/v194_validation/rocket01/surface.json',
      ROOT/'outputs/v194_validation/native_audit01/report.json',ROOT/'outputs/v194_validation/rocket01/manifest.json',
      ROOT/'outputs/v220_validation/build02/dependencies.json',ROOT/'outputs/v167_validation/build01/source_tree.json',
      ROOT/'outputs/v222_validation/rocket01/case.cfg',ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json',
      ROOT/'outputs/v226_validation/matlab02/rocket_faces.csv',ROOT/'docs/validation/v212/evidence/design01/regions.json',
      ROOT/'outputs/v222_validation/process_engine/report.json',ROOT/'outputs/v220_validation/build02/SU2_CFD.exe']
    from preparePhaseEulerV228 import bound_sources
    sources,_=bound_sources();locks += [Path(v['path']) for v in sources.values()]
    if stage=='audit01':locks += [p for p in (BASE/'prepare03').iterdir() if p.is_file()]
    return guard.supervise(stage,[sys.executable,str(script),stage],120,128,locks,BASE,'offline')
if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
