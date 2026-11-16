"""模块职责：绑定非相似工程链的成熟方法、实际源、固定MATLAB预算及独立审核。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v252_validation';DOC=ROOT/'docs/validation/v252';CLOSED=True
OWN=['solvePhaseAwareNonsimilarStation.m','marchPhaseAwareLaminarProfiles.m','phaseAwareNoseMarchInputs.m',
 'mapPhaseAwareMarchToSurface.m','tests/testPhaseNonsimilarV252.m']
LIB=['phaseAwareAirState.m','phaseAwareAirCaloric.m','phaseAwarePayloadSHA.m','phaseAwareMeshIdentity.m',
 'buildPhaseAwareNoseViewerData.m','data/thermo/phase_aware_dry_air_v1.json']
INPUT=[ROOT/'outputs/v104_validation/baseline02/current_case.mat',ROOT/'outputs/v226_validation/matlab02/rocket_faces.csv',
 ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json',ROOT/'outputs/v227_validation/audit02/heating.json']
TEXT=ROOT/'outputs/v52_validation/research/tp3271_text.txt'

def setup():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    return guard

def prepare():
    """模块一：固定同相态沿程方程与成熟解法决定，不重启冻结v61或旧CFD。"""
    out=BASE/'method01';out.mkdir(exist_ok=False)
    text=TEXT.read_text(encoding='utf8');start=text.index('PDF PAGE 5');end=text.index('PDF PAGE 9')
    (out/'tp3271_equations_excerpt.txt').write_text(text[start:end],encoding='utf8')
    bindings={name:sha(ROOT/name) for name in OWN+LIB}
    for name in ['regionalFlow/nonsimilar/advanceNonsimilarLaminarStation.m','regionalFlow/nonsimilar/nonsimilarLaminarResidual.m',
      'solvePhaseAwareLaminarProfile.m','solvePhaseAwareNoseHeating.m','docs/NONSIMILAR_BOUNDARY_LAYER_DRAFT.md',
      'docs/LABELLED_PRESSURE_SURFACE_FLOW_V96.md','scripts/tetraStreamtubes.py']:
        bindings[name]=sha(ROOT/name)
    for name,h in read(ROOT/'outputs/v176_validation/protected_initial.json').items():assert sha(ROOT/name)==h
    save(out/'report.json',dict(passed=True,decision='adopt TP3271 laminar Levy-Lees + backward-xi + mature coupled bvp4c',
      full3DCrossflowAdopted=False,oldV61Frozen=True,sourceURL='https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf',
      sourceTextSHA256=sha(TEXT),equations=[6,7,19,25,27,54,56],
      equationAuthor='Hamilton, Millman and Greendyke, NASA TP-3271 (1992)',
      implementation='Own transformed residual / local installed MATLAB bvp4c; not the report tridiagonal algorithm',
      license='NASA report equations cited; no MathWorks or SciPy source redistributed; existing licensed MATLAB invoked',
      externalState='constant H/S along each supplied line; no cross-stream pressure-momentum equation',
      oldPhysicalN='Rejected for new production adoption: frozen custom Newton and prior integral-balance failure',
      researchOnly=True,newQualifiedAreaM2=0,codeBindings=bindings,inputSHA256={str(p):sha(p) for p in INPUT},
      fixedBranches=['2deg Tw300','1deg Tw300','2deg Tw pulse 15-35deg return300'],noNewCFD=True))

def main(mode):
    if CLOSED:raise SystemExit('v252 phase nonsimilar tasks closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);g=setup()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',DOC/'PLAN.md',TEXT]+INPUT+[ROOT/x for x in OWN+LIB]
    if mode=='prepare':return supervise(g,BASE,mode,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,256,locks,BASE,'offline')
    if mode in ('matlab01','matlab02'):
        out=BASE/mode;out.mkdir(exist_ok=False)
        for rel in OWN+LIB:
            dest=out/'executed_source'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,dest)
        if mode=='matlab02':assert not read(BASE/'process_matlab01/report.json')['passed']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_matlab*/report.json'))
        expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseNonsimilarV252('"+out.as_posix()+"');"
        return supervise(g,BASE,mode,['D:/matlab/bin/matlab.exe','-batch',expression],min(300,600-used),256,locks+[BASE/'method01/report.json'],ROOT,'matlab')
    raise ValueError('MODE_NOT_ADMITTED')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v252 phase nonsimilar tasks closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
