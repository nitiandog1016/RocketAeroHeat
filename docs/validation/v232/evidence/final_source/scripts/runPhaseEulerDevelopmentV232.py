"""模块职责：绑定修复对象，添加只读前缀/拒绝记录，并监督唯一真实外流发展。"""
from pathlib import Path
import difflib,json,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_BUILD
BASE=ROOT/'outputs/v232_validation';DOC=ROOT/'docs/validation/v232'
PREV=ROOT/'outputs/v231_validation';BASELINE=ROOT/'outputs/v230_validation'
CASE=BASE/'rocket01';BUILD=BASE/'build01';HEADER=ROOT/'native/su2PhaseAwareAir/EulerDevelopmentObservationV232.hpp'
CLOSED=True

def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False);CASE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    prev=read(PREV/'final01/report.json');assert prev['passed'] and prev['nativePrimitiveDomainRejectionFixed']
    assert sha(PREV/'recovery_build02/primitive.obj')==prev['primitiveObjectSHA256']
    deps=read(ROOT/'docs/validation/v231/evidence/actual_dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    original=ROOT/'native/su2PhaseAwareAir/rocketEulerRun.cpp';text=original.read_text(encoding='utf8');before=text
    text=text.replace('#include "rocketEulerFailure.hpp"','#include "rocketEulerFailure.hpp"\n#include "EulerDevelopmentObservationV232.hpp"')
    text=text.replace('    observation::failure_phase="native-iterate";',
      '    observation::verify_native_prefix(config[zone]->GetInnerIter());\n    observation::failure_phase="native-iterate";')
    line='    observation::observe(geometry[zone][inst][0],solver[zone][inst][0][FLOW_SOL],config[zone]);'
    assert text.count(line)==1;text=text.replace(line,line+'\n    observation::record_rejections_and_check(config[zone]);')
    text=text.replace('    observation::frame(g,f,0,false);','    observation::frame(g,f,0,false);\n    observation::load_baseline();')
    source=out/'rocketEulerObserved.cpp';source.write_text(text,encoding='utf8');shutil.copyfile(HEADER,out/HEADER.name)
    (out/'observation.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),
      fromfile='v231/rocketEulerRun.cpp',tofile='v232/rocketEulerObserved.cpp')),encoding='utf8')
    for n in ('case.cfg','mesh.su2'):shutil.copyfile(BASELINE/'rocket01'/n,CASE/n)
    for src,dst in [('physical_history.csv','baseline_physical.csv'),('history.csv','baseline_native.csv')]:
        shutil.copyfile(BASELINE/'rocket01'/src,CASE/dst)
    cmd=read(PREV/'recovery_build03/command.json')
    cmd=[str(source) if x==str(original) else x.replace(str(PREV/'recovery_build03/rocket.exe'),str(BUILD/'rocket.exe'))
        .replace(str(PREV/'recovery_build03/rocket.d'),str(BUILD/'rocket.d')) for x in cmd]
    save(BUILD/'command.json',cmd)
    filtered=[];i=0
    while i<len(cmd):
        if cmd[i] in ('-MF','-o'):i+=2;continue
        if cmd[i] in ('-MD','-Wl,--start-group','-Wl,--end-group') or cmd[i].endswith(('.a','.obj')):i+=1;continue
        filtered.append(cmd[i]);i+=1
    dep=out/'driver.d';dep_cmd=filtered+['-M','-MT','dep','-MF',str(dep)]
    subprocess.check_call(dep_cmd,cwd=str(OLD_BUILD),timeout=20)
    bindings={}
    for n in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
        if p==out/HEADER.name:assert h==sha(HEADER)
        elif p!=source:assert deps[str(p)]==h,str(p)
        bindings[str(p)]=h
    save(out/'dependencies.json',bindings);save(out/'dependency_command.json',dep_cmd)
    solver=ROOT/'outputs/v229_validation/source/SU2_CFD/src/solvers/CEulerSolver.cpp'
    expected=read(ROOT/'outputs/v229_validation/process_build/report.json')['inputSHA256'][str(solver)]
    assert sha(solver)==expected
    src=solver.read_text(encoding='utf8');a=src.index('void CEulerSolver::CommonPreprocessing');b=src.index('void CEulerSolver::Preprocessing',a)
    (out/'actual_count_source.cpp.txt').write_text(src[a:b],encoding='utf8')
    save(out/'report.json',dict(passed=True,driverSourceSHA256=sha(source),originalDriverSHA256=sha(original),
      primitiveObjectSHA256=sha(PREV/'recovery_build02/primitive.obj'),sourceCounterSHA256=expected,
      dependencies=len(bindings),physicsObjectsChanged=False,sourceAdaptation='observation only',flowExecuted=False))
    save(CASE/'contract.json',dict(role='unique-repaired-primitive-Euler-development',updatesLimit=800,wallSecondsLimit=600,
      sameConfigAsV230=True,meshSHA256=sha(CASE/'mesh.su2'),configSHA256=sha(CASE/'case.cfg'),
      expectedStates=[0,1]+list(range(50,801,50)),phaseDomainUnchanged=True,newInitialState='same compatible uniform state',
      prefixFirst30MustEqual=True,rejectionCountRole='last successful native preprocessing only; not all substage cumulative',
      physicalQualification=False,absoluteLogRMSLimit=-8,lastPairPressureAreaChangeLimit=.01,massEnergyBalanceLimit=1e-8,
      spatialFarfieldQualified=False,automaticRetries=0))

def locks():
    items=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',HEADER,
      PREV/'final01/report.json',PREV/'recovery_build03/command.json',PREV/'recovery_build02/primitive.obj']
    items += [Path(n) for n in read(ROOT/'docs/validation/v231/evidence/actual_dependencies.json')]
    for n in read(PREV/'recovery_build03/command.json'):
        if n.endswith('.a'):items.append(Path(n))
    return list(dict.fromkeys(items))

def main(stage):
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    inputs=locks();kind='offline';seconds=120;cwd=BASE
    if stage=='prepare':cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare']
    elif stage=='build':
        assert read(BASE/'prepare01/report.json')['passed'];cmd=read(BUILD/'command.json');kind='build';seconds=60
        inputs += [Path(n) for n in read(BASE/'prepare01/dependencies.json')]+[BUILD/'command.json']
    elif stage=='run':
        assert read(BASE/'process_build/report.json')['passed'] and not (BASE/'process_run').exists()
        inputs += [p for p in CASE.iterdir() if p.is_file()]+[BUILD/'rocket.exe']
        cmd=[str(BUILD/'rocket.exe'),'case.cfg','.'];cwd=CASE;kind='cfd';seconds=600
    else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise(stage,cmd,seconds,512,list(dict.fromkeys(inputs)),cwd,kind,False)

if __name__=='__main__':
    if CLOSED:raise SystemExit('v232 candidate closed; no retry or continuation')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
