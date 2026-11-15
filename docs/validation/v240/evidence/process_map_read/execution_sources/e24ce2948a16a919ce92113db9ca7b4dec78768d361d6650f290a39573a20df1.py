"""模块职责：恢复全局号到内部号映射，以保存末场实际调用新统计函数，不推进流动。"""
from pathlib import Path
import json,shutil,subprocess,sys
import runForebodyEulerV240 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;FIX=BASE/'fix01';CONTROL=BASE/'readback01'
CLOSED=False


def prepare():
    CONTROL.mkdir(exist_ok=False)
    previous=read(BASE/'process_run/report.json')['inputSHA256']
    header=run.NATIVE/'forebodyEulerV240.hpp';assert sha(FIX/'forebodyEulerV240.executed.hpp')==previous[str(header)]
    source=run.PREP/'forebodyEuler.cpp';text=source.read_text(encoding='utf8')
    text=text.replace('    region240::initialize(g,f,cfg);','    region240::initialize(g,f,cfg);\n    if(region240::surface_readback){region240::readback(g,f);return;}')
    anchor='if(argc!=3)throw std::runtime_error("CONFIG_AND_OUTPUT_REQUIRED");'
    assert text.count(anchor)==1;text=text.replace(anchor,'if(argc!=4||std::string(argv[3])!="--surface-readback")throw std::runtime_error("READBACK_ONLY_NO_FLOW"); region240::surface_readback=true;')
    target=FIX/'forebodyReadback.cpp';target.write_text(text,encoding='utf8')
    command=read(run.BUILD/'link_command.json')
    replacements={str(source):str(target),str(run.BUILD/'forebody.exe'):str(FIX/'readback.exe'),str(run.BUILD/'driver.d'):str(FIX/'readback.d')}
    command=[replacements.get(s,s) for s in command];save(FIX/'link_command.json',command)
    for name in ('case.cfg','mesh.su2','original_faces.csv'):shutil.copyfile(run.CASE/name,CONTROL/name)
    shutil.copyfile(run.CASE/'state_000200.csv',CONTROL/'saved_final.csv')
    inputs=read(run.PREP/'dependencies.json')
    for p,h in inputs.items():
        if Path(p)==header:continue
        assert sha(p)==h,p
    inputs[str(header)]=sha(header);inputs[str(target)]=sha(target)
    inputs[str(run.BUILD/'boundary.obj')]=sha(run.BUILD/'boundary.obj')
    save(FIX/'dependencies.json',inputs)
    save(FIX/'report.json',dict(passed=True,originalExecutedHeaderSHA256=sha(FIX/'forebodyEulerV240.executed.hpp'),
      currentCorrectedHeaderSHA256=sha(header),originalCFDReplayed=False,flowUpdatesAdmitted=0,numericalObjectsRecompiled=0))


def main(stage):
    if CLOSED:raise SystemExit('v240 mapping recovery closed; no flow replay')
    run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC
    size=run.guard.size_bytes;run.guard.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runForebodyEulerV240.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'REGIONAL_PLAN.md',
      run.NATIVE/'forebodyEulerV240.hpp',FIX/'forebodyEulerV240.executed.hpp',BASE/'process_run/report.json']
    if stage=='prepare':
        locks += [run.PREP/'forebodyEuler.cpp',run.PREP/'dependencies.json',run.BUILD/'link_command.json']+[run.CASE/p for p in ('case.cfg','mesh.su2','original_faces.csv','state_000200.csv')]
        return run.guard.supervise('map_prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],30,256,locks,BASE)
    for p,h in read(FIX/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
    if stage=='build':
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        locks += [FIX/'link_command.json',FIX/'dependencies.json']
        return run.guard.supervise('map_build',read(FIX/'link_command.json'),120-used,256,locks,run.OLD_BUILD,'build')
    if stage=='read':
        assert read(BASE/'process_map_build/report.json')['passed']
        locks += [FIX/'readback.exe']+[p for p in CONTROL.iterdir() if p.is_file()]
        return run.guard.supervise('map_read',[str(FIX/'readback.exe'),'case.cfg','.','--surface-readback'],15,256,locks,CONTROL,'query')
    raise ValueError('NO_ADMISSION')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v240 mapping recovery closed; no flow replay')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
