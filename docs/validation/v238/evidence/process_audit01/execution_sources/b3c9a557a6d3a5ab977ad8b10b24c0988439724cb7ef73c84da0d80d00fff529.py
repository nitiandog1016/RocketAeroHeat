"""模块职责：准确绑定成熟HLLC对象与必要控制，预算内隔离构建；不自动准入真实CFD。"""
from pathlib import Path
import json,os,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v238_validation';DOC=ROOT/'docs/validation/v238';STAGE=BASE/'prepare01';BUILD=BASE/'build01';CASE=BASE/'interface01'
CPP=ROOT/'native/su2PhaseAwareAir/hllcInterfaceV238.cpp';OLD=ROOT/'outputs/v235_validation';GAS=ROOT/'outputs/v220_validation'
OLDBUILD=ROOT/'tools/runtime/su2-thermal-build/build-v117';CLOSED=False


def configure():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes
    guard.size_bytes=lambda p:original(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0


def prepare():
    """模块一：编译器实际依赖、成熟归档成员和准确通量/工厂源一次绑定。"""
    for p in (STAGE,BUILD,CASE):p.mkdir(exist_ok=False)
    cmd=read(ROOT/'outputs/v237_validation/stage02/commands.json')['control']
    cmd=[str(CPP) if x.endswith('phaseStepControlV237.cpp') else str(BUILD/'interface.d') if x.endswith('control.d') else
         str(BUILD/'interface.exe') if x.endswith('control.exe') else x.replace(str(ROOT/'outputs/v237_validation/build01'),str(OLD/'build01')) if x.endswith('.a') else x for x in cmd]
    save(STAGE/'command.json',cmd);filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(STAGE/'inputs.d')],cwd=str(OLDBUILD),timeout=20)
    known=read(ROOT/'outputs/v237_validation/stage02/inputs.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');deps={}
    for raw in (STAGE/'inputs.d').read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(raw);p=(p if p.is_absolute() else OLDBUILD/p).resolve();h=sha(p)
        if p!=CPP:
            expected=known.get(str(p),tree.get(str(p)))
            if expected is None:expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/p.relative_to(GAS/'source'))]
            assert h==expected,str(p)
        deps[str(p)]=h
    for rows in read(OLD/'build01/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256'];deps[row['path']]=row['sha256']
    for row in read(ROOT/'outputs/v237_validation/outcome01/report.json')['sourceBindings']:
        assert sha(row['path'])==row['sha256'];deps[row['path']]=row['sha256']
    shutil.copyfile(OLD/'observe02/case.cfg',CASE/'case.cfg')
    deps[str(CASE/'case.cfg')]=sha(CASE/'case.cfg')
    save(STAGE/'inputs.json',deps);save(STAGE/'report.json',dict(passed=True,boundInputs=len(deps),numericalObjectsRecompiled=0,realCFDAdmitted=False))


def main(stage):
    if CLOSED:raise SystemExit('v238 execution closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',CPP]
    if stage=='prepare':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE)
    assert read(STAGE/'report.json')['passed']
    for p,h in read(STAGE/'inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
    locks += [STAGE/'inputs.json',STAGE/'command.json']
    if stage=='build':return guard.supervise(stage,read(STAGE/'command.json'),180,1024,locks,OLDBUILD,'build')
    if stage=='interface':
        assert read(BASE/'process_build/report.json')['passed'];locks.append(BUILD/'interface.exe')
        return guard.supervise(stage,[str(BUILD/'interface.exe'),'case.cfg','.'],60,1024,locks,CASE,'query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
