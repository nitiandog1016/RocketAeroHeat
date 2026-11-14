"""模块职责：准确绑定成熟HLLC对象与必要控制，预算内隔离构建；不自动准入真实CFD。"""
from pathlib import Path
import difflib,json,os,posixpath,re,shutil,subprocess,sys,time
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
    if stage=='repair_prepare':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_repair'],60,1024,locks,BASE)
    assert read(STAGE/'report.json')['passed']
    for p,h in read(STAGE/'inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
    locks += [STAGE/'inputs.json',STAGE/'command.json']
    if stage=='build':return guard.supervise(stage,read(STAGE/'command.json'),180,1024,locks,OLDBUILD,'build')
    if stage=='interface':
        assert read(BASE/'process_build/report.json')['passed'];locks.append(BUILD/'interface.exe')
        return guard.supervise(stage,[str(BUILD/'interface.exe'),'case.cfg','.'],60,1024,locks,CASE,'query')
    if stage in ('build2','interface2'):
        for p,h in read(BASE/'repair01/inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
        locks.append(BASE/'repair01/commands.json')
        if stage=='build2':
            remaining=180-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
            return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build2'],remaining,1024,locks,OLDBUILD,'build')
        assert read(BASE/'process_build2/report.json')['passed'];exe=BASE/'build02/interface.exe';locks.append(exe)
        return guard.supervise(stage,[str(exe),'case.cfg','.'],50,1024,locks,BASE/'interface02','query')
    raise ValueError('STAGE_NOT_ADMITTED')


def repair_prepare():
    """模块二：只修复独立一致态控制已证实的四处守恒压力导数符号。"""
    stage=BASE/'repair01';build=BASE/'build02';case=BASE/'interface02'
    for p in (stage,build,case):p.mkdir(exist_ok=False)
    source=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/numerics/flow/convection/hllc.cpp'
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');assert sha(source)==tree[str(source)]
    before=source.read_text(encoding='utf8');text=before
    for side in ('i','j'):
        old='dPI_dU[0] = Chi_%s - 0.5 * Kappa_%s * sq_vel_%s;'%(side,side,side)
        assert text.count(old)==2;text=text.replace(old,old.replace(' - 0.5',' + 0.5'))
    rel='SU2_CFD/src/numerics/flow/convection/hllc.cpp'
    text=re.sub(r'#include "(\.\./[^\"]+)"',lambda m:'#include "'+(GAS/'source'/posixpath.normpath(posixpath.join(posixpath.dirname(rel),m.group(1)))).as_posix()+'"',text)
    target=stage/'hllc.cpp';target.write_text(text,encoding='utf8')
    (stage/'hllc_pressure_derivative.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=rel,tofile='pressure-derivative/'+rel)),encoding='utf8')
    cmd=read(GAS/'prepare03/commands.json')['SU2_CFD/src/solvers/CEulerSolver.cpp'];cmd=list(cmd)
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(build/'hllc.obj')
    helper=[x.replace(str(BUILD),str(build)) if x.endswith(('.exe','.d')) else x.replace(str(OLD/'build01'),str(build)) if x.endswith('.a') else x for x in read(STAGE/'command.json')]
    save(stage/'commands.json',dict(hllc=cmd,helper=helper))
    filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x=='-o':i+=2;continue
        if x=='-c':i+=1;continue
        filtered.append(x);i+=1
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(stage/'inputs.d')],cwd=str(OLDBUILD),timeout=20)
    deps={};known=read(STAGE/'inputs.json')
    for raw in (stage/'inputs.d').read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(raw);p=(p if p.is_absolute() else OLDBUILD/p).resolve();h=sha(p)
        if p!=target:
            expected=known.get(str(p),tree.get(str(p)))
            if expected is None:expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/p.relative_to(GAS/'source'))]
            assert h==expected,str(p)
        deps[str(p)]=h
    shutil.copyfile(CASE/'case.cfg',case/'case.cfg');deps[str(case/'case.cfg')]=sha(case/'case.cfg')
    save(stage/'inputs.json',deps);save(stage/'report.json',dict(passed=True,onlyFourDerivativeSignsChanged=True,originalSourceSHA256=sha(source),patchedSourceSHA256=sha(target),newNumericalObjects=1))


def build2():
    """模块三：隔离编译一个有限修正对象，成熟其余对象全部复用。"""
    stage=BASE/'repair01';build=BASE/'build02';commands=read(stage/'commands.json');t=time.monotonic()
    subprocess.check_call(commands['hllc'],cwd=str(OLDBUILD));save(build/'compile.json',dict(hllcSeconds=time.monotonic()-t))
    ar=str(Path(commands['hllc'][0]).parent/'ar.exe');members={};count=0
    for name,rows in read(OLD/'build01/members.json').items():
        result=[]
        for r in rows:
            p=Path(r['path']);assert sha(p)==r['sha256']
            if 'convection_hllc.cpp.obj' in p.name:p=build/'hllc.obj';count+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(build/name)]+[r['path'] for r in result]);members[name]=result
    assert count==1;save(build/'members.json',members);subprocess.check_call(commands['helper'],cwd=str(OLDBUILD))
    save(build/'report.json',dict(passed=True,nativeObjectsChanged=1,helperSHA256=sha(build/'interface.exe')))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':prepare()
    elif sys.argv[1]=='task_repair':repair_prepare()
    elif sys.argv[1]=='task_build2':build2()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
