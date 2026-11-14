"""模块职责：绑定单连通生产循环准确源码，生成只读观测适配与固定构建输入。"""
from pathlib import Path
import json,re,shutil,subprocess,sys,difflib,posixpath
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v235_validation';DOC=ROOT/'docs/validation/v235';GAS=ROOT/'outputs/v220_validation'
OLD=ROOT/'tools/runtime/su2-thermal-build/source';BUILDOLD=OLD.parent/'build-v117'
CLOSED=False
NATIVE=ROOT/'native/su2PhaseAwareAir';STAGE=BASE/'stage01';BUILD=BASE/'build01'
SPECS=['SU2_CFD/src/drivers/CSinglezoneDriver.cpp','SU2_CFD/src/iteration/CFluidIteration.cpp',
 'SU2_CFD/src/integration/CIntegration.cpp','SU2_CFD/src/integration/CMultiGridIntegration.cpp',
 'SU2_CFD/src/integration/CSingleGridIntegration.cpp','SU2_CFD/src/variables/CNSVariable.cpp',
 'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
 'SU2_CFD/src/solvers/CTurbSolver.cpp','SU2_CFD/src/solvers/CTurbSASolver.cpp',
 'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl']


def sources():
    target=BASE/'source01';target.mkdir(exist_ok=False);known=read(GAS/'build02/dependencies.json')
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');bound=[]
    for rel in SPECS:
        expected=known.get(str(GAS/'source'/rel),tree.get(str(OLD/rel)));assert expected,rel
        candidates=[GAS/'source'/rel,ROOT/'outputs/v222_validation/source01'/rel,
          ROOT/'outputs/v167_validation/build01/executed_sources'/rel,
          ROOT/'outputs/v219_validation/source01/source'/rel,OLD/rel]
        p=next((p for p in candidates if p.is_file() and sha(p)==expected),None)
        assert p,('EXACT_EXECUTED_SOURCE_MISSING',rel)
        dest=target/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        bound.append(dict(relative=rel,source=str(p),sha256=expected,snapshot=str(dest)))
    save(target/'bindings.json',bound)
    save(target/'report.json',dict(passed=True,sourceFiles=len(bound),flowUpdates=0,nativeExecution=False))


def stage():
    """模块二：仅添加阶段/拒绝回调，直接绑定原v220头和成熟归档。"""
    STAGE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False);known=read(GAS/'build02/dependencies.json')
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');commands={};patch=[];changes=[]
    hook=(NATIVE/'SingleDomainHooksV235.hpp').as_posix()
    def insert(text,anchor,addition):
        assert text.count(anchor)==1,(anchor,text.count(anchor));return text.replace(anchor,anchor+addition)
    for rel in ('SU2_CFD/src/integration/CIntegration.cpp','SU2_CFD/src/variables/CNSVariable.cpp'):
        p=BASE/'source01'/rel;before=p.read_text(encoding='utf8');t=before
        def includes(m):
            path=posixpath.normpath(posixpath.join(posixpath.dirname(rel),m.group(1)));q=GAS/'source'/path
            assert q.is_file() and sha(q)==known.get(str(q),tree.get(str(OLD/path))),str(q)
            return '#include "'+q.as_posix()+'"'
        t=re.sub(r'#include "(\.\./[^\"]+)"',includes,t)
        t='#include "'+hook+'"\n'+t
        if rel.endswith('CIntegration.cpp'):
            call=lambda label:'\n  if (rah_single_observer::stage) rah_single_observer::stage("'+label+'",geometry,solver_container,config,RunTime_EqSystem);\n'
            anchor='  solver_container[MainSolver]->Source_Residual(geometry, solver_container, numerics, config, iMesh);'
            t=insert(t,anchor,call('after_sources'))
            anchor='  /*--- Modification of the system on the whole domain, like for a strong BC. ---*/'
            assert t.count(anchor)==1;t=t.replace(anchor,call('after_weak')+anchor)
            anchor='  solver_container[MainSolver]->Impose_Fixed_Values(geometry, config);'
            t=insert(t,anchor,call('before_strong'))
            t=insert(t,'  //AD::ResumePreaccumulation(pausePreacc);',call('after_strong'))
            k=t.index('void CIntegration::Time_Integration');head,tail=t[:k],t[k:]
            tail=insert(tail,'  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);',call('before_time'))
            anchor='  }\n\n}\n\nvoid CIntegration::SetDualTime_Geometry';assert tail.count(anchor)==1
            tail=tail.replace(anchor,'  }\n'+call('after_time')+'\n}\n\nvoid CIntegration::SetDualTime_Geometry');t=head+tail
        else:
            anchor='      throw;\n    }\n    if (air_state_rejected)';assert t.count(anchor)==1
            t=t.replace(anchor,'      if (rah_single_observer::rejection) rah_single_observer::rejection(this,iPoint,trial.data(),restored.data(),false,error.what());\n'+anchor)
            t=insert(t,'    RightVol = false;','\n    if (rah_single_observer::rejection) rah_single_observer::rejection(this,iPoint,trial.data(),restored.data(),true,rejected_reason.c_str());\n')
        target=STAGE/Path(rel).name;target.write_text(t,encoding='utf8')
        patch+=list(difflib.unified_diff(before.splitlines(True),t.splitlines(True),fromfile=rel,tofile='observed/'+rel))
        changes.append(dict(source=str(p),originalSHA256=sha(p),patched=str(target),patchedSHA256=sha(target)))
        cmd=read(GAS/'prepare03/commands.json')['SU2_CFD/src/variables/CNSVariable.cpp']
        cmd=list(cmd);cmd[cmd.index('-c')+1]=str(target);obj=BUILD/(Path(rel).stem+'.obj');cmd[cmd.index('-o')+1]=str(obj)
        commands[rel]=cmd+['-MD','-MF',str(obj)+'.d']
    helper=read(ROOT/'outputs/v234_validation/build03/command.json')
    helper=[str(NATIVE/'singleDomainRunV235.cpp') if x.endswith('multiZoneCompatibilityV234.cpp') else
      x.replace(str(ROOT/'outputs/v234_validation/build03'),str(BUILD)).replace(str(GAS/'build02'),str(BUILD)) for x in helper]
    assert str(NATIVE/'singleDomainRunV235.cpp') in helper;commands['helper']=helper
    (STAGE/'observation.patch').write_text(''.join(patch),encoding='utf8');save(STAGE/'changes.json',changes);save(STAGE/'commands.json',commands)
    bindings={}
    for j,cmd in enumerate(commands.values()):
        a=[];i=0
        while i<len(cmd):
            x=cmd[i]
            if x in ('-o','-MF'):i+=2;continue
            if x in ('-c','-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            a.append(x);i+=1
        dep=STAGE/('unit%d.d'%j);subprocess.check_call(a+['-M','-MT','dep','-MF',str(dep)],cwd=str(BUILDOLD),timeout=20)
        for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(name);p=(p if p.is_absolute() else BUILDOLD/p).resolve();h=sha(p)
            if p.parent not in (STAGE,NATIVE):assert h==known.get(str(p),tree.get(str(p))),str(p)
            bindings[str(p)]=h
    save(STAGE/'dependencies.json',bindings)
    save(STAGE/'report.json',dict(passed=True,compiledObjectsPlanned=2,completeNativeLoopPreserved=True,
      callbacksOnly=True,sourceFiles=len(changes),dependencies=len(bindings),nativeExecution=False,flowUpdates=0))


def main(stage):
    if CLOSED:raise SystemExit('v235 preparation closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      GAS/'build02/dependencies.json',ROOT/'outputs/v167_validation/build01/source_tree.json']
    if stage=='sources':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_sources'],120,128,locks,BASE)
    if stage=='stage':
        locks += [NATIVE/'singleDomainRunV235.cpp',NATIVE/'SingleDomainHooksV235.hpp']
        locks += [Path(r['snapshot']) for r in read(BASE/'source01/bindings.json')]
        return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_stage'],120,128,locks,BASE)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v235 preparation closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_sources':sources()
    elif sys.argv[1]=='task_stage':stage()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
