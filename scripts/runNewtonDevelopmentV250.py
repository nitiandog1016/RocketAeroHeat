"""模块职责：绑定一个原生NK观察对象，监督同初态的唯一有限发展及真实终态容量。"""
from pathlib import Path
import difflib
import os
import re
import shutil
import subprocess
import sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise

BASE=ROOT/'outputs/v250_validation';DOC=ROOT/'docs/validation/v250'
CASE=BASE/'rocket01';STAGE=BASE/'stage01';BUILD=BASE/'build01';CLOSED=True
PREV=ROOT/'outputs/v249_validation';GAS=ROOT/'outputs/v220_validation'
WD=ROOT/'tools/runtime/su2-thermal-build/build-v117';NATIVE=ROOT/'native/su2PhaseAwareAir'
DRIVER=NATIVE/'newtonDevelopmentV250.cpp';HOOK=NATIVE/'NewtonHooksV250.hpp'
ORIGINAL=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/integration/CNewtonIntegration.cpp'
REFERENCE=ROOT/'outputs/v248_validation/rocket01/state_000200.bin'
FACES=ROOT/'outputs/v236_validation/rocket01/original_faces.txt'
MEMBERS=ROOT/'outputs/v245_validation/build01/members.json'


def setup_guard():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes
    guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    return guard


def prepare():
    """模块一：原算式逐字保留，只插可移除回调；实际输入/最坏保存数量先算容量。"""
    assert read(PREV/'audit03/report.json')['nativeControlCompleted']
    method=read(PREV/'method01/report.json')
    expected=next(r['sha256'] for r in method['nativeSourceExcerpts'] if Path(r['path'])==ORIGINAL)
    assert sha(ORIGINAL)==expected
    for p in (CASE,STAGE,BUILD):p.mkdir(exist_ok=False)
    old=ORIGINAL.read_text(encoding='utf8')
    includes={
      '../../include/integration/CNewtonIntegration.hpp':GAS/'source/SU2_CFD/include/integration/CNewtonIntegration.hpp',
      '../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp':GAS/'source/Common/include/linear_algebra/CMatrixVectorProduct.hpp'}
    text=old
    for rel,p in includes.items():assert p.is_file();text=text.replace('"'+rel+'"','"'+p.as_posix()+'"')
    normalized=text;insertions=[]
    def insert(anchor,addition,before=False,count=1):
        nonlocal text
        assert text.count(anchor)==count,anchor
        text=text.replace(anchor,addition+anchor if before else anchor+addition,1)
        insertions.append(addition)
    def hook(label,dr='nullptr',mag='0.0'):
        return '\n  if(rah_newton_observer::event) rah_newton_observer::event("'+label+'",geometry,solvers,config,'+dr+','+mag+',finDiffStep);\n'
    insert('using Scalar = CNewtonIntegration::Scalar;','\n#include "'+HOOK.as_posix()+'"\n')
    insert('  /*--- Save the default integration scheme, and force to explicit if required. ---*/',
      '\n  if(rah_newton_observer::event) rah_newton_observer::event(type==ResEvalType::DEFAULT?"default_begin":"explicit_begin",geometry,solvers,config,nullptr,0.0,finDiffStep);\n',True)
    insert('  Space_Integration(geometry, solvers, numerics[FLOW_SOL], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS);',
      '\n  if(rah_newton_observer::event) rah_newton_observer::event(type==ResEvalType::DEFAULT?"default_end":"explicit_end",geometry,solvers,config,nullptr,0.0,finDiffStep);\n')
    insert('  solvers[FLOW_SOL]->CompleteImplicitIteration(geometry, solvers, config);',hook('update_begin','&linSysSol','1.0'),True)
    insert('  solvers[FLOW_SOL]->CompleteImplicitIteration(geometry, solvers, config);',hook('update_end'))
    insert('  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);',hook('accepted_preprocessing_begin'),True)
    insert('  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);',hook('accepted_preprocessing_end'))
    insert('  }\n  END_SU2_OMP_PARALLEL',hook('flow_complete'),True)
    insert('  PerturbSolution(u, factor);',hook('jv_begin','&u','factor'),True)
    insert('  CSysMatrixComms::Complete(v, geometry, config);',hook('jv_end'),count=2)
    recovered=text
    for addition in reversed(insertions):assert recovered.count(addition)==1;recovered=recovered.replace(addition,'',1)
    assert recovered==normalized
    patched=STAGE/'CNewtonIntegration_observed.cpp';patched.write_text(text,encoding='utf8')
    shutil.copyfile(ORIGINAL,STAGE/'CNewtonIntegration_original.cpp')
    (STAGE/'observer.patch').write_text(''.join(difflib.unified_diff(normalized.splitlines(True),text.splitlines(True),fromfile='native',tofile='observed')),encoding='utf8')
    cmd=list(read(ROOT/'outputs/v235_validation/stage01/commands.json')['SU2_CFD/src/integration/CIntegration.cpp'])
    cmd[cmd.index('-c')+1]=str(patched);cmd[cmd.index('-o')+1]=str(BUILD/'CNewtonIntegration.obj');cmd[cmd.index('-MF')+1]=str(BUILD/'newton.d')
    helper=[str(DRIVER) if x==str(NATIVE/'newtonControlV249.cpp') else x.replace(str(PREV/'build01'),str(BUILD)).replace(str(ROOT/'outputs/v245_validation/build01'),str(BUILD)) for x in read(PREV/'stage01/command.json')]
    save(STAGE/'commands.json',dict(object=cmd,helper=helper))
    known=read(PREV/'stage01/dependencies.json');known.update(read(ROOT/'outputs/v245_validation/stage01/dependencies.json'));deps={}
    for j,command in enumerate((cmd,helper)):
        args=[];i=0
        while i<len(command):
            x=command[i]
            if x in ('-o','-MF'):i+=2;continue
            if x in ('-c','-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            args.append(x);i+=1
        dep=STAGE/('deps%d.d'%j);subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=str(WD),timeout=25)
        for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
            if p not in (DRIVER,HOOK,patched):assert h==known[str(p)],str(p)
            deps[str(p)]=h
    save(STAGE/'dependencies.json',deps)
    for p in (DRIVER,HOOK):shutil.copyfile(p,STAGE/p.name)
    cfg=(PREV/'native01/case.cfg').read_text(encoding='utf8')
    for key,value in (('ITER','200'),('OUTPUT_WRT_FREQ','200')):
        cfg,n=re.subn('^'+key+r'=.*$',key+'= '+value,cfg,count=1,flags=re.M);assert n==1
    (CASE/'case.cfg').write_text(cfg,encoding='utf8')
    for name in ('mesh.su2','input_restart.csv'):os.link(str(PREV/'native01'/name),str(CASE/name))
    shutil.copyfile(FACES,CASE/'original_faces.txt')
    schema=read(PREV/'native01/schema.json');schema.update(maximumUpdates=200,interval=50,
      residualTrend='DEFAULT/baseline after_strong flow; SA after_time, separately staged',
      newtonObservation='read-only callbacks; native physical recovery preserved; Jv recovery rejects candidate',
      normalTerminalSpatial='v248 zero-time flow then SA; six states and SA mut bit identity required')
    save(CASE/'schema.json',schema)
    frame=24+335726*34*8
    components=dict(mesh=(CASE/'mesh.su2').stat().st_size,inputRestart=(CASE/'input_restart.csv').stat().st_size,
      sevenFullFrames=7*frame,finalNativeCSV=(PREV/'native01/restart_flow.csv').stat().st_size,
      terminalTwoBinary=48+335726*15*8,buildAndSourceAllowance=128*1024**2,observationsAndAuditUIAllowance=128*1024**2)
    assert sum(components.values())<1536*1024**2
    save(STAGE/'report.json',dict(passed=True,originalNewtonSHA256=expected,numericalTextUnchangedAfterRemovingHooks=True,
      callbacks=len(insertions)-1,observedObjectsToCompile=1,matureObjectsToRecompile=0,initialSourceSHA256=sha(REFERENCE),
      maximumUpdates=200,capacityComponentsBytes=components,projectedLogicalBytes=sum(components.values()),limitBytes=1536*1024**2))


def build():
    """模块二：仅替换一个只读观察对象，其他成熟对象逐SHA复用为薄归档。"""
    commands=read(STAGE/'commands.json');subprocess.check_call(commands['object'],cwd=str(WD))
    ar=str(Path(commands['helper'][0]).parent/'ar.exe');new={};count=0
    for name,items in read(MEMBERS).items():
        rows=[]
        for item in items:
            p=Path(item['path']);assert sha(p)==item['sha256']
            if p.name=='integration_CNewtonIntegration.cpp.obj':p=BUILD/'CNewtonIntegration.obj';count+=1
            rows.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in rows],cwd=str(WD));new[name]=rows
    assert count==1;save(BUILD/'members.json',new)
    subprocess.check_call(commands['helper'],cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,objectsCompiled=1,readOnlyObserverObject=True,
      matureObjectsRecompiled=0,engineSHA256=sha(BUILD/'preflight.exe')))


def main(mode):
    if CLOSED:raise SystemExit('v250 bounded NK development closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);g=setup_guard()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',DOC/'PLAN.md',
      DRIVER,HOOK,ORIGINAL,REFERENCE,FACES,MEMBERS,PREV/'audit03/report.json',PREV/'method01/report.json',
      PREV/'stage01/command.json',PREV/'stage01/dependencies.json',PREV/'native01/case.cfg',PREV/'native01/mesh.su2',
      PREV/'native01/input_restart.csv',PREV/'native01/restart_flow.csv',PREV/'native01/schema.json',
      ROOT/'outputs/v235_validation/stage01/commands.json',ROOT/'outputs/v245_validation/stage01/dependencies.json']
    if mode=='prepare':return supervise(g,BASE,mode,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,1536,locks,BASE,'offline')
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for items in read(MEMBERS).values():
        for row in items:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    locks += [STAGE/'commands.json',STAGE/'dependencies.json',STAGE/'report.json',CASE/'case.cfg',CASE/'mesh.su2',CASE/'input_restart.csv',CASE/'original_faces.txt',CASE/'schema.json']
    if mode=='build':return supervise(g,BASE,mode,[sys.executable,str(Path(__file__).resolve()),'task_build'],120,1536,list(dict.fromkeys(locks)),BASE,'build')
    if mode=='run':
        assert read(BASE/'process_build/report.json')['passed'] and read(BASE/'process_build/terminal_capacity.json')['passed']
        for items in read(BUILD/'members.json').values():
            for row in items:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
        locks += [BUILD/'preflight.exe',BUILD/'members.json']
        return supervise(g,BASE,'rocket',[str(BUILD/'preflight.exe'),'case.cfg','.',str(REFERENCE),'original_faces.txt'],1800,1536,list(dict.fromkeys(locks)),CASE,'cfd')
    raise ValueError('MODE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v250 bounded NK development closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
