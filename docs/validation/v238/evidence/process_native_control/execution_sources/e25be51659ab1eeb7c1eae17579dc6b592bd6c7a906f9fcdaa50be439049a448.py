"""模块职责：修复控制通过后，构建隔离HLLC配置并执行一次实际工厂/flow-SA整体控制。"""
from pathlib import Path
import difflib,json,os,posixpath,re,shutil,subprocess,sys,time
import runHLLCCompatibilityV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;STAGE=BASE/'native_prepare01';BUILD=BASE/'build03';CASE=BASE/'native01'
GAS=run.GAS;OLD=run.OLD;CLOSED=False


def prepare():
    assert read(BASE/'audit02/report.json')['passed']
    for p in (STAGE,BUILD,CASE):p.mkdir(exist_ok=False)
    src=GAS/'source/Common/src/CConfig.cpp';assert sha(src)==read(GAS/'build02/dependencies.json')[str(src)]
    before=src.read_text(encoding='utf8');a=before.index('  if (Kind_FluidModel == RAH_PHASE_AIR) {');b=before.index('\n  }',a)+4
    block=before[a:b];assert block.count('Kind_Upwind_Flow != UPWIND::ROE')==1
    block=block.replace('Kind_Upwind_Flow != UPWIND::ROE','Kind_Upwind_Flow != UPWIND::HLLC').replace('NS/SA, ROE, resolved wall','NS/SA, HLLC-V238, resolved wall')
    text=before[:a]+block+before[b:];rel='Common/src/CConfig.cpp'
    text=re.sub(r'#include "(\.\./[^\"]+)"',lambda m:'#include "'+(GAS/'source'/posixpath.normpath(posixpath.join(posixpath.dirname(rel),m.group(1)))).as_posix()+'"',text)
    target=STAGE/'CConfig.cpp';target.write_text(text,encoding='utf8')
    (STAGE/'config.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=rel,tofile='hllc-v238/'+rel)),encoding='utf8')
    original=ROOT/'native/su2PhaseAwareAir/singleDomainRunV235.cpp';before=original.read_text(encoding='utf8');text=before
    assert sha(original)==read(ROOT/'outputs/v237_validation/stage02/inputs.json')[str(original)]
    text='#include "numerics/flow/convection/hllc.hpp"\n#include "PhaseStepAcceptanceV237.hpp"\n#include "output/COutput.hpp"\n'+text
    text=text.replace('SingleDomainDriverV235','HLLCDriverV238')
    anchor='   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");'
    assert text.count(anchor)==1
    text=text.replace(anchor,'''   // 检查实际工厂及同一气体/SA组合，不靠配置字符串推断派发。
   auto*** num=numerics_container[0][0][0];
   if(cfg->GetKind_Upwind_Flow()!=UPWIND::HLLC||geo->GetnPoint()!=24||cfg->GetnInner_Iter()!=5||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_TERM])||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_BOUND_TERM]))
      throw std::runtime_error("HLLC_V238_NATIVE_DISPATCH");
   std::ofstream identity(folder+"/native_identity.json");
   identity<<"{\\\"schemeID\\\":\\\"phase-air-general-hllc-density-jacobian-fix-v238\\\",\\\"physicalGasID\\\":15,\\\"actualInteriorBoundaryHLLC\\\":true,\\\"sameSAAirIdentity\\\":true,\\\"physicalQualification\\\":false}\\n";identity.close();
'''+anchor)
    anchor='   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;'
    assert text.count(anchor)==1
    text=text.replace(anchor,'''   // 显式原生终点保存，不设置或伪造收敛标志。
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced||completed!=5||rah_phase_step_v237::checked!=5*24)throw std::runtime_error("HLLC_V238_COMPLETE_CONTROL");
   std::ofstream acceptance(folder+"/acceptance_summary.json");
   acceptance<<"{\\\"checked\\\":"<<rah_phase_step_v237::checked<<",\\\"limited\\\":"<<rah_phase_step_v237::limited<<",\\\"forcedFinalNativeOutput\\\":true}\\n";acceptance.close();
'''+anchor)
    driver=STAGE/'nativeHLLCV238.cpp';driver.write_text(text,encoding='utf8')
    (STAGE/'driver.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=original.name,tofile=driver.name)),encoding='utf8')
    command=list(read(GAS/'prepare03/commands.json')[rel]);command[command.index('-c')+1]=str(target);command[command.index('-o')+1]=str(BUILD/'CConfig.obj')
    helper=read(ROOT/'outputs/v237_validation/stage02/commands.json')['driver']
    helper=[str(driver) if x.endswith('restartBranchV237.cpp') else x.replace(str(ROOT/'outputs/v237_validation/build01'),str(BUILD)) if x.endswith(('.a','.exe','.d')) else x for x in helper]
    commands=dict(config=command,helper=helper);save(STAGE/'commands.json',commands)
    deps={};known=read(ROOT/'outputs/v237_validation/stage02/inputs.json');known.update(read(BASE/'prepare01/inputs.json'));tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    for name,cmd in commands.items():
        filtered=[];i=0
        while i<len(cmd):
            x=cmd[i]
            if x in ('-o','-MF'):i+=2;continue
            if x in ('-c','-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            filtered.append(x);i+=1
        dep=STAGE/(name+'.d');subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(run.OLDBUILD),timeout=20)
        for raw in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(raw);p=(p if p.is_absolute() else run.OLDBUILD/p).resolve();h=sha(p)
            if p not in (target,driver):
                expected=known.get(str(p),tree.get(str(p)))
                if expected is None:expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/p.relative_to(GAS/'source'))]
                assert h==expected,str(p)
            deps[str(p)]=h
    for rows in read(BASE/'build02/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'];deps[r['path']]=r['sha256']
    guarded=ROOT/'outputs/v237_validation/build01/CEulerSolver.obj'
    record=next(r for rows in read(ROOT/'outputs/v237_validation/build01/members.json').values() for r in rows if r['path']==str(guarded))
    assert sha(guarded)==record['sha256'];deps[str(guarded)]=sha(guarded)
    cfg=(OLD/'observe02/case.cfg').read_text(encoding='utf8');assert cfg.count('CONV_NUM_METHOD_FLOW= ROE')==1
    cfg=cfg.replace('CONV_NUM_METHOD_FLOW= ROE','CONV_NUM_METHOD_FLOW= HLLC')
    # 不复制旧小网格；配置引用原准确文件。
    lines=[]
    for line in cfg.splitlines():
        if line.startswith('MESH_FILENAME='):
            raw=line.split('=',1)[1].strip();p=(OLD/'observe02'/raw).resolve();deps[str(p)]=sha(p)
            line='MESH_FILENAME= '+os.path.relpath(str(p),str(CASE)).replace('\\','/')
        lines.append(line)
    (CASE/'case.cfg').write_text('\n'.join(lines)+'\n',encoding='utf8');deps[str(CASE/'case.cfg')]=sha(CASE/'case.cfg')
    save(STAGE/'inputs.json',deps);save(STAGE/'report.json',dict(passed=True,schemeID='phase-air-general-hllc-density-jacobian-fix-v238',physicalGasIdentity=15,
        oldROEPackUnchanged=True,configObjectsToCompile=1,hllcAndGuardObjectsReused=True,realCFDAdmitted=False))


def build():
    cmd=read(STAGE/'commands.json');t=time.monotonic();subprocess.check_call(cmd['config'],cwd=str(run.OLDBUILD))
    save(BUILD/'compile.json',dict(configSeconds=time.monotonic()-t));ar=str(Path(cmd['config'][0]).parent/'ar.exe');members={};changed=0
    for name,rows in read(BASE/'build02/members.json').items():
        result=[]
        for r in rows:
            p=Path(r['path']);assert sha(p)==r['sha256']
            if 'CConfig' in p.name and p.suffix=='.obj':p=BUILD/'CConfig.obj';changed+=1
            if p.name=='CEulerSolver.obj':p=ROOT/'outputs/v237_validation/build01/CEulerSolver.obj';changed+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result]);members[name]=result
    assert changed==2,changed;save(BUILD/'members.json',members)
    subprocess.check_call(cmd['helper'],cwd=str(run.OLDBUILD));save(BUILD/'report.json',dict(passed=True,driverSHA256=sha(BUILD/'driver.exe'),schemeID='phase-air-general-hllc-density-jacobian-fix-v238'))


def main(stage):
    if CLOSED:raise SystemExit('v238 native execution closed; follow NEXT_DECISION')
    run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',BASE/'audit02/report.json']
    if stage=='prepare':return run.guard.supervise('native_prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,1024,locks,BASE)
    assert read(STAGE/'report.json')['passed']
    for p,h in read(STAGE/'inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
    locks+=[STAGE/'inputs.json',STAGE/'commands.json']
    if stage=='build':
        remaining=180-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        return run.guard.supervise('native_build',[sys.executable,str(Path(__file__).resolve()),'task_build'],remaining,1024,locks,run.OLDBUILD,'build')
    if stage=='control':
        assert read(BUILD/'report.json')['passed'];locks.append(BUILD/'driver.exe')
        remaining=60-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='query')
        return run.guard.supervise('native_control',[str(BUILD/'driver.exe'),'case.cfg','.','observe','5'],remaining,1024,locks,CASE,'query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 native execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
