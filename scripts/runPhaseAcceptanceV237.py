"""模块职责：隔离相态更新守卫、回放真实事件，控制唯一真实候选准入；不修改共享求解器。"""
from pathlib import Path
import difflib,json,os,re,shutil,subprocess,sys,time,posixpath
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v237_validation';DOC=ROOT/'docs/validation/v237';STAGE=BASE/'stage01';BUILD=BASE/'build01'
NATIVE=ROOT/'native/su2PhaseAwareAir';OLD=ROOT/'outputs/v235_validation';GAS=ROOT/'outputs/v220_validation'
OLDBUILD=ROOT/'tools/runtime/su2-thermal-build/build-v117';CLOSED=True


def configure():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes
    guard.size_bytes=lambda p:original(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0


def prepare():
    """模块一：准确旧源码只插入相态域守卫，继承原生松弛上限；生成明确补丁与依赖。"""
    STAGE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=True);(BASE/'control01').mkdir(exist_ok=True)
    src=OLD/'source01/SU2_CFD/src/solvers/CEulerSolver.cpp';known=read(OLD/'stage01/dependencies.json')
    assert sha(src)==read(GAS/'build02/dependencies.json')[str(GAS/'source/SU2_CFD/src/solvers/CEulerSolver.cpp')]
    before=src.read_text(encoding='utf8');text=before;rel='SU2_CFD/src/solvers/CEulerSolver.cpp'
    def absolute(m):return '#include "'+(GAS/'source'/posixpath.normpath(posixpath.join(posixpath.dirname(rel),m.group(1)))).as_posix()+'"'
    text=re.sub(r'#include "(\.\./[^\"]+)"',absolute,text)
    text='#include "'+(NATIVE/'PhaseStepAcceptanceV237.hpp').as_posix()+'"\n'+text
    anchor='    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);';assert text.count(anchor)==1
    text=text.replace(anchor,'''    // 相态守卫只缩减原生松弛步长，五量共用同一个系数，不回写或裁剪状态。
    if (config->GetKind_FluidModel() == RAH_PHASE_AIR) {
      if (nDim != 3 || nVar != 5 || omp_get_max_threads() != 1)
        throw std::runtime_error("PHASE_STEP_NATIVE_IDENTITY");
      localUnderRelaxation = rah_phase_step_v237::native_limit(config->GetInnerIter()+1, iPoint,
          nodes->GetSolution(iPoint), LinSysSol.GetBlock(iPoint), localUnderRelaxation);
    }
'''+anchor)
    target=STAGE/'CEulerSolver.cpp';target.write_text(text,encoding='utf8')
    (STAGE/'native_guard.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),text.splitlines(True),fromfile=rel,tofile='phase-acceptance/'+rel)),encoding='utf8')
    original=ROOT/'native/su2PhaseAwareAir/restartBranchV236.cpp';driver=original.read_text(encoding='utf8')
    assert sha(original)==read(ROOT/'outputs/v236_validation/build01/inputs.json')[str(original)]
    driver=driver.replace('V236','V237').replace('branch236','branch237').replace('reset-v236','reset-v237')
    driver='#include "'+(NATIVE/'PhaseStepAcceptanceV237.hpp').as_posix()+'"\n'+driver
    insert='''
// 模块四：只记录实际受限更新；包括真实线性增量、原生上限和全局ID。
namespace accept237 {
std::ofstream log;
void record(unsigned long update,unsigned long point,const rah_phase_step_v237::State& old,
 const rah_phase_step_v237::State& delta,double upper,const rah_phase_step_v237::Result& r){
 if(rah_phase_step_v237::limited>100000)throw std::runtime_error("PHASE_STEP_EVENT_LIMIT");
 log<<update<<','<<point<<','<<observe::geo->nodes->GetGlobalIndex(point)<<','<<upper<<','<<r.alpha<<','<<r.reductions;
 for(const auto* a:{&old,&delta,&r.accepted})for(double v:*a)log<<','<<v;log<<'\\n';log.flush();
 if(!log)throw std::runtime_error("PHASE_STEP_EVENT_WRITE");
}
}
'''
    driver=driver.replace('int main(int argc,char** argv){try{',insert+'\nint main(int argc,char** argv){try{')
    anchor='SU2_MPI::Init(&argc,&argv);BranchDriverV237 driver'
    assert driver.count(anchor)==1
    driver=driver.replace(anchor,'''accept237::log.open(observe::folder+"/accepted_steps.csv");accept237::log<<std::setprecision(17)<<"update,localID,globalID,nativeUpper,acceptedAlpha,reductions";
 for(auto prefix:{"old","delta","accepted"})for(int j=0;j<5;++j)accept237::log<<','<<prefix<<j;accept237::log<<'\\n';
 rah_phase_step_v237::observer=accept237::record;
 '''+anchor)
    anchor='   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;'
    assert driver.count(anchor)==1
    driver=driver.replace(anchor,'''   // 收益停止不能依赖收敛标志；调用成熟强制输出接口，保存原生最终检查点。
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced)throw std::runtime_error("FORCED_FINAL_OUTPUT_MISSING");
   std::ofstream acceptReport(folder+"/acceptance_summary.json");acceptReport<<"{\\\"checked\\\":"<<rah_phase_step_v237::checked<<",\\\"limited\\\":"<<rah_phase_step_v237::limited<<",\\\"forcedFinalNativeOutput\\\":true}\\n";
'''+anchor)
    (STAGE/'restartBranchV237.cpp').write_text(driver,encoding='utf8')
    (STAGE/'driver.patch').write_text(''.join(difflib.unified_diff(original.read_text(encoding='utf8').splitlines(True),driver.splitlines(True),fromfile=original.name,tofile='restartBranchV237.cpp')),encoding='utf8')
    commands={};cmd=read(GAS/'prepare03/commands.json')[rel];cmd=list(cmd);cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(BUILD/'CEulerSolver.obj')
    commands['euler']=cmd+['-I'+str(GAS/'source/SU2_CFD/include'),'-MD','-MF',str(BUILD/'CEulerSolver.obj.d')]
    for name,cpp in [('control',NATIVE/'phaseStepControlV237.cpp'),('driver',STAGE/'restartBranchV237.cpp')]:
        cmd=read(ROOT/'outputs/v236_validation/build01/command.json');cmd=[str(cpp) if x.endswith('restartBranchV236.cpp') else
          str(BUILD/(name+'.d')) if x.endswith('branch.d') else str(BUILD/(name+'.exe')) if x.endswith('branch.exe') else
          x.replace(str(OLD/'build01'),str(BUILD)) if x.endswith('.a') else x for x in cmd];commands[name]=cmd+['-I'+str(NATIVE)]
    save(STAGE/'commands.json',commands)
    dependencies()


def dependencies():
    """恢复准确头来源映射，不覆盖已生成源码、首失败或命令。"""
    commands=read(STAGE/'commands.json');known=read(GAS/'build02/dependencies.json')
    if '-I'+str(NATIVE) not in commands['driver']:
        shutil.copyfile(STAGE/'commands.json',STAGE/'commands.beforeIncludeFix.json')
        for name in ('driver','control'):commands[name].append('-I'+str(NATIVE))
        save(STAGE/'commands.json',commands)
    deps={};tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');previous=read(ROOT/'outputs/v236_validation/build01/inputs.json')
    for name,cmd in commands.items():
        filtered=[];i=0
        while i<len(cmd):
            x=cmd[i]
            if x in ('-MF','-o'):i+=2;continue
            if x in ('-MD','-c','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
            filtered.append(x);i+=1
        dep=STAGE/(name+'.d');subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(OLDBUILD),timeout=20)
        for raw in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(raw);p=(p if p.is_absolute() else OLDBUILD/p).resolve();h=sha(p)
            if p.parent not in (STAGE,NATIVE):
                expected=previous.get(str(p),known.get(str(p),tree.get(str(p))))
                if expected is None:
                    rel=p.relative_to(GAS/'source');expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/rel)]
                assert h==expected,str(p)
            elif p.name in ('PhaseAwareAir.hpp','PhaseAirTable.hpp'):assert h==sha(GAS/'source/SU2_CFD/include/fluid'/p.name)
            deps[str(p)]=h
    for rows in read(OLD/'build01/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'];deps[r['path']]=r['sha256']
    save(STAGE/'inputs.json',deps)
    cfg=ROOT/'outputs/v235_validation/observe02/case.cfg';shutil.copyfile(cfg,BASE/'control01/case.cfg')
    save(STAGE/'report.json',dict(passed=True,sourceSHA256=sha(OLD/'source01/SU2_CFD/src/solvers/CEulerSolver.cpp'),oldDriverSHA256=sha(NATIVE/'restartBranchV236.cpp'),boundInputs=len(deps),numericalObjectsPlanned=1,
       guard='feasibility-only halving within native local relaxation',newNativeCalls=0,physicalQualification=False))


def build():
    commands=read(STAGE/'commands.json');receipts=[]
    t=time.monotonic();subprocess.check_call(commands['euler'],cwd=str(OLDBUILD));receipts.append(dict(unit='euler',seconds=time.monotonic()-t))
    save(BUILD/'compile_receipts.json',receipts);ar=str(Path(commands['euler'][0]).parent/'ar.exe');members={};changed=0
    for name,rows in read(OLD/'build01/members.json').items():
        result=[]
        for row in rows:
            p=Path(row['path']);assert sha(p)==row['sha256']
            if p.name=='CEulerSolver.obj':p=BUILD/'CEulerSolver.obj';changed+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result],cwd=str(OLDBUILD));members[name]=result
    assert changed==1;save(BUILD/'members.json',members)
    for name in ('control','driver'):
        t=time.monotonic();subprocess.check_call(commands[name],cwd=str(OLDBUILD));receipts.append(dict(unit=name,seconds=time.monotonic()-t));save(BUILD/'compile_receipts.json',receipts)
    save(BUILD/'report.json',dict(passed=True,changedNativeObjects=1,controlSHA256=sha(BUILD/'control.exe'),driverSHA256=sha(BUILD/'driver.exe')))


def repair():
    """仅修复控制薄驱动的const指针列表；原成功数值对象不重编，历史输入保持。"""
    old=BASE/'source_history/phaseStepControlV237.first.cpp';now=NATIVE/'phaseStepControlV237.cpp'
    assert now.read_text(encoding='utf8')==old.read_text(encoding='utf8').replace(
        'for(const auto* values:{&old,&delta,&accepted.accepted})',
        'for(const State* values:std::array<const State*,3>{&old,&delta,&accepted.accepted})')
    bound=read(STAGE/'inputs.beforeConstPointerFix.json');assert bound[str(now)]==sha(old)
    for name,h in bound.items():
        if name!=str(now):assert sha(name)==h,name
    bound[str(now)]=sha(now);save(STAGE/'inputs.json',bound)
    save(STAGE/'repair.json',dict(onlyConstPointerListChanged=True,previousSHA256=sha(old),currentSHA256=sha(now)))


def build2():
    """复用首构建已经成功的CEuler对象和薄归档，只恢复两个薄驱动。"""
    for archive,rows in read(BUILD/'members.json').items():
        assert (BUILD/archive).is_file()
        for row in rows:assert sha(row['path'])==row['sha256'],row['path']
    receipts=read(BUILD/'compile_receipts.first.json');assert len(receipts)==1 and receipts[0]['unit']=='euler'
    commands=read(STAGE/'commands.json')
    for name in ('control','driver'):
        t=time.monotonic();subprocess.check_call(commands[name],cwd=str(OLDBUILD));receipts.append(dict(unit=name,seconds=time.monotonic()-t));save(BUILD/'compile_receipts.json',receipts)
    save(BUILD/'report.json',dict(passed=True,changedNativeObjects=1,eulerObjectReusedAfterThinCompileFailure=True,
        controlSHA256=sha(BUILD/'control.exe'),driverSHA256=sha(BUILD/'driver.exe')))


def main(stage):
    global STAGE
    if CLOSED:raise SystemExit('v237 execution closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);configure()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',NATIVE/'PhaseStepAcceptanceV237.hpp',NATIVE/'phaseStepControlV237.cpp']
    if stage in ('prepare','prepare2'):return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],120,1024,locks,BASE)
    STAGE=BASE/'stage02'
    if stage=='repair':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_repair'],30,1024,locks+[STAGE/'inputs.beforeConstPointerFix.json'],BASE)
    if stage in ('dependencies','dependencies2'):return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_dependencies'],120,1024,locks+[STAGE/'CEulerSolver.cpp',STAGE/'restartBranchV237.cpp'],BASE)
    assert read(STAGE/'report.json')['passed']
    for n,h in read(STAGE/'inputs.json').items():assert sha(n)==h,n;locks.append(Path(n))
    locks += [STAGE/'commands.json',STAGE/'inputs.json']
    if stage=='build':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build'],180,1024,locks,BASE,'build')
    if stage=='build2':
        remaining=180-read(BASE/'process_build/report.json')['elapsedSeconds']
        return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build2'],remaining,1024,locks+[BUILD/'members.json'],BASE,'build')
    if stage=='control':
        assert read(BUILD/'report.json')['passed']
        case=BASE/'control01';events=ROOT/'outputs/v236_validation/rocket01/rejections.csv'
        locks += [BUILD/'control.exe',case/'case.cfg',events]
        return guard.supervise(stage,[str(BUILD/'control.exe'),'case.cfg',str(events),'.'],60,1024,locks,case,'query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v237 execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_prepare2':STAGE=BASE/'stage02';prepare()
    elif sys.argv[1]=='task_dependencies':STAGE=BASE/'stage02';dependencies()
    elif sys.argv[1]=='task_build':STAGE=BASE/'stage02';build()
    elif sys.argv[1]=='task_repair':STAGE=BASE/'stage02';repair()
    elif sys.argv[1]=='task_build2':STAGE=BASE/'stage02';build2()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
