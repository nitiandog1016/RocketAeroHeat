"""模块职责：原失败输入的唯一具名异常/阶段观察，保留首源码/二进制及缺失回执。"""
from pathlib import Path
import json,shutil,sys
import runMultiZoneCompatibilityV234 as base
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=base.BASE;OBS=BASE/'observer01';BUILD=BASE/'build02';DOC=base.DOC


def prepare():
    OBS.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    original=base.CPP.read_text(encoding='utf8');assert sha(base.CPP)==read(BASE/'process_build/admission.json')['inputSHA256'][str(base.CPP)]
    shutil.copyfile(base.CPP,OBS/'originalHelper.cpp');code=original
    header='''
// 仅记录未捕获异常和确切调用阶段，禁止将异常退出当作成功。
namespace failure_observer {
const char* phase="constructor";int scenario=-1,zone=-1;
void failed(){std::string message="unknown";try{auto p=std::current_exception();if(p)std::rethrow_exception(p);}
catch(const std::exception& e){message=e.what();}catch(...){message="non-standard exception";}
std::ofstream f("termination.txt");f<<"phase="<<phase<<"\\nscenario="<<scenario<<"\\nzone="<<zone<<"\\nwhat="<<message<<"\\n";f.flush();std::_Exit(86);}
}
'''
    code=code.replace('struct Zone {',header+'\nstruct Zone {',1)
    code=code.replace('SU2_MPI::Init(&argc,&argv);if(argc!=4)', 'std::set_terminate(failure_observer::failed);SU2_MPI::Init(&argc,&argv);if(argc!=4)',1)
    code=code.replace('for(int scenario=0;scenario<4;++scenario){','for(int scenario=0;scenario<4;++scenario){failure_observer::scenario=scenario;failure_observer::phase="install";',1)
    code=code.replace('driver.exchange();for(unsigned z=0;z<nz;++z)gradients', 'failure_observer::phase="exchange-and-gradient";driver.exchange();for(unsigned z=0;z<nz;++z)gradients',1)
    code=code.replace('auto z=driver.zone(iz);auto* fn=', 'failure_observer::zone=iz;failure_observer::phase="dump-state";auto z=driver.zone(iz);auto* fn=',1)
    for needle,phase in [
      ('clear(z);auto** nf=', 'interface'),
      ('clear(z);z.flow->Upwind_Residual', 'interior-flow'),
      ('z.sa->Upwind_Residual', 'interior-SA'),
      ('z.flow->BC_Riemann', 'farfield-flow'),('z.sa->BC_Riemann','farfield-SA'),
      ('z.flow->BC_Isothermal_Wall','thermal-wall-flow'),('z.sa->BC_Isothermal_Wall','thermal-wall-SA'),
      ('z.flow->Pressure_Forces','surface-loads')]:
        assert code.count(needle)==1,needle;code=code.replace(needle,'failure_observer::phase="'+phase+'";'+needle)
    source=OBS/'controlObserved.cpp';source.write_text(code,encoding='utf8')
    cmd=read(BASE/'build01/command.json');cmd=[str(source) if x==str(base.CPP) else x.replace(str(BASE/'build01/control.d'),str(BUILD/'control.d')).replace(str(BASE/'build01/control.exe'),str(BUILD/'control.exe')) for x in cmd]
    save(BUILD/'command.json',cmd)
    out=BASE/'capture01';out.mkdir()
    for n in ('mesh.su2','case.cfg'):shutil.copyfile(BASE/'single01'/n,out/n)
    save(OBS/'report.json',dict(passed=True,originalSourceSHA256=sha(OBS/'originalHelper.cpp'),observedSourceSHA256=sha(source),
      originalBinarySHA256=sha(BASE/'build01/control.exe'),originalNativeReceiptMissing=True,firstQueryBudgetReservedSeconds=60,
      configurationSHA256=sha(out/'case.cfg'),manufacturedStatesUnchanged=True,numericalMethodsUnchanged=True))


def safe_wait_procs(original):
    def wrapped(procs,timeout=None,callback=None):
        try:return original(procs,timeout=timeout,callback=callback)
        except base.guard.psutil.AccessDenied:
            gone=[];alive=[]
            for p in procs:
                try:
                    if p.is_running():alive.append(p)
                    else:gone.append(p)
                except base.guard.psutil.NoSuchProcess:gone.append(p)
                except base.guard.psutil.AccessDenied:alive.append(p)
            return gone,alive
    return wrapped


def main(stage):
    g=base.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda _:size(BASE)
    g.psutil.wait_procs=safe_wait_procs(g.psutil.wait_procs)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runMultiZoneCompatibilityV234.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage=='prepare_observer':return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,128,locks,BASE)
    deps=read(BASE/'prepare01/dependencies.json');locks += [Path(n) for n in deps]+[OBS/'controlObserved.cpp',BUILD/'command.json']
    for group in read(base.GAS/'build02/members.json').values():
        for row in group:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    locks += [base.GAS/'build02/libSU2core.a',base.GAS/'build02/libSU2Common.a']
    if stage=='build_observer':
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_build*/report.json'))
        return g.supervise(stage,read(BUILD/'command.json'),120-used,128,list(dict.fromkeys(locks)),base.OLD_BUILD,'build')
    if stage=='capture':
        out=BASE/'capture01';locks += [BUILD/'control.exe',out/'case.cfg',out/'mesh.su2']
        return g.supervise(stage,[str(BUILD/'control.exe'),'case.cfg','.','single'],20,128,list(dict.fromkeys(locks)),out,'recovery-query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
