"""模块职责：使用已有边界参考速度修复制造夹具，独立保存并监督整体接口控制。"""
from pathlib import Path
import shutil,subprocess,sys
import runMultiZoneCompatibilityV234 as base
from observeMultiZoneFailureV234 import safe_wait_procs
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=base.BASE;BUILD=BASE/'build03';PREP=BASE/'recovery_prepare01';DOC=base.DOC
CLOSED=True


def prepare():
    PREP.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    source=BASE/'observer01/controlObserved.cpp';text=source.read_text(encoding='utf8')
    old='double v[3]={40.+(affine?4*x[1]:0),30.+(affine?3*x[2]:0),-20.+(affine?2*x[0]:0)};'
    new='double v[3]={2*std::sqrt(ref.a2)+(affine?4*x[1]:0),.4*std::sqrt(ref.a2)+(affine?3*x[2]:0),.3*std::sqrt(ref.a2)+(affine?2*x[0]:0)};'
    assert text.count(old)==1;text=text.replace(old,new)
    base.CPP.write_text(text,encoding='utf8');shutil.copyfile(base.CPP,PREP/'control.cpp')
    for mode in ('single','multi'):
        out=BASE/(mode+'02');out.mkdir()
        for p in (BASE/(mode+'01')).iterdir():
            if p.suffix in ('.cfg','.su2') or p.name=='mesh_mapping.json':shutil.copyfile(p,out/p.name)
    cmd=read(BASE/'build01/command.json');cmd=[x.replace(str(BASE/'build01/control.d'),str(BUILD/'control.d')).replace(str(BASE/'build01/control.exe'),str(BUILD/'control.exe')) for x in cmd]
    save(BUILD/'command.json',cmd);filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    dep=PREP/'inputs.d';subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(base.OLD_BUILD),timeout=20)
    olddeps=read(BASE/'prepare01/dependencies.json');deps={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else base.OLD_BUILD/p).resolve();h=sha(p)
        if p!=base.CPP:assert olddeps[str(p)]==h,str(p)
        deps[str(p)]=h
    save(PREP/'dependencies.json',deps)
    save(PREP/'report.json',dict(passed=True,originalObserverSHA256=sha(source),newSourceSHA256=sha(base.CPP),
      referenceVelocity='v221: (2a,0.4a,0.3a), a at p1200Pa/T350.123K',originalSubsonicFailurePreserved=True,
      boundaryAndGasAndNativeObjectsUnchanged=True,manufacturedFixtureChanged=True,flowUpdates=0))


def main(stage):
    if CLOSED:raise SystemExit('v234 complete control closed; use NEXT_DECISION')
    g=base.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda _:size(BASE)
    g.psutil.wait_procs=safe_wait_procs(g.psutil.wait_procs)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runMultiZoneCompatibilityV234.py',ROOT/'scripts/observeMultiZoneFailureV234.py',
      ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage=='prepare_control':
        locks += [BASE/'observer01/controlObserved.cpp',BASE/'capture01/termination.txt',BASE/'prepare01/dependencies.json']
        return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,128,locks,BASE)
    deps=read(PREP/'dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    locks += [Path(n) for n in deps]+[BUILD/'command.json']
    for group in read(base.GAS/'build02/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    locks += [base.GAS/'build02/libSU2core.a',base.GAS/'build02/libSU2Common.a']
    if stage=='build_control':
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_build*/report.json'))
        return g.supervise(stage,read(BUILD/'command.json'),120-used,128,list(dict.fromkeys(locks)),base.OLD_BUILD,'build')
    if stage in ('single2','multi2'):
        mode=stage[:-1];out=BASE/(mode+'02');locks += [BUILD/'control.exe']+[p for p in out.iterdir() if p.is_file()]
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='recovery-query')
        assert used<60
        return g.supervise(stage,[str(BUILD/'control.exe'),'case.cfg','.',mode],60-used,128,list(dict.fromkeys(locks)),out,'recovery-query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v234 complete control closed; use NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
