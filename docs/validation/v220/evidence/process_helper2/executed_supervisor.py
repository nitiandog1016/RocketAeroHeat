"""模块职责：按v220事前分组预算监督隔离准备/编译/控制，不含CFD入口。"""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time
import psutil
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v220_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf8')

def supervise(stage,argv,kind,limit,output_mib,locks,cwd):
    folder=BASE/('process_'+stage);folder.mkdir(parents=True,exist_ok=False)
    records=[read(p) for p in BASE.glob('process_*/report.json')]
    total={'offline':120,'build':480,'query':60}[kind]
    used=sum(r['elapsedSeconds'] for r in records if r['budgetKind']==kind)
    limit=min(limit,total-used);assert limit>0
    before={str(p):sha(p) for p in locks}
    shutil.copyfile(__file__,folder/'executed_supervisor.py')
    if stage in ('prepare','prepare2','prepare3'):shutil.copyfile(ROOT/'scripts/prepareSU2PhaseAirV220.py',folder/'executed_task.py')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf8')
    start=time.monotonic();peak=0;vpeak=0;family={};stop=None
    save(folder/'admission.json',dict(argv=argv,budgetKind=kind,limitSeconds=limit,memoryBytes=4*1024**3,
        outputMiB=output_mib,inputSHA256=before,cwd=str(cwd),hardAddressSpaceLimit=False))
    with (folder/'stdout.txt').open('wb') as log:
        process=subprocess.Popen(argv,cwd=str(cwd),env=env,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        parent=psutil.Process(process.pid);parent.cpu_affinity(parent.cpu_affinity()[:1]);family[parent.pid]=parent
        try:
            while process.poll() is None:
                try:
                    for child in parent.children(recursive=True):family[child.pid]=child
                except psutil.NoSuchProcess:pass
                rss=0;vms=0
                for child in family.values():
                    try:mi=child.memory_info();rss+=mi.rss;vms+=mi.vms
                    except psutil.NoSuchProcess:pass
                peak=max(peak,rss);vpeak=max(vpeak,vms)
                size=0
                for top in BASE.iterdir():
                    if top.name=='source':
                        if kind=='offline':
                            for name in ('Common','SU2_CFD'):
                                part=top/name
                                if part.exists():size+=sum(p.stat().st_size for p in part.rglob('*') if p.is_file())
                        continue
                    if kind=='query' and top.name=='build':continue
                    if top.is_dir():size+=sum(p.stat().st_size for p in top.rglob('*') if p.is_file())
                    elif top.is_file():size+=top.stat().st_size
                if time.monotonic()-start>limit:stop='time';break
                if max(rss,vms)>4*1024**3:stop='memory';break
                if size>output_mib*1024**2:stop='output';break
                time.sleep(.05)
        finally:
            if process.poll() is None:
                for child in reversed(list(family.values())):
                    try:child.kill()
                    except psutil.NoSuchProcess:pass
            code=process.wait(timeout=5);_,alive=psutil.wait_procs(list(family.values()),timeout=5)
    stable=all(sha(Path(n))==h for n,h in before.items())
    report=dict(passed=code==0 and stop is None and stable and not alive,status='terminal',budgetKind=kind,
        argv=argv,exitCode=code,elapsedSeconds=time.monotonic()-start,limitSeconds=limit,stopReason=stop,
        peakObservedGroupRSSBytes=peak,peakObservedGroupVirtualBytes=vpeak,outputMiB=output_mib,
        inputSHA256=before,inputsStable=stable,liveGroupAfterCleanup=[p.pid for p in alive],cpuLimit=1,
        hardAddressSpaceLimit=False,CFDRuns=0)
    save(folder/'report.json',report)
    print(json.dumps({k:report[k] for k in ('passed','elapsedSeconds','peakObservedGroupRSSBytes','exitCode','stopReason')}))
    return report['passed']

def main(stage):
    BASE.mkdir(exist_ok=True)
    locks=[Path(__file__).resolve(),ROOT/'docs/validation/v220/PLAN.md',ROOT/'scripts/prepareSU2PhaseAirV220.py',
        ROOT/'native/su2PhaseAwareAir/PhaseAwareAir.hpp',ROOT/'native/su2PhaseAwareAir/CSU2PhaseAwareAir.hpp',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'outputs/v167_validation/build01/source_tree.json',
        ROOT/'outputs/v167_validation/build01/objects_after.json']
    if stage in ('prepare','prepare2','prepare3'):
        return supervise(stage,[sys.executable,str(ROOT/'scripts/prepareSU2PhaseAirV220.py')],'offline',60,512,locks,BASE)
    from buildSU2PhaseAirV220 import compile_commands
    assert read(BASE/'prepare03/report.json')['passed']
    locks += [ROOT/'scripts/buildSU2PhaseAirV220.py',BASE/'prepare03/commands.json',BASE/'prepare03/changes.json']
    if stage=='dependencies':
        return supervise(stage,[sys.executable,str(ROOT/'scripts/buildSU2PhaseAirV220.py'),'dependencies'],'build',60,256,locks,BASE)
    if stage.startswith('compile') and stage[-1:] in '01234' and len(stage)==8:
        inputs=read(BASE/'dependencies01/inputs.json')
        assert all(sha(Path(n))==h for n,h in inputs.items())
        locks += [Path(n) for n in inputs]
        return supervise(stage,compile_commands()[int(stage[-1])],'build',120,256,locks,BASE)
    if stage=='link':
        for i in range(5):assert read(BASE/('process_compile%d/report.json'%i))['passed']
        return supervise(stage,[sys.executable,str(ROOT/'scripts/buildSU2PhaseAirV220.py'),'link'],'build',30,256,locks,BASE)
    if stage in ('helper','query','helper2','query2'):
        from prepareSU2PhaseAirV220 import SOURCE,BUILD,OLD_BUILD,BOOST
        assert read(BASE/'link01/report.json')['passed']
        dependencies=read(BASE/'dependencies01/inputs.json')
        assert all(sha(Path(n))==h for n,h in dependencies.items())
        locks += [Path(n) for n in dependencies]+[ROOT/'native/su2PhaseAwareAir/probe.cpp']
        suffix='02' if stage.endswith('2') else ''
        if stage.startswith('helper'):
            compiler=compile_commands()[0][0]
            command=[compiler,'-std=c++17','-O2','-fopenmp','-DNDEBUG','-march=x86-64','-ffast-math',
                '-fno-finite-math-only','-DEIGEN_DONT_PARALLELIZE','-I'+str(SOURCE/'SU2_CFD/include'),
                '-I'+str(SOURCE/'SU2_CFD/include/fluid'),'-I'+str(OLD_BUILD/'Common/src'),
                '-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),'-I'+str(BOOST),
                '-MD','-MF',str(BUILD/('probe'+suffix+'.d')),str(ROOT/'native/su2PhaseAwareAir/probe.cpp'),
                '-Wl,--start-group',str(BUILD/'libSU2core.a'),str(BUILD/'libSU2Common.a'),'-Wl,--end-group',
                '-o',str(BUILD/('phase_probe'+suffix+'.exe'))]
            return supervise(stage,command,'build',120,256,locks,BASE)
        assert read(BASE/('process_helper2' if suffix else 'process_helper')/'report.json')['passed']
        out=BASE/('control02' if suffix else 'control01');out.mkdir(exist_ok=False)
        locks += [BUILD/('phase_probe'+suffix+'.exe')]
        return supervise(stage,[str(BUILD/('phase_probe'+suffix+'.exe')),str(out)],'query',30,32,locks,out)
    raise ValueError('Stage not yet admitted')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
