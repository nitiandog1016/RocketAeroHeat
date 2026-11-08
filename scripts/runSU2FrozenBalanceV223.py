"""模块职责：唯一完整末场的只读空间算子审核；不启动流动发展或重建引擎。"""
from pathlib import Path
import json,shutil,subprocess,sys,os
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,sha,read,save
BASE=ROOT/'outputs/v223_validation';DOC=ROOT/'docs/validation/v223';BUILD=BASE/'build01';CASE=BASE/'frozen01'
CPP=ROOT/'native/su2PhaseAwareAir/frozenBalance.cpp';CLOSED=True
guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
real_size=guard.size_bytes;guard.size_bytes=lambda unused:real_size(BASE)

def command():
    old=ROOT/'outputs/v221_validation/build'
    argv=read(ROOT/'outputs/v221_validation/process_build2/report.json')['argv']
    return [str(CPP) if a==str(ROOT/'native/su2PhaseAwareAir/compatibility.cpp') else
        str(BUILD/'frozen.d') if a==str(old/'compatibility.d') else
        str(BUILD/'frozen.exe') if a==str(old/'compatibility.exe') else a for a in argv]

def prepare():
    BUILD.mkdir(exist_ok=False);CASE.mkdir(exist_ok=False)
    assert read(BASE/'inspection01/report.json')['passed'] and read(BASE/'localize01/report.json')['passed']
    origin=ROOT/'outputs/v222_validation/rocket01';cfg=(origin/'case.cfg').read_text(encoding='utf8')
    replacements={'RESTART_SOL= NO':'RESTART_SOL= YES',
        'MESH_FILENAME= mesh.su2':'MESH_FILENAME= '+(origin/'mesh.su2').as_posix(),
        'SOLUTION_FILENAME= restart_flow':'SOLUTION_FILENAME= '+(origin/'restart_flow_000199').as_posix()}
    for a,b in replacements.items():assert cfg.count(a)==1;cfg=cfg.replace(a,b)
    (CASE/'case.cfg').write_text(cfg,encoding='utf8')
    known=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in known.items():assert sha(n)==h,n
    original=command();filtered=[];i=0
    while i<len(original):
        a=original[i]
        if a in ('-MF','-o'):i+=2;continue
        if a in ('-MD','-Wl,--start-group','-Wl,--end-group') or a.endswith('.a'):i+=1;continue
        filtered.append(a);i+=1
    dep=BUILD/'inputs.d';oldbuild=ROOT/'tools/runtime/su2-thermal-build/build-v117'
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(oldbuild),timeout=20)
    names=dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split();inputs={}
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');new=ROOT/'outputs/v220_validation/source';old=ROOT/'tools/runtime/su2-thermal-build/source'
    for n in names:
        p=Path(n);p=(p if p.is_absolute() else oldbuild/p).resolve();h=sha(p);inputs[str(p)]=h
        if str(p) in known:assert h==known[str(p)]
        elif p!=CPP:
            try:rel=p.relative_to(new)
            except ValueError:continue
            assert h==tree[str(old/rel)],str(p)
    for p in (origin/'mesh.su2',origin/'restart_flow_000199.csv',CASE/'case.cfg',CPP,
              ROOT/'outputs/v220_validation/build02/libSU2core.a',ROOT/'outputs/v220_validation/build02/libSU2Common.a'):
        inputs[str(p)]=sha(p)
    assert inputs[str(origin/'restart_flow_000199.csv')]==read(ROOT/'outputs/v222_validation/audit01/report.json')['outputSHA256']['restart_flow_000199.csv']
    save(BUILD/'inputs.json',inputs);save(BUILD/'command.json',original)
    save(BUILD/'prepare_report.json',dict(passed=True,inputs=len(inputs),nativeQueries=0,flowUpdates=0))

def main(stage):
    if CLOSED:raise SystemExit('v223 frozen entry closed')
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',CPP]
    if stage=='prepare':return guard.supervise('prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],30,256,locks,BASE,'offline')
    assert read(BUILD/'prepare_report.json')['passed'];inputs=read(BUILD/'inputs.json')
    for n,h in inputs.items():assert sha(n)==h,n
    locks += [Path(n) for n in inputs]+[BUILD/'inputs.json']
    if stage=='build':return guard.supervise('build',command(),60,256,locks,BASE,'build')
    assert read(BASE/'process_build/report.json')['passed'];locks+=[BUILD/'frozen.exe']
    if stage in ('run2','query2'):
        out=BASE/'frozen02';cfg=out/'case.cfg'
        if stage=='run2':
            previous=read(BASE/'process_query/report.json')
            assert previous['exitCode']==1 and 'SetConfig_Parsing' in (BASE/'process_query/stdout.txt').read_text(encoding='utf8')
            assert not (CASE/'native.json').exists() and not (BASE/'process_query2').exists()
            out.mkdir(exist_ok=False);text=(CASE/'case.cfg').read_text(encoding='utf8')
            origin=ROOT/'outputs/v222_validation/rocket01'
            for name in ('mesh.su2','restart_flow_000199'):
                absolute=(origin/name).as_posix();relative=os.path.relpath(str(origin/name),str(out)).replace('\\','/')
                assert text.count(absolute)==1;text=text.replace(absolute,relative)
            cfg.write_text(text,encoding='utf8')
        locks += [cfg]
        if stage=='run2':return guard.supervise('parent2',[sys.executable,str(Path(__file__).resolve()),'query2'],150,256,locks,BASE,'parent')
        remaining=120-read(BASE/'process_query/report.json')['elapsedSeconds']
        return guard.supervise('query2',[str(BUILD/'frozen.exe'),str(cfg),str(out)],remaining,256,locks,out,'query',True)
    if stage=='run':
        assert not (BASE/'process_query').exists()
        return guard.supervise('parent',[sys.executable,str(Path(__file__).resolve()),'query'],150,256,locks,BASE,'parent')
    if stage=='query':return guard.supervise('query',[str(BUILD/'frozen.exe'),str(CASE/'case.cfg'),str(CASE)],120,256,locks,CASE,'query',True)
    raise ValueError('STAGE_NOT_ADMITTED')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v223 frozen entry closed')
    if sys.argv[1]=='task_prepare':prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
