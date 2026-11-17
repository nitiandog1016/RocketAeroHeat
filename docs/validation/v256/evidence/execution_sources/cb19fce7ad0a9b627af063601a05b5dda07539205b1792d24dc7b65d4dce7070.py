"""模块职责：准确原生对象复用、薄编F(x)适配器及固定完整域接口控制。"""
from pathlib import Path
import sys,os,re,shutil,subprocess
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v256_validation';DOC=ROOT/'docs/validation/v256';CASE=BASE/'rocket01';STAGE=BASE/'stage01';BUILD=BASE/'build01';CLOSED=False
PREV=ROOT/'outputs/v250_validation';WD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
DRIVER=ROOT/'native/su2PhaseAwareAir/fullDomainResidualV256.cpp';MEMBERS=PREV/'build01/members.json'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'


def prepare():
    for p in [CASE,STAGE,BUILD]:p.mkdir(exist_ok=False)
    members=read(MEMBERS)
    for rows in members.values():
        for r in rows:assert sha(r['path'])==r['sha256']
    cfg=(PREV/'rocket01/case.cfg').read_text(encoding='utf8')
    changes={'NEWTON_KRYLOV':'NO','MGLEVEL':'0','ITER':'1','MESH_FILENAME':(PREV/'rocket01/mesh.su2').as_posix(),
        'SOLUTION_FILENAME':(PREV/'rocket01/input_restart').as_posix()}
    for key,value in changes.items():cfg,n=re.subn('^'+key+r'=.*$',key+'= '+value,cfg,count=1,flags=re.M);assert n==1
    (CASE/'case.cfg').write_text(cfg,encoding='utf8')
    command=[str(DRIVER) if p==str(ROOT/'native/su2PhaseAwareAir/newtonDevelopmentV250.cpp') else p.replace(str(PREV/'build01'),str(BUILD)) for p in read(PREV/'stage01/commands.json')['helper']]
    save(STAGE/'command.json',command);args=[];i=0
    while i<len(command):
        x=command[i]
        if x in ['-o','-MF']:i+=2;continue
        if x in ['-MD','-Wl,--start-group','-Wl,--end-group'] or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    dep=STAGE/'deps.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=WD,timeout=30)
    known=read(PREV/'stage01/dependencies.json');dependencies={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
        if p!=DRIVER:assert h==known[str(p)],str(p)
        dependencies[str(p)]=h
    save(STAGE/'dependencies.json',dependencies)
    # 小证据绑定原文的明确限制，及实际成熟热壁/SA刷新调用；不重跑工程计算。
    tp=ROOT/'docs/validation/v252/evidence/method01/tp3271_equations_excerpt.txt'
    cr=ROOT/'outputs/v52_validation/research/cr112180_text.txt'
    text=cr.read_text(encoding='utf8');excerpt=text[text.index('PDF PAGE 41'):text.index('PDF PAGE 43')]
    (STAGE/'cr112180_scope_excerpt.txt').write_text(excerpt,encoding='utf8');shutil.copyfile(tp,STAGE/'tp3271_equations_excerpt.txt')
    physics=[SOURCE/'SU2_CFD/src/solvers/CTurbSASolver.cpp',SOURCE/'SU2_CFD/src/solvers/CNSSolver.cpp']
    sourceTree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    for p in physics:assert sha(p)==sourceTree[str(p)],str(p)
    save(STAGE/'method.json',dict(decision='reject present shock-mass-table-only near-wall profile reconstruction; implement mature connected-domain residual provider',
        CR112180PrintedPages=[36,37],TP3271Equations=[17,18,19],reason='shock-shape approximation explicitly not detailed near-wall layer; equation19 uses constant edge entropy',
        rawSourceSHA256={str(p):sha(p) for p in [tp,cr]+physics},KINSOLCandidateCommit='8e17876d3b4d682b4098684b07a85b005a122f81',
        KINSOLExecuted=False,KINSOLInstalled=False,externalLicense='BSD-3-Clause; LICENSE stored with source excerpts',
        originalPhysicalDomainUnchanged=True,artificialOutletCreated=False,weakEnergyWallPreserved=True))
    save(STAGE/'report.json',dict(passed=True,physicsObjectsRecompiled=0,nodes=335726,cells=705490,originalFaces=19752,
        expectedWallNodes=9878,expectedUnknowns=6*335726-4*9878,oldInitialSource='v248/200 through v250 input_restart',
        initialSourceIsNotBoundary=True,newNonlinearIterationsAdmitted=0,newTimeStepsAdmitted=0,physicalQualification=False,newQualifiedAreaM2=0))


def build():
    command=read(STAGE/'command.json');ar=str(Path(command[0]).parent/'ar.exe');members=read(MEMBERS)
    for name,rows in members.items():
        for r in rows:assert sha(r['path'])==r['sha256']
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in rows],cwd=WD)
    save(BUILD/'members.json',members);subprocess.check_call(command,cwd=WD)
    save(BUILD/'report.json',dict(passed=True,physicsObjectsRecompiled=0,adapterSHA256=sha(DRIVER),engineSHA256=sha(BUILD/'preflight.exe')))


def main(stage):
    if CLOSED:raise SystemExit('v256 fixed complete residual tasks closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original=guard.size_bytes;guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',DRIVER,ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        MEMBERS,PREV/'stage01/commands.json',PREV/'stage01/dependencies.json',PREV/'rocket01/case.cfg',PREV/'rocket01/mesh.su2',PREV/'rocket01/input_restart.csv',
        ROOT/'docs/validation/v252/evidence/method01/tp3271_equations_excerpt.txt',ROOT/'outputs/v52_validation/research/cr112180_text.txt',
        ROOT/'outputs/v167_validation/build01/source_tree.json',SOURCE/'SU2_CFD/src/solvers/CTurbSASolver.cpp',SOURCE/'SU2_CFD/src/solvers/CNSSolver.cpp']
    if stage=='prepare':argv=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];seconds=120;kind='offline';cwd=ROOT
    else:
        for p,h in read(STAGE/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
        for rows in read(MEMBERS).values():
            for r in rows:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
        locks += [STAGE/'command.json',STAGE/'dependencies.json',STAGE/'report.json',CASE/'case.cfg']
        if stage in ['build','build2']:
            if stage=='build2':assert not read(BASE/'process_build/report.json')['passed']
            argv=[sys.executable,str(Path(__file__).resolve()),'task_build'];seconds=120 if stage=='build' else 60;kind='build';cwd=ROOT
        elif stage=='control':
            assert (BUILD/'report.json').exists();locks += [BUILD/'preflight.exe',BUILD/'members.json']
            argv=[str(BUILD/'preflight.exe'),'case.cfg','.'];seconds=180;kind='control';cwd=CASE
        else:raise ValueError('STAGE_NOT_ADMITTED')
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix.lower() in ['.py','.cpp','.hpp','.h','.inl','.md','.cfg','.json','.txt','.c'] and p.stat().st_size<2*1024**2:
            target=snap/(sha(p)+p.suffix)
            if not target.exists():shutil.copyfile(p,target)
    return supervise(guard,BASE,stage,argv,seconds,512,list(dict.fromkeys(locks)),cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v256 fixed complete residual tasks closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
