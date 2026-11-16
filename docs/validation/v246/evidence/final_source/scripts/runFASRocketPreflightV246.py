"""模块职责：原真实域固定两层预检的准确输入、薄构建和唯一零流动初始化监督。"""
from pathlib import Path
import shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v246_validation';DOC=ROOT/'docs/validation/v246';CLOSED=True
PREV=ROOT/'outputs/v245_validation';GAS=ROOT/'outputs/v220_validation';OLD=ROOT/'tools/runtime/su2-thermal-build/source'
WD=OLD.parent/'build-v117';NATIVE=ROOT/'native/su2PhaseAwareAir/fasRocketPreflightV246.cpp'
BUILD=BASE/'build02';CASE=BASE/'initial01';STAGE=BASE/'stage02';MESH=ROOT/'outputs/v222_validation/rocket01/mesh.su2'


def prepare():
    """模块一：原细网格逐字复制，冻结单粗层配置；所有原数值对象复用。"""
    for p in (STAGE,BUILD,CASE):p.mkdir(exist_ok=False)
    assert sha(MESH)=='00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181'
    shutil.copyfile(MESH,CASE/'mesh.su2');cfg={}
    for line in (MESH.parent/'case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(MGLEVEL='1',MGCYCLE='V_CYCLE',MG_MIN_MESHSIZE='500',MG_PRE_SMOOTH='(1, 1)',
       MG_POST_SMOOTH='(0, 0)',MG_CORRECTION_SMOOTH='(0, 0)',MG_DAMP_RESTRICTION='0.5',MG_DAMP_PROLONGATION='0.5',
       MG_SMOOTH_EARLY_EXIT='NO',MG_SMOOTH_OUTPUT='YES',ITER='5',OUTPUT_FILES='(RESTART_ASCII)',OUTPUT_WRT_FREQ='5')
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    cmd=list(read(PREV/'stage01/commands.json')['helper']);cmd[cmd.index('-MF')+1]=str(BUILD/'preflight.d')
    cmd[cmd.index('-o')+1]=str(BUILD/'preflight.exe');cmd[cmd.index(str(ROOT/'native/su2PhaseAwareAir/fasCompatibilityV245.cpp'))]=str(NATIVE)
    save(STAGE/'command.json',cmd)
    # 依赖输出只用于准确绑定，不运行求解或重编成熟对象。
    args=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-o','-MF'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    dep=STAGE/'inputs.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=str(WD),timeout=30)
    known=read(PREV/'stage01/dependencies.json');gas=read(GAS/'build02/dependencies.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');deps={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
        if p!=NATIVE:
            expected=known.get(str(p),gas.get(str(p),tree.get(str(p))))
            if expected is None:expected=tree[str(OLD/p.relative_to(GAS/'source'))]
            assert h==expected,str(p)
        deps[str(p)]=h
    for rows in read(PREV/'build01/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256']
    save(STAGE/'dependencies.json',deps)
    src=[ROOT/'outputs/v235_validation/source01/SU2_CFD/src/drivers/CSinglezoneDriver.cpp',
      ROOT/'outputs/v235_validation/source01/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
      PREV/'stage01/CMultiGridIntegration.cpp',PREV/'stage01/CMultiGridGeometry.cpp',PREV/'stage01/CNSSolver.cpp']
    save(STAGE/'native_sources.json',{str(p):sha(p) for p in src})
    save(STAGE/'report.json',dict(passed=True,dependencies=len(deps),meshSHA256=sha(CASE/'mesh.su2'),flowUpdatesAdmitted=0,
       finePoints=335726,fineCells=705490,wallFaces=19752,productionMinimumCoarsePoints=500,physicalQualification=False))


def build():
    """模块二：仅薄入口；不构建物性、通量、边界、聚合或积分对象。"""
    cmd=read(STAGE/'command.json');subprocess.check_call(cmd,cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,thinEntryOnly=True,numericalObjectsCompiled=0,engineSHA256=sha(BUILD/'preflight.exe')))


def recover():
    """只修访问私有包装的薄入口错误；初配置和网格不变，首构建回执保留。"""
    STAGE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False);old=BASE/'stage01'
    deps=read(old/'dependencies.json');deps[str(NATIVE)]=sha(NATIVE);save(STAGE/'dependencies.json',deps)
    cmd=read(old/'command.json');cmd=[x.replace(str(BASE/'build01'),str(BUILD)) for x in cmd];save(STAGE/'command.json',cmd)
    for n in ('native_sources.json','report.json'):shutil.copyfile(old/n,STAGE/n)


def main(stage):
    if CLOSED:raise SystemExit('v246 real preflight closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',NATIVE,
        MESH,MESH.parent/'case.cfg',PREV/'stage01/commands.json',PREV/'stage01/dependencies.json',PREV/'build01/members.json']
    if stage=='prepare':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,512,locks,BASE)
    if stage=='recover':return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_recover'],30,512,locks+[BASE/'stage01/dependencies.json',BASE/'stage01/command.json'],BASE)
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for rows in read(PREV/'build01/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    locks += [STAGE/'command.json',STAGE/'dependencies.json',CASE/'mesh.su2',CASE/'case.cfg']
    if stage=='build':
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='build')
        return guard.supervise('build2',[sys.executable,str(Path(__file__).resolve()),'task_build'],60-used,512,list(dict.fromkeys(locks)),BASE,'build')
    if stage=='initialize':
        assert read(BASE/'process_build2/report.json')['passed'];locks += [BUILD/'preflight.exe']
        return guard.supervise(stage,[str(BUILD/'preflight.exe'),'case.cfg','.'],120,512,list(dict.fromkeys(locks)),CASE,'query')
    raise ValueError('TASK_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v246 real preflight closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_recover':recover()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
