"""模块职责：准确绑定原成熟对象，监督唯一200步真实FAS发展；不重启旧场或扫描参数。"""
from pathlib import Path
import os,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v247_validation';DOC=ROOT/'docs/validation/v247';CASE=BASE/'rocket01'
STAGE=BASE/'stage01';BUILD=BASE/'build01';PREV=ROOT/'outputs/v246_validation'
NATIVE=ROOT/'native/su2PhaseAwareAir/fasDevelopmentV247.cpp'
WD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
MESH=ROOT/'outputs/v222_validation/rocket01/mesh.su2'
UNIFORM=PREV/'control_observe01/state_000000.bin'
FACES=ROOT/'outputs/v236_validation/rocket01/original_faces.txt'
CLOSED=False


def prepare():
    """模块一：冻结原网格/初场身份，生成显式数组语义和薄入口构建命令。"""
    assert read(PREV/'control_audit02/report.json')['completeNativeFiveStepCompatibilityPassed']
    for p in (CASE,STAGE,BUILD):p.mkdir(exist_ok=False)
    assert sha(MESH)=='00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181'
    os.link(str(MESH),str(CASE/'mesh.su2'));shutil.copyfile(FACES,CASE/'original_faces.txt')
    cfg={}
    for line in (PREV/'control_observe01/case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(ITER='200',OUTPUT_WRT_FREQ='200',SCREEN_WRT_FREQ_INNER='10')
    (CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    columns='globalID x y z volume rho rhou rhov rhow rhoE_internal nuTilde T u v w p rho_primitive H_internal a mu mut_flow k Cp mut_SA wallDistance flowRelaxation flowCFL fasCorrection0 fasCorrection1 fasCorrection2 fasCorrection3 fasCorrection4_internal lastFineSAImplicitRHS reserved0'.split()
    save(CASE/'schema.json',dict(format='RAH247LE',columns=columns,initialWorkArraysUnknown=True,internalEnergyOffsetJkg=300000,
        primitiveStage='fine completed flow, prior to latest SA mut refresh',fineWorkArray='FAS prolongation correction, NOT flow residual',
        residualTrend='stages.csv fine after_time: pre-smoothing before coarse correction',
        terminalResidual=dict(format='RAH247RS',columns='globalID Rmass Rmx Rmy Rmz REinternal RSA'.split(),
          stage='after saved complete cycle, flow preprocessing then native space; SA preprocessing then native space; no time advance',
          sideEffects='primitive/gradient/Jacobian/RHS caches refreshed; conserved/SA/mut identity checked; saved frame/CSV unchanged'),
        coarsePrimitiveStage='pre-update cache; decode current state from conserved variables',physicalQualification=False))
    cmd=read(PREV/'control_stage01/command.json')
    cmd=[str(NATIVE) if x==str(ROOT/'native/su2PhaseAwareAir/fasRocketControlV246.cpp') else x.replace(str(PREV/'build03'),str(BUILD)) for x in cmd]
    save(STAGE/'command.json',cmd);args=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-o','-MF'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    dep=STAGE/'inputs.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=str(WD),timeout=30)
    known=read(PREV/'control_stage01/dependencies.json');deps={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
        if p!=NATIVE:assert h==known[str(p)],str(p)
        deps[str(p)]=h
    save(STAGE/'dependencies.json',deps);shutil.copyfile(NATIVE,STAGE/NATIVE.name)
    save(STAGE/'report.json',dict(passed=True,updates=200,newUniformInitial=True,physicalQualification=False,
        points=335726,cells=705490,faces=19752,areaM2=1.9646292565715553,onlyConfigChanges=['ITER','OUTPUT_WRT_FREQ','SCREEN_WRT_FREQ_INNER'],numericalObjectsToCompile=0))


def build():
    subprocess.check_call(read(STAGE/'command.json'),cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'preflight.exe'),thinEntryOnly=True,numericalObjectsCompiled=0))


def setup_guard():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes
    guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    return guard


def main(stage):
    if CLOSED:raise SystemExit('v247 bounded development closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);g=setup_guard()
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',NATIVE,MESH,UNIFORM,FACES,
      PREV/'control_stage01/command.json',PREV/'control_stage01/dependencies.json',PREV/'control_audit02/report.json',
      PREV/'control_observe01/case.cfg',ROOT/'outputs/v245_validation/build01/members.json']
    if stage=='prepare':return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_prepare'],120,1536,locks,BASE)
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for group in read(ROOT/'outputs/v245_validation/build01/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    locks += [STAGE/'command.json',STAGE/'dependencies.json',CASE/'case.cfg',CASE/'schema.json',CASE/'mesh.su2',CASE/'original_faces.txt']
    if stage=='build':return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build'],120,1536,list(dict.fromkeys(locks)),BASE,'build')
    if stage=='run':
        assert read(BASE/'process_build/report.json')['passed'];locks += [BUILD/'preflight.exe']
        return g.supervise('rocket',[str(BUILD/'preflight.exe'),'case.cfg','.',str(UNIFORM),'original_faces.txt'],1800,1536,list(dict.fromkeys(locks)),CASE,'cfd')
    raise ValueError('TASK_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v247 bounded development closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
