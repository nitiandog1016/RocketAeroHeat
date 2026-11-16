"""模块职责：绑定已有原生NK对象，准备并监督唯一一次真实状态完整控制。"""
from pathlib import Path
import os,shutil,subprocess,sys
import prepareOuterStencilV249 as task
from terminalCapacityV249 import supervise as terminal_supervise
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=task.BASE;DOC=task.DOC;CASE=BASE/'native01';BUILD=BASE/'build01';STAGE=BASE/'stage01';CLOSED=True
PREV=ROOT/'outputs/v248_validation';WD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
NATIVE=ROOT/'native/su2PhaseAwareAir/newtonControlV249.cpp';REFERENCE=PREV/'rocket01/state_000200.bin'


def prepare():
    decision=read(BASE/'method01/report.json');assert decision['passed'] and not decision['productionAdopted']
    for p in (CASE,BUILD,STAGE):p.mkdir(exist_ok=False)
    cfg={}
    for line in (PREV/'rocket01/case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    for name in list(cfg):
        if name.startswith('MG'):del cfg[name]
    changes=dict(NEWTON_KRYLOV='YES',NEWTON_KRYLOV_IPARAM='(0, 3, 2, 0)',NEWTON_KRYLOV_DPARAM='(0.0, 0.1, -3.0, 1e-4, 1.0)',
      MGLEVEL='0',RESTART_SOL='YES',SOLUTION_FILENAME='input_restart',ITER='1',OUTPUT_WRT_FREQ='1',SCREEN_WRT_FREQ_INNER='1')
    cfg.update(changes);(CASE/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    for name,src in [('mesh.su2',PREV/'rocket01/mesh.su2'),('input_restart.csv',PREV/'rocket01/restart_flow.csv')]:os.link(str(src),str(CASE/name))
    save(CASE/'schema.json',dict(format='RAH235LE',columns='globalID x y z volume rho rhou rhov rhow rhoEinternal SA T u v w p rhoPrimitive Hinternal a mu mutFlow k Cp mutSA wallDistance flowRelaxation flowCFL flowWorking0 flowWorking1 flowWorking2 flowWorking3 flowWorking4 saWorking reserved'.split(),
      initialWorkUnknown=True,workingArrayIsFinalResidual=False,initialStateSource='v248200',nativeNewtonStartupDisabled=True,physicalQualification=False))
    cmd=read(PREV/'stage02/command.json');cmd=[str(NATIVE) if x==str(ROOT/'native/su2PhaseAwareAir/fasRecoveryDevelopmentV248.cpp') else x.replace(str(PREV/'build01'),str(BUILD)) for x in cmd]
    save(STAGE/'command.json',cmd);args=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-o','-MF'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    dep=STAGE/'inputs.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(dep)],cwd=str(WD),timeout=30)
    known=read(PREV/'stage02/dependencies.json');deps={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();h=sha(p)
        if str(p) in known:assert h==known[str(p)],str(p)
        deps[str(p)]=h
    save(STAGE/'dependencies.json',deps);shutil.copyfile(NATIVE,STAGE/NATIVE.name)
    save(STAGE/'report.json',dict(passed=True,method='native Newton-Krylov',changes=changes,updates=1,
      referenceSHA256=sha(REFERENCE),sixStateMustMatchBeforeRun=True,numericalObjectsToCompile=0))


def build():
    subprocess.check_call(read(STAGE/'command.json'),cwd=str(WD))
    save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'preflight.exe'),thinOnly=True,numericalObjectsCompiled=0))


def main(mode):
    if CLOSED:raise SystemExit('v249 native control closed; follow NEXT_DECISION')
    g=task.setup_guard();locks=[Path(__file__).resolve(),Path(task.__file__),ROOT/'scripts/terminalCapacityV249.py',NATIVE,DOC/'PLAN.md',DOC/'NATIVE_CONTROL_PLAN.md',
      BASE/'method01/report.json',REFERENCE,PREV/'rocket01/case.cfg',PREV/'rocket01/mesh.su2',PREV/'rocket01/restart_flow.csv',
      PREV/'stage02/command.json',PREV/'stage02/dependencies.json',ROOT/'outputs/v245_validation/build01/members.json']
    if mode=='prepare':return terminal_supervise(g,BASE,'prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,512,locks,BASE,'offline')
    for name,h in read(STAGE/'dependencies.json').items():assert sha(name)==h,name;locks.append(Path(name))
    for members in read(ROOT/'outputs/v245_validation/build01/members.json').values():
        for item in members:assert sha(item['path'])==item['sha256'];locks.append(Path(item['path']))
    locks += [STAGE/'command.json',STAGE/'dependencies.json',CASE/'case.cfg',CASE/'mesh.su2',CASE/'input_restart.csv',CASE/'schema.json']
    if mode=='build':return terminal_supervise(g,BASE,'build',[sys.executable,str(Path(__file__).resolve()),'task_build'],120,512,list(dict.fromkeys(locks)),BASE,'build')
    if mode=='run':
        assert read(BASE/'process_build/report.json')['passed'];locks.append(BUILD/'preflight.exe')
        return terminal_supervise(g,BASE,'native',[str(BUILD/'preflight.exe'),'case.cfg','.',str(REFERENCE)],180,512,list(dict.fromkeys(locks)),CASE,'cfd')
    raise ValueError('MODE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v249 native control closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
