"""模块职责：在真实几何预检通过后，监督唯一固定五步观察开/关完整FAS控制。"""
from pathlib import Path
import os,subprocess,sys
import runFASRocketPreflightV246 as prep
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=prep.BASE;DOC=prep.DOC;STAGE=BASE/'control_stage01';BUILD=BASE/'build03';CLOSED=False
NATIVE=ROOT/'native/su2PhaseAwareAir/fasRocketControlV246.cpp'


def prepare():
    assert read(BASE/'audit01/report.json')['geometryPreflightPassed'];STAGE.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    cmd=read(prep.STAGE/'command.json');cmd=[str(NATIVE) if x==str(prep.NATIVE) else x.replace(str(prep.BUILD),str(BUILD)) for x in cmd]
    save(STAGE/'command.json',cmd);args=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-o','-MF'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        args.append(x);i+=1
    d=STAGE/'inputs.d';subprocess.check_call(args+['-M','-MT','dep','-MF',str(d)],cwd=str(prep.WD),timeout=30)
    known=read(prep.STAGE/'dependencies.json');known.update(read(prep.PREV/'stage01/dependencies.json'))
    old=read(ROOT/'outputs/v235_validation/stage01/dependencies.json');known.update(old)
    deps={}
    for name in d.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else prep.WD/p).resolve();h=sha(p)
        if p!=NATIVE:assert h==known[str(p)],str(p)
        deps[str(p)]=h
    save(STAGE/'dependencies.json',deps)
    cfg={}
    for line in (prep.CASE/'case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(WRT_RESTART_OVERWRITE='YES',SCREEN_WRT_FREQ_INNER='1')
    for mode in ('observe','off'):
        case=BASE/('control_'+mode+'01');case.mkdir(exist_ok=False)
        os.link(str(prep.MESH),str(case/'mesh.su2'))
        (case/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    save(STAGE/'report.json',dict(passed=True,updatesEach=5,sameFixedNumericalCombination=True,numericalObjectsToCompile=0,physicalQualification=False))


def build():
    subprocess.check_call(read(STAGE/'command.json'),cwd=str(prep.WD))
    save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'preflight.exe'),thinEntryOnly=True,numericalObjectsCompiled=0))


def main(stage):
    if CLOSED:raise SystemExit('v246 fixed real control closed; follow NEXT_DECISION')
    g=prep.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(prep.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',NATIVE,
      ROOT/'native/su2PhaseAwareAir/singleDomainRunV235.cpp',ROOT/'native/su2PhaseAwareAir/FASHooksV245.hpp',
      BASE/'audit01/report.json',prep.CASE/'case.cfg',prep.STAGE/'command.json',prep.MESH]
    if stage=='prepare':return g.supervise('control_prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,1536,locks,BASE)
    for n,h in read(STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for group in read(prep.PREV/'build01/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    locks += [STAGE/'command.json',STAGE/'dependencies.json',prep.PREV/'build01/members.json']
    if stage=='build':return g.supervise('control_build',[sys.executable,str(Path(__file__).resolve()),'task_build'],60,1536,list(dict.fromkeys(locks)),BASE,'build')
    if stage in ('observe','off'):
        assert read(BASE/'process_control_build/report.json')['passed'];case=BASE/('control_'+stage+'01');locks += [case/'case.cfg',case/'mesh.su2',BUILD/'preflight.exe']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_control_*/report.json') if read(p)['budgetKind']=='cfd')
        return g.supervise('control_'+stage,[str(BUILD/'preflight.exe'),'case.cfg','.',stage],180-used,1536,list(dict.fromkeys(locks)),case,'cfd')
    raise ValueError('TASK_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v246 fixed real control closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
