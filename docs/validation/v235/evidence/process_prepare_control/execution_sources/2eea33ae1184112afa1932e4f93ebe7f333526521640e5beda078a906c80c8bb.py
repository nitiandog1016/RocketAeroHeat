"""模块职责：有限构建与完整原生flow-SA观察中性控制；不准入真实箭体运行。"""
from pathlib import Path
import json,shutil,subprocess,sys,time
import prepareSingleDomainV235 as prep
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=prep.BASE;DOC=prep.DOC;BUILD=prep.BUILD;CLOSED=False


def build():
    cmds=read(prep.STAGE/'commands.json');receipts=[]
    for name,cmd in cmds.items():
        if name=='helper':continue
        t=time.monotonic();subprocess.check_call(cmd,cwd=str(prep.BUILDOLD));receipts.append(dict(unit=name,seconds=time.monotonic()-t));save(BUILD/'compile_receipts.json',receipts)
    ar=str(Path(cmds['helper'][0]).parent/'ar.exe');members=read(prep.GAS/'build02/members.json');changed=0;new={}
    for name,rows in members.items():
        result=[]
        for r in rows:
            p=Path(r['path']);assert sha(p)==r['sha256']
            if p.name=='CNSVariable.obj':p=BUILD/'CNSVariable.obj';changed+=1
            elif p.name=='integration_CIntegration.cpp.obj':p=BUILD/'CIntegration.obj';changed+=1
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/name)]+[r['path'] for r in result],cwd=str(prep.BUILDOLD));new[name]=result
    assert changed==2;save(BUILD/'members.json',new)
    subprocess.check_call(cmds['helper'],cwd=str(prep.BUILDOLD))
    save(BUILD/'report.json',dict(passed=True,helperSHA256=sha(BUILD/'control.exe'),changedObjects=2,
      fullEngineBuilt=False,numericalFormulaChanges=0,flowUpdates=0))


def prepare_control():
    """同一个既有小壳，真实5次flow-SA更新；只比较观察开/关，不安装制造场或扫描参数。"""
    original=ROOT/'outputs/v221_validation/control01';text=(original/'case.cfg').read_text(encoding='utf8')
    cfg={k.strip():v.strip() for k,v in (s.split('=',1) for s in text.splitlines() if '=' in s)}
    cfg.update(ITER='5',OUTPUT_WRT_FREQ='5',SCREEN_WRT_FREQ_INNER='1',WRT_RESTART_OVERWRITE='YES',CONV_STARTITER='10')
    for name in ('off01','observe01'):
        case=BASE/name;case.mkdir(exist_ok=False);shutil.copyfile(original/'mesh.su2',case/'mesh.su2')
        (case/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    save(BASE/'control_contract.json',dict(sameMeshSHA256=sha(BASE/'off01/mesh.su2'),sameConfigSHA256=sha(BASE/'off01/case.cfg'),
      updatesEach=5,existingPoints=24,changedOnlyObservationCallback=True,fullRocketAdmitted=False))


def main(stage):
    if CLOSED:raise SystemExit('v235 native entry closed; follow NEXT_DECISION')
    g=prep.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/prepareSingleDomainV235.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',prep.STAGE/'commands.json']
    for n,h in read(prep.STAGE/'dependencies.json').items():assert sha(n)==h,n;locks.append(Path(n))
    for group in read(prep.GAS/'build02/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
    if stage=='build':return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_build'],120,128,list(dict.fromkeys(locks)),BASE,'build')
    if stage=='prepare_control':
        locks += [ROOT/'outputs/v221_validation/control01/case.cfg',ROOT/'outputs/v221_validation/control01/mesh.su2']
        return g.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_control'],120,128,list(dict.fromkeys(locks)),BASE)
    if stage in ('off','observe'):
        assert read(BASE/'process_build/report.json')['passed'] and read(BUILD/'report.json')['passed']
        case=BASE/(stage+'01');locks += [BUILD/'control.exe',case/'case.cfg',case/'mesh.su2']
        used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='query')
        return g.supervise(stage,[str(BUILD/'control.exe'),'case.cfg','.',stage,'5'],60-used,128,list(dict.fromkeys(locks)),case,'query')
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v235 native entry closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_build':build()
    elif sys.argv[1]=='task_control':prepare_control()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
