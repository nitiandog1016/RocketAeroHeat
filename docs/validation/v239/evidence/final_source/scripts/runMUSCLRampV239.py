"""模块职责：复用成熟对象，登记渐进调度控制及唯一真实原域候选的完整资源合同。"""
from pathlib import Path
from types import SimpleNamespace
import difflib,json,os,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v239_validation';DOC=ROOT/'docs/validation/v239';OLD=ROOT/'outputs/v238_validation'
STAGE=BASE/'prepare01';BUILD=BASE/'build01';CASE=BASE/'rocket01';CONTROL=BASE/'control01'
CPP=ROOT/'native/su2PhaseAwareAir/rampBranchV239.cpp';GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'
MESH=ROOT/'outputs/v222_validation/rocket01/mesh.su2';OLDBUILD=ROOT/'tools/runtime/su2-thermal-build/build-v117';CLOSED=True
RAMP=dict(RAMP_MUSCL='YES',KIND_MUSCL_RAMP='SMOOTH_FUNCTION',RAMP_MUSCL_COEFF='(0, 1, 100)',RAMP_MUSCL_POWER='1')


def configure():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    guard.size_bytes=lambda p:sum(q.stat().st_size for q in BASE.rglob('*') if q.is_file()) if Path(p)==BASE or Path(p).name.startswith('process_') else 0


run=SimpleNamespace(configure=configure,guard=guard)


def config(path):return {k.strip():v.strip() for k,v in (l.split('=',1) for l in path.read_text(encoding='utf8').splitlines() if '=' in l)}
def write_config(path,c):path.write_text('\n'.join(k+'= '+v for k,v in c.items())+'\n',encoding='utf8')


def prepare():
    for p in (STAGE,BUILD,CONTROL):p.mkdir(exist_ok=False)
    oldInputs=read(OLD/'rocket_prepare01/inputs.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    common=OLD/'rocket_prepare01/commonBranchV238.hpp';source=ROOT/'native/su2PhaseAwareAir/uniformHLLCBranchV238.cpp'
    assert sha(common)==oldInputs[str(common)] and sha(source)==oldInputs[str(source)]
    before=common.read_text(encoding='utf8');assert before.count('if(completed!=200&&completed!=300)return;')==1
    after=before.replace('if(completed!=200&&completed!=300)return;','if(completed!=300&&completed!=400)return;')
    (STAGE/'commonBranchV239.hpp').write_text(after,encoding='utf8')
    (STAGE/'benefit_windows.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='commonBranchV238.hpp',tofile='commonBranchV239.hpp')),encoding='utf8')
    prefix=source.read_text(encoding='utf8').split('class UniformHLLCDriverV238 final:')[0].replace('"commonBranchV238.hpp"','"commonBranchV239.hpp"')
    (STAGE/'hllcCoreV239.hpp').write_text(prefix,encoding='utf8')
    cmd=read(OLD/'rocket_prepare01/commands.json')['helper'];cmd=list(cmd)
    cmd[cmd.index(str(source))]=str(CPP);cmd[cmd.index('-MF')+1]=str(BUILD/'driver.d');cmd[cmd.index('-o')+1]=str(BUILD/'driver.exe')
    cmd=[x.replace(str(OLD/'rocket_prepare01'),str(STAGE)) if x.startswith('-I') else x for x in cmd]
    save(STAGE/'command.json',cmd);filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(STAGE/'inputs.d')],cwd=str(OLDBUILD),timeout=20)
    bound={}
    for raw in (STAGE/'inputs.d').read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(raw);p=(p if p.is_absolute() else OLDBUILD/p).resolve();h=sha(p)
        if p!=CPP and p.parent!=STAGE:
            expected=oldInputs.get(str(p),tree.get(str(p)))
            if expected is None:expected=tree[str(ROOT/'tools/runtime/su2-thermal-build/source'/p.relative_to(ROOT/'outputs/v220_validation/source'))]
            assert h==expected,str(p)
        bound[str(p)]=h
    for name,members in read(OLD/'build04/members.json').items():
        bound[str(OLD/'build04'/name)]=sha(OLD/'build04'/name)
        for r in members:assert sha(r['path'])==r['sha256'];bound[r['path']]=r['sha256']
    bindings=read(OLD/'outcome02/report.json')['sourceBindings']
    for r in bindings:assert sha(r['path'])==r['sha256'];bound[r['path']]=r['sha256']
    # 只读控制与生产使用相同真实ramp日程；不以压缩时间的日程代替生产端点。
    c=config(OLD/'native01/case.cfg');c.update(RAMP);c.update(ITER='103',SCREEN_WRT_FREQ_INNER='103',OUTPUT_WRT_FREQ='103')
    write_config(CONTROL/'case.cfg',c)
    mesh=(CONTROL/c['MESH_FILENAME']).resolve();assert mesh.exists()
    for p in (CONTROL/'case.cfg',mesh,common,source,DOC/'PLAN.md',OLD/'build04/members.json'):bound[str(p)]=sha(p)
    save(STAGE/'inputs.json',bound);save(STAGE/'report.json',dict(passed=True,boundInputs=len(bound),numericalObjectsRecompiled=0,
      controlFlowUpdates=0,expectedFullOrderFirstUpdate=102,realCFDAdmitted=False))


def rocket_prepare():
    assert read(BASE/'control_audit01/report.json')['passed'];CASE.mkdir(exist_ok=False)
    old=config(OLD/'rocket01/case.cfg');c=dict(old);c.update(RAMP);write_config(CASE/'case.cfg',c)
    for name in ('uniform_six_state.txt','original_faces.txt'):shutil.copyfile(OLD/'rocket01'/name,CASE/name)
    contract=read(OLD/'rocket_contract.json');contract.update(branch='phase-air-hllc-native-muscl-ramp-v239',
      sourceObjectsRecompiledThisVersion=0,fixedNumericsPreservedExceptNamedHLLC=False,
      fixedNumericsPreservedExceptNativeMUSCLRamp=True,configChanges={k:dict(old=old.get(k),new=v) for k,v in c.items() if old.get(k)!=v},
      expectedFirstFullOrderUpdate=102,checkUpdates=[300,400],adjustmentUpdates=101,rampSchedule=RAMP,uniformFreestreamInitial=True,
      fixedNumericsPreserved=False,fixedFullOrderNumericsPreserved=True,
      initialStateSource='v238 same uniform freestream',sourceCheckpoint=None,
      realInitialUniformSixMustMatchV238=True,controlAuditSHA256=sha(BASE/'control_audit01/report.json'))
    # 本轮唯一变化是原生空间阶次日程及同算子窗口位置，旧数值与物性身份保持。
    assert set(contract['configChanges'])==set(RAMP)
    assert sha(MESH)==contract['meshSHA256'] and sha(GEOMETRY)==contract['geometrySHA256']
    save(BASE/'rocket_contract.json',contract)
    bound=read(STAGE/'inputs.json')
    for p,h in bound.items():assert sha(p)==h,p
    for p in (MESH,GEOMETRY,CASE/'case.cfg',CASE/'uniform_six_state.txt',CASE/'original_faces.txt',DOC/'ROCKET_PLAN.md',
              BASE/'rocket_contract.json',BASE/'control_audit01/report.json',BUILD/'driver.exe',Path(__file__).resolve()):bound[str(p)]=sha(p)
    save(BASE/'rocket_inputs.json',bound);save(BASE/'rocket_admission.json',dict(passed=True,unique=True,
      contractSHA256=sha(BASE/'rocket_contract.json'),inputsSHA256=sha(BASE/'rocket_inputs.json'),physicalQualification=False))


def main(stage):
    if CLOSED:raise SystemExit('v239 execution closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',CPP]
    if stage=='prepare':return guard.supervise('prepare',[sys.executable,str(Path(__file__).resolve()),'task_prepare'],60,1024,locks,BASE)
    for p,h in read(STAGE/'inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
    locks += [STAGE/'inputs.json',STAGE/'command.json']
    if stage=='build':return guard.supervise('build',read(STAGE/'command.json'),60,1024,locks,OLDBUILD,'build')
    if stage=='control':
        assert read(BASE/'process_build/report.json')['passed'];locks.append(BUILD/'driver.exe')
        return guard.supervise('control',[str(BUILD/'driver.exe'),'case.cfg','.','--schedule-control'],30,1024,locks,CONTROL,'query')
    if stage=='admit':
        locks += [DOC/'ROCKET_PLAN.md',BASE/'control_audit01/report.json',OLD/'rocket_contract.json',OLD/'rocket01/case.cfg',OLD/'rocket01/uniform_six_state.txt',OLD/'rocket01/original_faces.txt']
        return guard.supervise('admit',[sys.executable,str(Path(__file__).resolve()),'task_admit'],30,1024,locks,BASE)
    if stage=='run':
        assert read(BASE/'rocket_admission.json')['passed'] and not (BASE/'process_rocket').exists()
        for p,h in read(BASE/'rocket_inputs.json').items():assert sha(p)==h,p;locks.append(Path(p))
        locks += [BASE/'rocket_admission.json',BASE/'rocket_inputs.json']
        return guard.supervise('rocket',[str(BUILD/'driver.exe'),'case.cfg','.','uniform_six_state.txt','original_faces.txt'],2400,1024,list(dict.fromkeys(locks)),CASE,'cfd',True)
    raise ValueError('STAGE_NOT_ADMITTED')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v239 execution closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_admit':rocket_prepare()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
