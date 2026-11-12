"""模块职责：准确源绑定后，仅将物性域拒绝接到原生变量回退；隔离编译和控制。"""
from pathlib import Path
import csv,difflib,json,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from capturePhaseEulerFailureV231 import ROOT,BASE,DOC,PREV,read,save,sha
from prepareSU2PhaseAirV220 import OLD_SOURCE,OLD_BUILD
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy

BUILD=BASE/'recovery_build02';NATIVE=ROOT/'native/su2PhaseAwareAir'
PREP=BASE/'recovery_prepare04'
CONTROL=NATIVE/'primitiveRecoveryControl.cpp';HEADER=NATIVE/'PhaseAirPrimitiveRecovery.hpp'
CLOSED=False

def prepare():
    out=PREP;out.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=True)
    assert read(BASE/'audit01/report.json')['originalFirst30TrajectoryByteIdentical']
    old_tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    old=OLD_SOURCE/'SU2_CFD/src/variables/CEulerVariable.cpp';assert sha(old)==old_tree[str(old)]
    text=old.read_text(encoding='utf8');original=text
    text=text.replace('#include "../../include/variables/CEulerVariable.hpp"','#include "variables/CEulerVariable.hpp"')
    text=text.replace('#include "../../include/fluid/CFluidModel.hpp"','#include "fluid/CFluidModel.hpp"\n#include "PhaseAirPrimitiveRecovery.hpp"')
    needle='  FluidModel->SetTDState_rhoe(density, staticEnergy);';assert text.count(needle)==2
    text=text.replace(needle,'  // 具名物性域拒绝进入下方原生Solution_Old回退，合法状态不改值。\n  const bool rah_valid_state = rah_phase_recovery::set_trial_state(FluidModel, density, staticEnergy);',1)
    for lhs in ('bool check_dens  =','bool check_press =','bool check_sos   =','bool check_temp  ='):
        assert text.count(lhs)==1;text=text.replace(lhs,lhs+' !rah_valid_state ||')
    staged=out/'CEulerVariable.cpp';staged.write_text(text,encoding='utf8')
    (out/'primitive_recovery.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),
        fromfile='SU2-8.6/CEulerVariable.cpp',tofile='phase-euler-native-primitive-rejection-v1/CEulerVariable.cpp')),encoding='utf8')
    shutil.copyfile(old,out/'CEulerVariable.original.cpp');shutil.copyfile(HEADER,out/HEADER.name);shutil.copyfile(CONTROL,out/CONTROL.name)
    sources={str(old):sha(old)}
    for rel in ('SU2_CFD/src/integration/CMultiGridIntegration.cpp',):
        p=OLD_SOURCE/rel;assert sha(p)==old_tree[str(p)];sources[str(p)]=sha(p)
    p=ROOT/'outputs/v229_validation/source/SU2_CFD/src/solvers/CEulerSolver.cpp'
    expected=read(ROOT/'outputs/v229_validation/process_build/report.json')['inputSHA256'][str(p)]
    assert sha(p)==expected;sources[str(p)]=sha(p)
    save(out/'sources.json',sources)
    common=read(PREV/'build02/command.json');i=common.index('-MD');flags=common[:i]
    flags=[x.replace('-O2','-O3') for x in flags]+['-D_FILE_OFFSET_BITS=64','-I'+str(NATIVE)]
    archives=common[common.index('-Wl,--start-group'):common.index('-o')]
    compile_cmd=flags+['-MD','-MF',str(BUILD/'primitive.d'),'-c',str(staged),'-o',str(BUILD/'primitive.obj')]
    commands={'object':compile_cmd}
    for role in ('original','recovered'):
        commands[role]=(flags+['-MD','-MF',str(BUILD/(role+'.d')),str(CONTROL)]+
            ([str(BUILD/'primitive.obj')] if role=='recovered' else [])+archives+['-o',str(BUILD/(role+'.exe'))])
    commands['driver']=(flags+['-MD','-MF',str(BUILD/'rocket.d'),str(NATIVE/'rocketEulerRun.cpp'),str(BUILD/'primitive.obj')]+
        archives+['-o',str(BUILD/'rocket.exe')])
    save(BUILD/'commands.json',commands)
    for name,source in [('object',staged),('control',CONTROL)]:
        dep=out/(name+'.d');subprocess.check_call(flags+['-M','-MT','dep','-MF',str(dep),str(source)],cwd=str(OLD_BUILD),timeout=15)
        bindings={};previous=read(PREV/'build02/dependencies.json')
        for n in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
            if p==out/HEADER.name:assert h==sha(HEADER)
            elif p not in (staged,HEADER,CONTROL):assert previous[str(p)]==h,str(p)
            bindings[str(p)]=h
        save(out/(name+'_dependencies.json'),bindings)
    raw=np.loadtxt(str(BASE/'capture01/termination_raw_state.csv'),delimiter=',',skiprows=1)
    saved=np.loadtxt(str(PREV/'rocket01/state_000001.csv'),delimiter=',',skiprows=1)
    trial=raw[raw[:,0]==17337,1:6][0];prior=saved[saved[:,0]==17337,4:9][0]
    np.savetxt(str(out/'point17337.csv'),np.array([prior,trial]),delimiter=',',header='rho,rhou,rhov,rhow,rhoE',comments='',fmt='%.17g')
    finish()

def finish():
    # 仅完成已生成准备件的配置核对；既有控制目录不覆盖、不重生成源或依赖。
    out=PREP;staged=out/'CEulerVariable.cpp';old=out/'CEulerVariable.original.cpp'
    for name in ('object','control'):
        for n,h in read(out/(name+'_dependencies.json')).items():assert sha(n)==h,n
    for role in ('original','recovered'):
        dst=BASE/('control_'+role);dst.mkdir(exist_ok=True)
        assert not (dst/'control.json').exists()
        for n in ('case.cfg','mesh.su2'):
            source=ROOT/'outputs/v229_validation/control02'/n
            if (dst/n).exists():assert sha(dst/n)==sha(source)
            else:shutil.copyfile(source,dst/n)
    save(out/'report.json',dict(passed=True,identity='phase-euler-native-primitive-rejection-v1',
      sourceSHA256=sha(staged),originalSourceSHA256=sha(old),helperSHA256=sha(CONTROL),headerSHA256=sha(HEADER),
      actualBadPointID=17337,oldReference='saved update1 of same global node; not original step30 Solution_Old',
      phaseDomainUnchanged=True,sourcePatchOnly=True,buildExecuted=False))

def audit():
    out=BASE/'recovery_audit01';out.mkdir(exist_ok=False)
    def rows(role):
        with (BASE/('control_'+role)/'primitive.csv').open() as f:return {r['case']:r for r in csv.DictReader(f)}
    orig=rows('original');fixed=rows('recovered');assert set(orig)==set(fixed) and len(orig)==8
    data=np.loadtxt(str(PREP/'point17337.csv'),delimiter=',',skiprows=1);old,bad=data
    assert orig['valid']==fixed['valid'] and orig['cold_negative_physical_e']==fixed['cold_negative_physical_e']
    def values(row):return np.array([float(row['U%d'%i]) for i in range(5)])
    for name in ('actual_raw_invalid','high_energy','negative_density','nonfinite','joint_domain'):
        assert orig[name]['returned']=='-1' and orig[name]['exception'].startswith('RAH_PHASE_')
        assert fixed[name]['returned']=='0' and fixed[name]['exception']=='none'
        assert np.array_equal(values(fixed[name]),old),name
    for role in (orig,fixed):
        assert role['invalid_old']['returned']=='-1' and role['invalid_old']['exception']=='RAH_PHASE_ENERGY_ENTROPY_DOMAIN'
        assert np.array_equal(values(role['invalid_old']),bad)
    air=PhaseAwareViscousAir();errors=[]
    for name,row in fixed.items():
        if row['exception']!='none':continue
        U=values(row);rho=U[0];v=U[1:4]/rho;e=U[4]/rho-.5*(v@v)-300000
        T=float(inverse_energy(air,np.array([e]))[0]);p=rho*air.R*T;q=air.evaluate(p,T)
        expected=np.array([p,T,air.R*T,rho*air.R/(float(q['cpJkgK'])-air.R)])
        actual=np.array([float(row[k]) for k in ('p','T','dpdrho','dpde')])
        errors.append(float(np.max(abs(actual-expected)/np.maximum(1,abs(expected)))))
        if name=='cold_negative_physical_e':assert e<0 and row['returned']=='1'
    assert max(errors)<2e-9,max(errors)
    for role in ('original','recovered'):assert read(BASE/('control_'+role)/'control.json')['flowUpdates']==0
    save(out/'report.json',dict(passed=True,validTwoCasesByteIdentical=True,rejectedTrialsRecoveredByNativeOldState=5,
      invalidOldStillRejected=True,otherGasPropagated=True,numericalErrorsNotRecovered=True,
      physicalNegativeInternalEnergyPreserved=True,maxIndependentPrimitiveRelativeError=max(errors),
      newFlowUpdates=0,fullRocketFixedNotExecuted=True,scientificAcceptance=False,newQualifiedAreaM2=0))
    print(json.dumps(read(out/'report.json')))

def driver_parity():
    # 变量对象沿用上游O3；观测薄驱动恢复v230的O2，避免混入观测累加编译差异。
    target=BASE/'recovery_build03';target.mkdir(exist_ok=False)
    cmd=read(BUILD/'commands.json')['driver']
    cmd=[x.replace('-O3','-O2').replace(str(BUILD/'rocket.d'),str(target/'rocket.d'))
          .replace(str(BUILD/'rocket.exe'),str(target/'rocket.exe')) for x in cmd if x!='-D_FILE_OFFSET_BITS=64']
    save(target/'command.json',cmd)
    save(target/'preparation.json',dict(passed=True,driverOptimization='O2',primitiveObjectOptimization='O3',
      originalUnexecutedO3DriverSHA256=sha(BUILD/'rocket.exe'),primitiveObjectSHA256=sha(BUILD/'primitive.obj'),
      flowExecuted=False))

def main(stage):
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/capturePhaseEulerFailureV231.py',ROOT/'scripts/runSU2RocketPhaseV222.py',
      DOC/'PLAN.md',HEADER,CONTROL,PREV/'build02/command.json',PREV/'build02/dependencies.json']
    locks += [Path(n) for n in read(PREV/'build02/dependencies.json')]
    kind='offline';seconds=120;cwd=BASE
    if stage in ('prepare','prepare2','prepare3','prepare4'):cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare']
    elif stage=='finish':
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_finish']
        locks += [p for p in PREP.iterdir() if p.is_file()]
    elif stage=='driver_parity':
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_driver_parity']
        locks += [BUILD/'commands.json',BUILD/'primitive.obj',BUILD/'rocket.exe']
    elif stage=='driver2':
        cmd=read(BASE/'recovery_build03/command.json');kind='build';seconds=30
        locks += [BASE/'recovery_build03/command.json',BUILD/'primitive.obj']
    elif stage in ('object','object2','original','recovered','driver'):
        assert read(PREP/'report.json')['passed']
        cmd=read(BUILD/'commands.json')['object' if stage=='object2' else stage];kind='build';seconds=30 if stage=='object2' else 60
        if stage in ('recovered','driver'):assert read(BASE/'process_object2/report.json')['passed'];locks += [BUILD/'primitive.obj']
        locks += [PREP/'CEulerVariable.cpp',BUILD/'commands.json']
        for n in ('object','control'):locks += [Path(p) for p in read(PREP/(n+'_dependencies.json'))]
    elif stage.startswith('control_'):
        role=stage[8:];assert role in ('original','recovered') and read(BASE/('process_'+role)/'report.json')['passed']
        cwd=BASE/stage;cmd=[str(BUILD/(role+'.exe')),'case.cfg','.',str(PREP/'point17337.csv')]
        locks += [BUILD/(role+'.exe'),cwd/'case.cfg',cwd/'mesh.su2',PREP/'point17337.csv']
        kind='query';seconds=15
    elif stage=='audit':
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_audit']
        locks += [p for d in ('control_original','control_recovered') for p in (BASE/d).iterdir() if p.is_file()]
        locks += [ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/phaseAwareDryAir.py']
    else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise('recovery_'+stage if stage in ('prepare','prepare2','prepare3','prepare4','audit') else stage,
      cmd,seconds,256,list(dict.fromkeys(locks)),cwd,kind,False)

if __name__=='__main__':
    if CLOSED:raise SystemExit('v231 primitive recovery controls closed; no automatic new run')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_finish':finish()
    elif sys.argv[1]=='task_driver_parity':driver_parity()
    elif sys.argv[1]=='task_audit':audit()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
