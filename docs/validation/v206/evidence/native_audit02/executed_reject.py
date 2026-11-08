"""原生构建后的有限恢复：只读已有正向结果，继续未执行的拒绝检查。"""
from pathlib import Path
import json,shutil,time,sys
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from phaseAwareDryAir import PhaseAwareDryAir,MODEL,DOMAIN,SOURCE,pressure_ceiling,dew_pressure

BASE=ROOT/'outputs/v206_validation'
OLD=BASE/'native01'
OUT=BASE/'native_audit02'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def local(p):
    return ROOT/p.split('RocketAeroHeat/',1)[1] if 'RocketAeroHeat/' in p else Path(p)

def prepare():
    """模块一：Windows已有scipy完成只读独立审核，绝不再次运行原生正向。"""
    start=time.monotonic();OUT.mkdir(exist_ok=False);shutil.copyfile(__file__,OUT/'executed_prepare.py')
    progress=read(OLD/'progress.json');gas=PhaseAwareDryAir()
    for key in ['library','helper','positive']:
        receipt=progress[key]
        if not receipt['passed'] or receipt['liveGroupAfterCleanup']:raise ValueError('既有阶段未通过')
        for p,h in receipt['inputSHA256'].items():
            if p.startswith('/opt/') or p.startswith('/usr/'):continue
            if sha(local(p))!=h:raise ValueError('既有准确输入改变 '+p)
    lines=(OLD/'positive/stdout.txt').read_text().splitlines()
    rows=np.array([[float(v) for v in line.split()[1:]] for line in lines if line.startswith('ROW ')])
    cases=[r.split() for r in (OLD/'positive.txt').read_text().splitlines()]
    target=np.array([float(r[2]) if r[0]=='T' else float(gas.inverse(float(r[4]),float(r[2]),'e')) for r in cases])
    if rows.shape!=(10,13):raise ValueError('正向行列数改变')
    terr=float(abs(rows[:,3]-target).max());eerr=float(np.max(abs(rows[:,4]-rows[:,1])/np.maximum(1,abs(rows[:,1]))))
    values=gas.evaluate(rows[:,5],rows[:,3]);keys=['cpJkgK','cvJkgK','enthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+6]-values[k])/np.maximum(1,abs(values[k])))) for i,k in enumerate(keys)}
    perror=float(np.max(abs(rows[:,12]/pressure_ceiling(rows[:,3])-1)))
    if terr>2e-6 or eerr>2e-10 or max(errors.values())>2e-9 or perror>2e-14:raise ValueError('正向误差未过')
    elo=float(gas.evaluate(1,60)['internalEnergyJkg']);ehi=float(gas.evaluate(1,2000)['internalEnergyJkg'])
    cases_bad=[('e-low','E 10 %.17g 100 10'%(elo-1),1,'RAH_PHASE_ENERGY_OUTSIDE_TABLE'),
        ('e-high','E 10 %.17g 100 10'%(ehi+1),1,'RAH_PHASE_ENERGY_OUTSIDE_TABLE'),
        ('e-nan','E 10 nan 100 10',1,'RAH_PHASE_INVERSE_INPUT'),('old-p','E -1 0 300 1',1,'RAH_PHASE_INVERSE_INPUT'),
        ('seed-low','T 10 80 59 10',1,'RAH_PHASE_INVERSE_INPUT'),('final-p-zero','T 10 80 100 0',10,'RAH_PHASE_PROBE_FINAL_REJECT'),
        ('final-p-high','T 10 1200 100 1000001',10,'RAH_PHASE_PROBE_FINAL_REJECT'),
        ('final-dew','T 10 60 300 %.17g'%(float(dew_pressure(60))*.5),10,'RAH_PHASE_PROBE_FINAL_REJECT')]
    for label,line,_,_ in cases_bad:(OUT/(label+'.txt')).write_text(line+'\n')
    report=dict(passed=True,model=MODEL,domain=DOMAIN,positiveStates=10,temperatureMaxAbsErrorK=terr,
        energyMaxScaledError=eerr,propertyMaxScaledErrors=errors,pressureCeilingMaxRelative=perror,
        v205InputConditionalRootK=float(rows[-1,3]),conditionalRootIsNotSavedPhysicalState=True,
        registrationLibrarySHA256=sha(OLD/'library/librahPhaseAwareAir.so'),
        queryBinarySHA256=sha(OLD/'probe/rahPhaseAwareAirProbe'),tableSHA256=sha(OLD/'mixture.dict'),
        positiveLogSHA256=sha(OLD/'positive/stdout.txt'),cases=cases_bad,elapsedSeconds=time.monotonic()-start,
        originalPlanPassed=False,originalFailure='post-positive Python ModuleNotFoundError: scipy',
        referenceResourcePlanPassed=False,referencePeakRSSBytes=None)
    save(OUT/'positive_audit.json',report);print(json.dumps(report))

def reject():
    """模块二：仅继续未运行的8拒绝，沿用既有二进制及准确源码。"""
    from runOpenFoamProcess import run
    start=time.monotonic();report=read(OUT/'positive_audit.json');progress=read(OLD/'progress.json')
    marker=OUT/'reject_admission.json'
    if marker.exists():raise ValueError('拒绝阶段已启动，禁止重复')
    binary=OLD/'probe/rahPhaseAwareAirProbe';libpath=OLD/'library/librahPhaseAwareAir.so'
    for p,h in progress['positive']['inputSHA256'].items():
        if sha(Path(p))!=h:raise ValueError('准确构建/查询输入改变 '+p)
    if sha(libpath)!=report['registrationLibrarySHA256'] or sha(binary)!=report['queryBinarySHA256']:raise ValueError('库改变')
    shutil.copyfile(__file__,OUT/'executed_reject.py')
    locks=[Path(p) for p in progress['positive']['inputSHA256']]+[Path(__file__).resolve(),OUT/'positive_audit.json',ROOT/'docs/validation/v206/RECOVERY_PLAN.md']
    save(marker,dict(totalSeconds=60,parentSeconds=75,memoryMiB=2048,outputMiB=96,automaticRetry=0,inputSHA256={str(p):sha(p) for p in locks}))
    rejects=[]
    for label,_,code,token in report.pop('cases'):
        sample=OUT/(label+'.txt');left=60-(time.monotonic()-start)
        if left<=0:raise TimeoutError('总查询预算耗尽')
        proc=run([binary,OLD/'mixture.dict',sample],OUT,OUT/label,min(8,left),locks+[sample,marker],96,2048)
        text=(OUT/label/'stdout.txt').read_text()
        ok=proc['nativeExitCode']==code and token in text and not proc['liveGroupAfterCleanup'] and all(proc['inputStable'].values())
        rejects.append(dict(case=label,passed=bool(ok),expectedExitCode=code,process=proc));save(OUT/'reject_progress.json',rejects)
        if not ok:raise ValueError('拒绝不符 '+label)
    report.update(rejections=rejects,actualSolverDispatchVerified=False,actualPerStepGuardVerified=False,
        oldPressureSeparatedFromFinalPressure=True,CFDRuns=0,physicalQualification=False,
        libraryBuild=progress['library'],probeBuild=progress['helper'],positive=progress['positive'],
        rejectionElapsedSeconds=time.monotonic()-start)
    save(OUT/'report.json',report);print(json.dumps(dict(passed=True,rejections=len(rejects),seconds=report['rejectionElapsedSeconds'])))

if __name__=='__main__':
    if sys.argv[1:]==['prepare']:prepare()
    elif sys.argv[1:]==['reject']:reject()
    else:raise ValueError('用prepare只读审核或reject有限继续；不重建。')
