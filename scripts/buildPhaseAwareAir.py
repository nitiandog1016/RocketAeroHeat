"""相态有界空气的隔离原生注册/每步守卫构建及限定查询。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from phaseAwareDryAir import PhaseAwareDryAir,MODEL,DOMAIN,SOURCE,pressure_ceiling,dew_pressure

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    from runOpenFoamProcess import run
    start=time.monotonic();base=ROOT/'outputs/v206_validation';out=base/'native01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');gas=PhaseAwareDryAir()
    for name in ['reference01','python_tests01']:
        if not read(base/name/'report.json')['passed']:raise ValueError('离线前置未过')
    library=out/'library';probe=out/'probe'
    shutil.copytree(ROOT/'native/openfoamPhaseAwareAir',library);shutil.copytree(ROOT/'native/openfoamPhaseAwareAirProbe',probe)
    (out/'mixture.dict').write_text(gas.dictionary())
    # 模块一：绑定实际原生调度、热物性源码及成熟求根头，不改变其内容。
    native_names=['src/OpenFOAM/db/Time/Time.C','src/OpenFOAM/db/functionObjects/functionObject/functionObject.H',
        'applications/solvers/foamRun/foamRun.C','applications/modules/shockFluid/thermophysicalPredictor.C',
        'applications/modules/shockFluid/pressureCorrector.C']
    original=[Path('/opt/openfoam13')/n for n in native_names]
    for p,n in zip(original,native_names):
        q=out/'native_source'/n;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    boost=read(ROOT/'outputs/v203_validation/boost01/report.json')
    headers=[Path(r['path']) for r in boost['sources']]
    for p,r in zip(headers,boost['sources']):
        if sha(p)!=r['sha256']:raise ValueError('Boost准确头改变')
    names=['scripts/buildPhaseAwareAir.py','scripts/phaseAwareDryAir.py','docs/validation/v206/NATIVE_PLAN.md',
        'outputs/v206_validation/reference01/report.json','outputs/v206_validation/python_tests01/report.json']
    locks=[*[ROOT/n for n in names],SOURCE,out/'mixture.dict',*original,*headers]
    sources=[p for d in [library,probe] for p in d.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(model=MODEL,domain=DOMAIN,buildSeconds=100,querySeconds=60,totalSeconds=160,parentSeconds=180,
        memoryMiB=2048,outputMiB=96,threads=1,automaticRetry=0,inputSHA256={str(p):sha(p) for p in locks+sources}))
    lib=run(['wmake','libso'],library,out/'build_library',100,locks+sources,96,2048)
    save(out/'progress.json',dict(library=lib))
    if not lib['passed']:raise RuntimeError('原生库构建失败，保留进程回执')
    helper=run(['wmake'],probe,out/'build_probe',100-lib['elapsedSeconds'],locks+sources,96,2048)
    save(out/'progress.json',dict(library=lib,helper=helper))
    if not helper['passed']:raise RuntimeError('原生查询器构建失败')
    binary=probe/'rahPhaseAwareAirProbe';libpath=library/'librahPhaseAwareAir.so'
    # 模块二：旧p到新T不是最终PT状态，探针显式分开pOld和pFinal。
    cases=[('T',100000,60,300,10),('T',10,60.01,2000,10),('T',10,66.25,2000,10),
        ('T',1000,80,60,1000),('T',1000,99.994,180,1000),('T',1197.0279007383137,226.50908361133003,60,1197.0279007383137),
        ('T',1000000,2000,60,1000000),('T',1000,130.25,60,1000),('T',100000,1200.75,60,100000),
        ('E',69.986503943620136,-227400.12220525689,100,69.986503943620136)]
    sample=out/'positive.txt';sample.write_text('\n'.join(' '.join(str(x) for x in r) for r in cases)+'\n')
    query_start=time.monotonic();positive=run([binary,out/'mixture.dict',sample],out,out/'positive',15,locks+[binary,libpath,sample],96,2048)
    save(out/'progress.json',dict(library=lib,helper=helper,positive=positive))
    if not positive['passed']:raise RuntimeError('正向原生查询失败')
    rows=np.array([[float(v) for v in line.split()[1:]] for line in (out/'positive/stdout.txt').read_text().splitlines() if line.startswith('ROW ')])
    if len(rows)!=len(cases):raise ValueError('正向行数错误')
    target=np.array([r[2] if r[0]=='T' else float(gas.inverse(r[4],r[2],'e')) for r in cases])
    terr=float(abs(rows[:,3]-target).max());eerr=float(np.max(abs(rows[:,4]-rows[:,1])/np.maximum(1,abs(rows[:,1]))))
    values=gas.evaluate(rows[:,5],rows[:,3]);keys=['cpJkgK','cvJkgK','enthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+6]-values[k])/np.maximum(1,abs(values[k])))) for i,k in enumerate(keys)}
    perror=float(np.max(abs(rows[:,12]/pressure_ceiling(rows[:,3])-1)))
    if terr>2e-6 or eerr>2e-10 or max(errors.values())>2e-9 or perror>2e-14:raise ValueError('原生跨表示误差失败')
    elo=float(gas.evaluate(1,60)['internalEnergyJkg']);ehi=float(gas.evaluate(1,2000)['internalEnergyJkg'])
    cases_bad=[('e-low','E 10 %.17g 100 10'%(elo-1),1,'RAH_PHASE_ENERGY_OUTSIDE_TABLE'),
        ('e-high','E 10 %.17g 100 10'%(ehi+1),1,'RAH_PHASE_ENERGY_OUTSIDE_TABLE'),
        ('e-nan','E 10 nan 100 10',1,'RAH_PHASE_INVERSE_INPUT'),('old-p','E -1 0 300 1',1,'RAH_PHASE_INVERSE_INPUT'),
        ('seed-low','T 10 80 59 10',1,'RAH_PHASE_INVERSE_INPUT'),('final-p-zero','T 10 80 100 0',10,'RAH_PHASE_PROBE_FINAL_REJECT'),
        ('final-p-high','T 10 1200 100 1000001',10,'RAH_PHASE_PROBE_FINAL_REJECT'),
        ('final-dew','T 10 60 300 %.17g'%(float(dew_pressure(60))*.5),10,'RAH_PHASE_PROBE_FINAL_REJECT')]
    rejects=[]
    for label,line,code,token in cases_bad:
        sample=out/(label+'.txt');sample.write_text(line+'\n')
        left=60-(time.monotonic()-query_start)
        if left<=0:raise TimeoutError('查询总预算耗尽')
        proc=run([binary,out/'mixture.dict',sample],out,out/label,min(8,left),locks+[binary,libpath,sample],96,2048)
        text=(out/label/'stdout.txt').read_text()
        if proc['nativeExitCode']!=code or token not in text or proc['liveGroupAfterCleanup'] or not all(proc['inputStable'].values()):
            raise ValueError('原生拒绝不符 '+label)
        rejects.append(dict(case=label,expectedExitCode=code,process=proc))
    report=dict(passed=True,model=MODEL,domain=DOMAIN,libraryBuild=lib,probeBuild=helper,positive=positive,
        positiveStates=len(rows),temperatureMaxAbsErrorK=terr,energyMaxScaledError=eerr,propertyMaxScaledErrors=errors,
        pressureCeilingMaxRelative=perror,rejections=rejects,oldPressureSeparatedFromFinalPressure=True,
        v205InputConditionalRootK=float(rows[-1,3]),conditionalRootIsNotSavedPhysicalState=True,
        actualSolverDispatchVerified=False,actualPerStepGuardVerified=False,CFDRuns=0,engineBuilds=0,
        registrationLibrarySHA256=sha(libpath),queryBinarySHA256=sha(binary),tableSHA256=sha(out/'mixture.dict'),
        physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','positiveStates','temperatureMaxAbsErrorK','energyMaxScaledError','v205InputConditionalRootK','elapsedSeconds']}))

if __name__=='__main__':main()
