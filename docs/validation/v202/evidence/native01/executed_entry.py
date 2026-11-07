"""有预算地构建原生模板注册组合/物性查询器，并验证同表状态及越界拒绝。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from nistIntegratedAir import AIR,evaluate,contract,R

BASE=ROOT/'outputs/v202_validation'


def main():
    from runOpenFoamProcess import run
    start=time.monotonic();out=BASE/'native01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    review=BASE/'source_review02/report.json';offline=BASE/'model_check01/report.json'
    for p in [review,offline]:assert json.loads(p.read_text())['passed']
    records=json.loads(review.read_text())['files']+json.loads((BASE/'source_supplement01/sources.json').read_text())
    source_locks=[]
    for row in records:
        p=Path('/opt/openfoam13')/row['path'];assert sha(p)==row['sha256'];source_locks.append(p)
    source_locks.append(Path('/opt/openfoam13/etc/controlDict'))
    library=out/'library';probe=out/'probe'
    shutil.copytree(ROOT/'native/openfoamNistAir',library);shutil.copytree(ROOT/'native/openfoamNistProbe',probe)
    (out/'mixture.dict').write_text(AIR.dictionary());save(out/'thermo_contract.json',contract())
    locks=[Path(__file__),ROOT/'scripts/nistIntegratedAir.py',ROOT/'data/thermo/nist_air_shomate_v1.json',
        ROOT/'docs/validation/v202/MODEL_PLAN.md',review,offline,*source_locks]
    libsrc=[p for p in library.rglob('*') if p.is_file()];probesrc=[p for p in probe.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(buildSeconds=120,parentSeconds=150,cpuCount=1,memoryMiB=2048,outputMiB=96,
        querySeconds=20,boundaryRejectionSecondsEach=10,automaticContinuation=0,inputSHA256={str(p):sha(p) for p in locks+libsrc+probesrc}))
    lib=run(['wmake','libso'],library,out/'build_library',100,locks+libsrc,96,2048)
    save(out/'build_progress.json',dict(library=lib));assert lib['passed'],'原生组合注册库构建失败'
    remaining=120-lib['elapsedSeconds'];assert remaining>0
    helper=run(['wmake'],probe,out/'build_probe',remaining,locks+probesrc,96,2048)
    save(out/'build_progress.json',dict(library=lib,probe=helper));assert helper['passed'],'原生物性查询器构建失败'
    rng=np.random.default_rng(202);ts=np.r_[100,100.01,150,180,199.99,200,226.50908361133003,298.15,499.999,500,500.001,699.999,700,700.001,999.999,1000,1000.001,1999.999,2000,2000.001,2499.99,2500,rng.uniform(101,2499,36)]
    ps=np.resize(np.array([100.,1197.0279007383137,5000.,100000.]),len(ts))
    samples=out/'samples.txt';samples.write_text('\n'.join('%.17g %.17g %.17g'%(p,t,t if t in (100,2500) else min(2499,max(101,t*1.000001))) for p,t in zip(ps,ts))+'\n')
    binary=probe/'rahNistProbe';libraryfile=library/'librahNistAir.so'
    q=run([binary,out/'mixture.dict',samples],out,out/'query',20,locks+[binary,libraryfile,out/'mixture.dict',samples],96,2048)
    assert q['passed'];text=(out/'query/stdout.txt').read_text();assert 'FOAM FATAL' not in text
    meta=np.array([float(x) for x in next(s for s in text.splitlines() if s.startswith('META ')).split()[1:]])
    assert np.allclose(meta,[R,0,298.15,100000],rtol=0,atol=2e-12)
    rows=np.array([[float(x) for x in s.split()[1:]] for s in text.splitlines() if s.startswith('ROW ')]);assert len(rows)==len(ts)
    g=evaluate(rows[:,0],rows[:,1]);keys=['cpJkgK','cvJkgK','enthalpyJkg','sensibleEnthalpyJkg','internalEnergyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(abs(rows[:,i+2]-g[k]).max()/max(1,abs(g[k]).max())) for i,k in enumerate(keys)}
    # 每个状态也独立缩放，避免高温大值掩盖低温/参考点错误。
    local_errors={k:float(np.max(abs(rows[:,i+2]-g[k])/np.maximum(1,abs(g[k])))) for i,k in enumerate(keys)}
    inverse=float(abs(rows[:,-1]-rows[:,1]).max())
    assert max(local_errors.values())<=2e-9 and inverse<=2e-5
    rejections=[]
    for name,t in [('below',99.99),('above',2500.01)]:
        sample=out/(name+'.txt');sample.write_text('5000 %.17g %.17g\n'%(t,t))
        r=run([binary,out/'mixture.dict',sample],out,out/('reject_'+name),10,locks+[binary,sample,out/'mixture.dict'],96,2048)
        log=(out/('reject_'+name)/'stdout.txt').read_text()
        assert r['nativeExitCode']!=0 and 'out of range' in log and not r['liveGroupAfterCleanup'] and all(r['inputStable'].values())
        rejections.append(r)
    result=dict(passed=True,libraryBuild=lib,probeBuild=helper,query=q,rejections=rejections,nativeRows=rows.tolist(),
        propertyMaxScaledErrors=local_errors,inverseTemperatureMaxAbsK=inverse,contract=contract(),
        registrationLibrarySHA256=sha(libraryfile),queryBinarySHA256=sha(binary),engineBuilds=0,adapterBuilds=1,helperBuilds=1,
        actualSolverDispatchVerified=False,CFDRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','propertyMaxScaledErrors','inverseTemperatureMaxAbsK','elapsedSeconds']}))


if __name__=='__main__':main()
