"""构建独立有界物性映射并以正向/拒绝样本核验，不修改引擎。"""
from pathlib import Path
import json,shutil,time,subprocess
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from nistIntegratedAir import AIR,evaluate


def main():
    from runOpenFoamProcess import run
    start=time.monotonic();base=ROOT/'outputs/v203_validation';out=base/'bounded01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    boost=json.loads((base/'boost01/report.json').read_text());scalar=json.loads((base/'scalar_audit02/report.json').read_text())
    assert boost['passed'] and scalar['passed'] and scalar['falseRejectionOfInRangeEnergy']
    assert subprocess.check_output(['dpkg-query','-W','-f=${Version}','libboost1.83-dev'],text=True)==boost['version']
    library=out/'library';probe=out/'probe'
    shutil.copytree(ROOT/'native/openfoamBoundedNist',library);shutil.copytree(ROOT/'native/openfoamBoundedNistProbe',probe)
    (out/'mixture.dict').write_text(AIR.dictionary())
    source=json.loads((base/'source01/report.json').read_text())
    native_sources=[Path('/opt/openfoam13')/r['path'] for r in source['files']]
    for p,row in zip(native_sources,source['files']):assert sha(p)==row['sha256']
    boosts=[Path(row['path']) for row in boost['sources']]
    for p,row in zip(boosts,boost['sources']):assert sha(p)==row['sha256']
    locks=[Path(__file__),ROOT/'scripts/nistIntegratedAir.py',ROOT/'data/thermo/nist_air_shomate_v1.json',
        ROOT/'docs/validation/v203/BOUNDED_INVERSE_PLAN.md',base/'boost01/report.json',base/'scalar_audit02/report.json',
        *native_sources,*boosts,out/'mixture.dict']
    libsrc=[p for p in library.rglob('*') if p.is_file()];probesrc=[p for p in probe.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(profile='rah-bounded-nist-es-v1',buildSeconds=120,parentSeconds=150,
        positiveSeconds=20,rejectionSecondsEach=10,maximumRejectionCases=8,queryTotalSeconds=100,
        memoryMiB=2048,outputMiB=96,automaticContinuation=0,inputSHA256={str(p):sha(p) for p in locks+libsrc+probesrc}))
    lb=run(['wmake','libso'],library,out/'build_library',100,locks+libsrc,96,2048)
    save(out/'build_progress.json',dict(library=lb));assert lb['passed'],'注册映射构建失败'
    pb=run(['wmake'],probe,out/'build_probe',120-lb['elapsedSeconds'],locks+libsrc+probesrc,96,2048)
    save(out/'build_progress.json',dict(library=lb,probe=pb));assert pb['passed'],'查询器构建失败'
    # 模块一：涵盖真实域内误拒控制、端点/接缝/随机状态和跨温域种子。
    cases=[(1197.0279007383137,row['trueTemperatureK'],row['seedK']) for row in source['cases']]
    rng=np.random.default_rng(203)
    for T in [499.999,500,500.001,699.999,700,700.001,1999.999,2000,2000.001,*rng.uniform(100.001,2499.999,32)]:
        cases.append((float(rng.choice([100,1000,100000])),float(T),float(rng.choice([100,180,1000,2500]))))
    sample=out/'positive.txt';sample.write_text('\n'.join('T %.17g %.17g %.17g'%row for row in cases)+'\n')
    binary=probe/'rahBoundedNistProbe';lib=library/'librahBoundedNist.so'
    q=run([binary,out/'mixture.dict',sample],out,out/'positive',20,locks+[binary,lib,sample],96,2048);assert q['passed'],'正向物性反算失败'
    text=(out/'positive/stdout.txt').read_text();rows=np.array([[float(x) for x in line.split()[1:]] for line in text.splitlines() if line.startswith('ROW ')])
    assert len(rows)==len(cases);target=np.array([x[1] for x in cases]);errT=float(abs(rows[:,3]-target).max())
    errE=float(np.max(abs(rows[:,4]-rows[:,1])/np.maximum(1,abs(rows[:,1]))));assert errT<=2e-5 and errE<=2e-10
    g=evaluate(rows[:,0],rows[:,3]);keys=['cpJkgK','cvJkgK','sensibleEnthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+5]-g[k])/np.maximum(1,abs(g[k])))) for i,k in enumerate(keys)};assert max(errors.values())<=2e-9
    elo=float(evaluate(1197.0279007383137,100)['internalEnergyJkg']);ehi=float(evaluate(1197.0279007383137,2500)['internalEnergyJkg'])
    rejects=[('energy_below',1197.,elo-1,180,'RAH_ES_OUTSIDE_SUPPORTED_RANGE'),('energy_above',1197.,ehi+1,180,'RAH_ES_OUTSIDE_SUPPORTED_RANGE'),
        ('energy_nan',1197.,'nan',180,'RAH_INVALID_THERMO_INPUT'),('pressure_zero',0,elo+1,180,'RAH_INVALID_THERMO_INPUT'),
        ('pressure_high',100001,elo+1,180,'RAH_INVALID_THERMO_INPUT'),('seed_low',1197.,elo+1,99.9,'RAH_INVALID_THERMO_INPUT'),
        ('seed_high',1197.,elo+1,2500.1,'RAH_INVALID_THERMO_INPUT'),('seed_nan',1197.,elo+1,'nan','RAH_INVALID_THERMO_INPUT')]
    rejected=[]
    for name,p,e,seed,token in rejects:
        sample=out/(name+'.txt');sample.write_text('E %s %s %s\n'%(p,e,seed))
        process=run([binary,out/'mixture.dict',sample],out,out/name,10,locks+[binary,lib,sample],96,2048)
        text=(out/name/'stdout.txt').read_text();assert process['nativeExitCode']==1 and token in text
        assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup'] and process['logExportPassed']
        rejected.append(dict(name=name,process=process))
    result=dict(passed=True,profile='rah-bounded-nist-es-v1',libraryBuild=lb,probeBuild=pb,positive=q,positiveStates=len(cases),
        inverseTemperatureMaxAbsK=errT,energyResidualMaxScaled=errE,propertyMaxScaledErrors=errors,rejections=rejected,
        gasModelUnchanged=True,tableSHA256=sha(out/'mixture.dict'),registrationLibrarySHA256=sha(lib),queryBinarySHA256=sha(binary),
        scalarInverseQualified=True,actualSolverDispatchVerified=False,engineBuilds=0,CFDRuns=0,physicalQualification=False,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','positiveStates','inverseTemperatureMaxAbsK','energyResidualMaxScaled','elapsedSeconds']}))


if __name__=='__main__':main()
