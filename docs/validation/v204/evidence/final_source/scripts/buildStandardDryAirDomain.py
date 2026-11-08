"""新具名气体的独立原生注册/查询；只做限定联合域和表示核验，无CFD。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from standardDryAir import StandardDryAir,MODEL,DOMAIN,pressure_ceiling

def main():
    raise ValueError('v204原生构建和查询已执行；首审核失败及独立恢复均保留，不重建或重跑。')
    from runOpenFoamProcess import run
    start=time.monotonic();base=ROOT/'outputs/v204_validation';out=base/'native01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    model=StandardDryAir();offline=json.loads((base/'standard01/report.json').read_text());assert offline['passed']
    final=json.loads((base/'domain_tests03/report.json').read_text());assert final['passed'] and final['finalSourceReferenceRowsPassed']==318
    library=out/'library';probe=out/'probe';shutil.copytree(ROOT/'native/openfoamStandardDryAir',library);shutil.copytree(ROOT/'native/openfoamStandardDryAirProbe',probe)
    (out/'mixture.dict').write_text(model.dictionary())
    boost=json.loads((ROOT/'outputs/v203_validation/boost01/report.json').read_text());sources=json.loads((ROOT/'outputs/v203_validation/source01/report.json').read_text())
    native=[Path('/opt/openfoam13')/row['path'] for row in sources['files']]
    for p,row in zip(native,sources['files']):assert sha(p)==row['sha256']
    headers=[Path(row['path']) for row in boost['sources']]
    for p,row in zip(headers,boost['sources']):assert sha(p)==row['sha256']
    locks=[Path(__file__),ROOT/'scripts/standardDryAir.py',ROOT/'data/thermo/standard_dry_air_cp_v1.json',
        ROOT/'docs/validation/v204/NATIVE_DOMAIN_PLAN.md',base/'standard01/report.json',base/'domain_tests03/report.json',out/'mixture.dict',*native,*headers]
    libsrc=[p for p in library.rglob('*') if p.is_file()];probesrc=[p for p in probe.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(model=MODEL,domain=DOMAIN,buildSeconds=100,querySeconds=60,totalSeconds=160,parentSeconds=180,
        threads=1,memoryMiB=2048,outputMiB=96,automaticContinuation=0,inputSHA256={str(p):sha(p) for p in locks+libsrc+probesrc}))
    lb=run(['wmake','libso'],library,out/'build_library',100,locks+libsrc,96,2048);save(out/'progress.json',dict(library=lb));assert lb['passed']
    pb=run(['wmake'],probe,out/'build_probe',100-lb['elapsedSeconds'],locks+libsrc+probesrc,96,2048);save(out/'progress.json',dict(library=lb,probe=pb));assert pb['passed']
    cases=[(10000,100,2000),(1000,100.01,180),(1000,298.15,2000),(1e6,2000,100),(1197.0279007383137,226.50908361133003,100),
        (100292.29488236419,1231.3523721818547,180),(702078.7681721675,1274.08470223713,180)]
    for T in [132.53,200.,250.,464.1588833612779]:cases.append((float(pressure_ceiling(T)),T,1000))
    for T in [100.02,150.,500.25,999.8,1500.25,1999.99]:cases.append((float(pressure_ceiling(T))*.99,T,2000 if T<1000 else 100))
    sample=out/'positive.txt';sample.write_text('\n'.join('T %.17g %.17g %.17g'%row for row in cases)+'\n')
    binary=probe/'rahStandardDryAirProbe';lib=library/'librahStandardDryAir.so'
    qstart=time.monotonic();q=run([binary,out/'mixture.dict',sample],out,out/'positive',15,locks+[binary,lib,sample],96,2048)
    save(out/'progress.json',dict(library=lb,probe=pb,positive=q));assert q['passed']
    rows=np.array([[float(v) for v in line.split()[1:]] for line in (out/'positive/stdout.txt').read_text().splitlines() if line.startswith('ROW ')])
    assert len(rows)==len(cases);Terror=float(np.max(abs(rows[:,3]-np.array([r[1] for r in cases]))))
    Eerror=float(np.max(abs(rows[:,4]-rows[:,1])/np.maximum(1,abs(rows[:,1]))));assert Terror<=2e-6 and Eerror<=2e-10
    d=model.evaluate(rows[:,0],rows[:,3]);names=['cpJkgK','cvJkgK','enthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+5]-d[k])/np.maximum(1,abs(d[k])))) for i,k in enumerate(names)}
    assert max(errors.values())<=2e-9
    assert np.array_equal(rows[:,11],pressure_ceiling(rows[:,3])),'原生/Python压力边界不同'
    elo=float(model.evaluate(1000,100)['internalEnergyJkg']);ehi=float(model.evaluate(1000,2000)['internalEnergyJkg'])
    rejects=[('energy-low','E',1000,elo-1,100,'RAH_STANDARD_ENERGY_OUTSIDE_TABLE'),('energy-high','E',1000,ehi+1,100,'RAH_STANDARD_ENERGY_OUTSIDE_TABLE'),
        ('energy-nan','E',1000,'nan',100,'RAH_STANDARD_INVALID_INPUT'),('p-zero','T',0,300,300,'RAH_STANDARD_INVALID_INPUT'),
        ('p-cap','T',1000001,1000,300,'RAH_STANDARD_INVALID_INPUT'),('cold-pressure','T',100000,100,100,'RAH_STANDARD_PT_OUTSIDE_DOMAIN'),
        ('seed-low','T',1000,300,99,'RAH_STANDARD_INVALID_INPUT'),('seed-high','T',1000,300,2001,'RAH_STANDARD_INVALID_INPUT'),
        ('rho-underflow','T',5e-324,300,300,'RAH_STANDARD_NUMERICAL_RANGE'),('one-ulp-pressure','T',float(np.nextafter(10000.,np.inf)),100,100,'RAH_STANDARD_PT_OUTSIDE_DOMAIN')]
    rejected=[]
    for label,mode,p,e,seed,token in rejects:
        assert time.monotonic()-qstart<60
        sample=out/(label+'.txt');sample.write_text('%s %s %s %s\n'%(mode,p,e,seed))
        proc=run([binary,out/'mixture.dict',sample],out,out/label,min(10,60-(time.monotonic()-qstart)),locks+[binary,lib,sample],96,2048)
        text=(out/label/'stdout.txt').read_text();assert proc['nativeExitCode']==1 and token in text,label
        assert all(proc['inputStable'].values()) and not proc['liveGroupAfterCleanup'] and proc['logExportPassed']
        rejected.append(dict(case=label,process=proc))
    assert time.monotonic()-start<160
    result=dict(passed=True,model=MODEL,domain=DOMAIN,libraryBuild=lb,probeBuild=pb,positive=q,positiveStates=len(cases),
        rejections=rejected,temperatureErrorMaxK=Terror,energyResidualMaxScaled=Eerror,propertyErrors=errors,
        pressureBoundaryBitwiseEqual=True,registrationLibrarySHA256=sha(lib),queryBinarySHA256=sha(binary),tableSHA256=sha(out/'mixture.dict'),
        nativeTemplateVerified=True,actualSolverDispatchVerified=False,engineBuilds=0,CFDRuns=0,physicalQualification=False,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','positiveStates','temperatureErrorMaxK','energyResidualMaxScaled','elapsedSeconds']}))

if __name__=='__main__':main()
