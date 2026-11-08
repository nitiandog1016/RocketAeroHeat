"""只读恢复联合域原生审核，保留首版错误预期；仅执行已登记的一个补充拒绝。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from standardDryAir import StandardDryAir,MODEL,DOMAIN,pressure_ceiling

def main():
    from runOpenFoamProcess import run
    start=time.monotonic();base=ROOT/'outputs/v204_validation';native=base/'native01';out=base/'native_audit02';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');reports={p.parent.name:json.loads(p.read_text()) for p in native.glob('*/report.json')}
    under=(native/'rho-underflow/stdout.txt').read_text();assert 'RAH_STANDARD_INVALID_INPUT e=-inf' in under and reports['rho-underflow']['nativeExitCode']==1
    sources=[]
    for name in ['HtoEthermo.H','perfectGasI.H']:
        p=(Path('/opt/openfoam13/src/thermophysicalModels/specie/lnInclude')/name).resolve();q=out/name;shutil.copyfile(p,q)
        sources.append(dict(path=str(p),sha256=sha(p),lines=[dict(line=i+1,text=x) for i,x in enumerate(p.read_text().splitlines()) if 'rho(' in x or 'return p/' in x]))
    assert 'p/EquationOfState::rho(p, T)' in (out/'HtoEthermo.H').read_text()
    assert 'p/(this->R()*T)' in (out/'perfectGasI.H').read_text()
    for label,r in reports.items():
        assert all(r['inputStable'].values()) and not r['liveGroupAfterCleanup'] and r['logExportPassed'],label
    rows=np.array([[float(x) for x in line.split()[1:]] for line in (native/'positive/stdout.txt').read_text().splitlines() if line.startswith('ROW ')])
    model=StandardDryAir();g=model.evaluate(rows[:,0],rows[:,3]);targets=np.array([float(s.split()[2]) for s in (native/'positive.txt').read_text().splitlines()])
    errT=float(np.max(abs(rows[:,3]-targets)));errE=float(np.max(abs(rows[:,4]-rows[:,1])/np.maximum(1,abs(rows[:,1]))))
    keys=['cpJkgK','cvJkgK','enthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+5]-g[k])/np.maximum(1,abs(g[k])))) for i,k in enumerate(keys)}
    assert errT<=2e-6 and errE<=2e-10 and max(errors.values())<=2e-9
    assert np.array_equal(rows[:,11],pressure_ceiling(rows[:,3]))
    # 模块一：先保存原计划失败状态，后续恢复绝不覆盖这个判断。
    save(native/'report.json',dict(passed=False,originalPlanAcceptance=False,reason='underflow rejected before expected final-density branch; tenth exact-PT one-ulp rejection was not run',
        positiveStates=17,executedRejections=9,physicalQualification=False,CFDRuns=0,engineBuilds=0))
    binary=native/'probe/rahStandardDryAirProbe';library=native/'library/librahStandardDryAir.so';table=native/'mixture.dict';sample=out/'resolved-pressure-gap.txt'
    sample.write_text('T 10000.01 100 100\n')
    oldtimes=sum(r['elapsedSeconds'] for k,r in reports.items() if k not in {'build_library','build_probe'})
    assert oldtimes<50
    process=run([binary,table,sample],out,out/'supplement',min(10,60-oldtimes),[binary,library,table,sample,Path(__file__),ROOT/'scripts/standardDryAir.py',ROOT/'docs/validation/v204/NATIVE_AUDIT_RECOVERY.md'],96,2048)
    assert process['nativeExitCode']==1 and 'RAH_STANDARD_PT_OUTSIDE_DOMAIN' in (out/'supplement/stdout.txt').read_text()
    assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup'] and process['logExportPassed']
    result=dict(passed=True,originalPlanAcceptance=False,mode='post-run-read-only-audit-plus-one-preregistered-rejection',model=MODEL,domain=DOMAIN,
        positiveStates=len(rows),executedRejections=10,priorRejectionExpectedTokenMatched=8,underflowRejectedEarlier=True,
        finalDensityFailureBranchExercised=False,originalOneUlpEnergyCaseNotRun=True,pressureBoundaryBitwiseEqual=True,
        temperatureMaxAbsErrorK=errT,energyMaxScaledResidual=errE,propertyMaxScaledErrors=errors,
        buildSeconds=sum(reports[k]['elapsedSeconds'] for k in ['build_library','build_probe']),querySeconds=oldtimes+process['elapsedSeconds'],
        sourceFiles=sources,supplement=process,registrationLibrarySHA256=sha(library),queryBinarySHA256=sha(binary),tableSHA256=sha(table),
        nativeTemplateVerified=True,actualSolverDispatchVerified=False,newQualifiedHeatAreaM2=0,physicalQualification=False,
        CFDRuns=0,engineBuilds=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','originalPlanAcceptance','temperatureMaxAbsErrorK','energyMaxScaledResidual','buildSeconds','querySeconds','elapsedSeconds']}))

if __name__=='__main__':main()
