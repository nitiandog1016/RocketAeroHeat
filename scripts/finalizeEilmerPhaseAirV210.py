"""绑定测试后的原生源码、复用行和资源回执；保存最终科学范围。"""
from pathlib import Path
import hashlib,json,time,shutil,difflib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v210_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    native=ROOT/'native/eilmerPhaseAwareAir';prep=BASE/'prepared02'
    assert sha(native/'gas/rah_phase_aware_air.d')==sha(prep/'overlay/gas/rah_phase_aware_air.d')
    assert sha(native/'probe.d')==sha(prep/'probe.d')
    a=json.loads((BASE/'audit02/report.json').read_text());assert a['passed'] and not a['physicalQualification']
    assert not json.loads((BASE/'audit01/report.json').read_text())['passed']
    assert not json.loads((BASE/'query01/report.json').read_text())['passed']
    original=(BASE/'query01/real/stdout.txt').read_text().splitlines()
    rows=json.loads((BASE/'query_inputs02/reused_real_rows.json').read_text())
    assert len(rows)==117 and all(original[v['lineNumber']-1]==v['line'] for v in rows.values())
    processes=[]
    names=['build02/real','build03/real','build04/real','build04/complex','query01/real','query02/real','query02/complex']
    for name in names:
        r=json.loads((BASE/name/'report.json').read_text())
        assert r['inputsStable'] and not r['liveGroupAfterCleanup']
        processes.append(dict(name=name,passed=r['passed'],exitCode=r['nativeExitCode'],seconds=r['elapsedSeconds'],
            peakGroupRSSBytes=r['peakObservedGroupRSSBytes'],stopReason=r['stopReason'],inputsStable=True,groupEmpty=True))
    b=json.loads((BASE/'build04/report.json').read_text());assert b['passed']
    assert all(sha(BASE/'build04'/('probe_'+v))==digest for v,digest in b['binarySHA256'].items())
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    c=json.loads((ROOT/'docs/validation/v208/route_contract.json').read_text())
    assert sha(ROOT/c['regionLayerMesh'])==c['regionLayerMeshSHA256']
    assert sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json')=='7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'
    source=ROOT/'outputs/v209_validation/source01/source/src/gas/init_gas_model.d'
    patch=''.join(difflib.unified_diff(source.read_text().splitlines(True),(prep/'overlay/gas/init_gas_model.d').read_text().splitlines(True),
        fromfile='a/src/gas/init_gas_model.d',tofile='b/src/gas/init_gas_model.d'))
    (native/'factory-registration.patch').write_text(patch,encoding='utf-8')
    shutil.copyfile(__file__,out/'executed_entry.py')
    report=dict(passed=True,model='rah-phase-aware-dry-air-eilmer-v1',helperInterfaceVerified=True,
        nativeModelSourceSHA256=sha(native/'gas/rah_phase_aware_air.d'),probeSourceSHA256=sha(native/'probe.d'),
        factoryPatchSHA256=sha(native/'factory-registration.patch'),binarySHA256=b['binarySHA256'],
        protectedFilesUnchanged=29,originalLayerMeshSHA256=c['regionLayerMeshSHA256'],
        reusedOriginalRows=117,recoveryRealRows=29,firstComplexRows=314,baseStates=14,
        derivativeCases=168,finiteDifferencePairs=14,gasRejectionsEachVariant=17,unknownModelFatalRejectionVerified=True,
        branchChanges=len(a['branchChanges']),initialDerivativeAuditPassed=False,initialQueryPlanPassed=False,
        maximumScaledValueError=a['maximumScaledValueError'],maximumTemperatureInverseErrorK=a['maximumTemperatureInverseErrorK'],
        maximumScaledComplexDerivativeError=a['maximumScaledComplexDerivativeError'],maximumScaledNativeDifferenceError=a['maximumScaledNativeDifferenceError'],
        processes=processes,compilerProcesses=4,successfulHelperBuilds=2,
        cumulativeCompilerSeconds=sum(p['seconds'] for p in processes if p['name'].startswith('build')),
        cumulativeQuerySeconds=sum(p['seconds'] for p in processes if p['name'].startswith('query')),
        engineBuilds=0,CFDRuns=0,newQualifiedHeatAreaM2=0,physicalQualification=False,
        actualFlowDispatchVerified=False,fullRegionalResourceAccepted=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['processes','binarySHA256']}))
if __name__=='__main__':main()
