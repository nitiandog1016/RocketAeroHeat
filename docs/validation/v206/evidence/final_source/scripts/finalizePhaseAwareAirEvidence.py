"""v206只读终验：绑定执行源、保留失败及资源监督缺口，不再次计算物性/CFD。"""
from pathlib import Path
import json,time,shutil
from prepareOpenFoamControl import ROOT,sha,save

BASE=ROOT/'outputs/v206_validation'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def main():
    start=time.monotonic();out=BASE/'final_audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    reference=read(BASE/'reference01/report.json');python=read(BASE/'python_tests01/report.json')
    native=read(BASE/'native_audit02/report.json');rocket=read(BASE/'rocket_audit01/report.json')
    matlab=read(BASE/'matlab02/report.json');matproc=read(BASE/'matlab02/process.json')
    reject=read(BASE/'failure_tests01/report.json');oldmat=read(BASE/'matlab01/process.json')
    for r in [reference,python,native,rocket,matlab,matproc,reject]:
        if not r['passed']:raise ValueError('已有审核未通过')
    if oldmat['passed'] or rocket['candidatePassed']:raise ValueError('历史失败不能改写')
    # 模块一：实际执行副本与原回执绑定。参考副本是事后保存，不能追称事前快照。
    bindings=[]
    for rel,digest in reference['sourceSHA256'].items():
        source=ROOT/rel
        if rel=='scripts/phaseAwareDryAir.py':source=BASE/'reference01/executed_phaseAwareDryAir.py'
        if sha(source)!=digest:raise ValueError('参考执行源无法绑定 '+rel)
        bindings.append(dict(source=source.relative_to(ROOT).as_posix(),sha256=digest))
    for rel,digest in python['sourceSHA256'].items():
        if sha(ROOT/rel)!=digest:raise ValueError('Python接口改变 '+rel)
    for rel,digest in matproc['inputSHA256'].items():
        if sha(ROOT/rel)!=digest:raise ValueError('MATLAB接口改变 '+rel)
    for sub in ['library','probe']:
        directory=ROOT/'native'/('openfoamPhaseAwareAir' if sub=='library' else 'openfoamPhaseAwareAirProbe')
        for p in directory.rglob('*'):
            if p.is_file() and p.name!='README.md':
                if sha(p)!=sha(BASE/'native01'/sub/p.relative_to(directory)):raise ValueError('实际原生源改变')
    protected={p:sha(ROOT/p.replace('\\','/'))==h for p,h in read(ROOT/'outputs/v176_validation/protected_initial.json').items()}
    if len(protected)!=29 or not all(protected.values()):raise ValueError('保护项改变')
    # 模块二：数值验证与计划执行分列，未经在线观测的峰值保持未知。
    resource=dict(numericalReferencePassed=True,resourcePlanAccepted=False,parentTimeoutSeconds=75,
        referenceElapsedSeconds=reference['elapsedSeconds'],onlineMemoryMonitorImplemented=False,
        onlineOutputMonitorImplemented=False,peakRSSBytes=None,rerunToRepairEvidence=False)
    save(out/'reference_resource_gap.json',resource)
    build_seconds=native['libraryBuild']['elapsedSeconds']+native['probeBuild']['elapsedSeconds']
    query_seconds=native['positive']['elapsedSeconds']+sum(r['process']['elapsedSeconds'] for r in native['rejections'])
    result=dict(passed=True,previousGoalTurn='progress',reference=resource,referenceExecutedSourceBindings=bindings,
        referenceSnapshotsCapturedAfterRun=True,PythonStates=python['positiveStates'],PythonRejections=len(python['rejections']),
        nativePositiveStates=native['positiveStates'],nativeRejections=len(native['rejections']),
        nativeBuildSeconds=build_seconds,nativeQuerySupervisionSeconds=query_seconds,
        nativeOriginalPlanPassed=False,nativeRecoveryPassed=True,MATLABOriginalFailed=True,MATLABRecoveryPassed=True,
        matlabTotalSupervisionSeconds=oldmat['elapsedSeconds']+matproc['elapsedSeconds'],
        failureIdentityRejections=len(reject['rejections']),closedCLIVerified=reject['closedCLIVerified'],
        originalWallFaces=rocket['originalWallFaces'],originalAreaM2=rocket['originalAreaM2'],
        CFDRuns=1,CFDCandidatePassed=False,CFDCandidateClosed=True,newEngineBuilds=0,adapterLibraryBuilds=1,probeBuilds=1,
        protected=protected,newQualifiedHeatAreaM2=0,physicalQualification=False,completeGoal=False,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','nativeBuildSeconds','nativeQuerySupervisionSeconds','matlabTotalSupervisionSeconds','elapsedSeconds']}))

if __name__=='__main__':main()
