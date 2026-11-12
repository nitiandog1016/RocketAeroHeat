"""模块职责：核对原区域索引/监督证据及保护文件，归档小证据，不复制母网格或启动数值任务。"""
from pathlib import Path
import ast,json,shutil,sys,time
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v233_validation';DOC=ROOT/'docs/validation/v233'
OWN=['scripts/phaseEulerRegionalVolume.py','scripts/inspectPhaseRegionalClosureV233.py',
 'scripts/finalizeRegionalClosureV233.py','tests/testPhaseRegionalVolumeV233.py',
 'docs/HYBRID_REGIONAL_CLOSURE_V233.md','docs/validation/v233/PLAN.md','docs/validation/v233/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：基线保护，接受明确保存的先前源码快照，不把开发前后版本混为同一次执行。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256'],row['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    assert set(records)=={'process_sources','process_regions','process_coupling','process_module_test','process_module_test2'}
    older=[]
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None
        assert r['passed']==(name!='process_module_test')
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            assert Path(p) in [ROOT/n for n in OWN],p
            snapshot=BASE/name/'execution_sources'/(h+Path(p).suffix)
            assert sha(snapshot)==h,p;older.append(dict(process=name,path=p,executionSnapshotSHA256=h,currentSHA256=sha(p)))
    a=read(BASE/'region01/report.json');b=read(BASE/'coupling01/report.json');s=read(BASE/'source01/report.json')
    assert a['passed'] and a['exactPartition'] and a['originalWallOwnershipExact']
    assert a['front']['cells']+a['aft']['cells']==64982 and a['front']['components']==a['aft']['components']==1
    assert a['frontInterface']['characteristicFaceCounts']=={'0':378,'1':33,'2':0,'3':0,'4':17,'5':11}
    assert not a['oneWaySupersonicFrontOutletPossible'] and not a['physicalQualification']
    assert b['passed'] and b['allRegionsConnected'] and b['allOriginalCellsAndWallsPreserved']
    assert [r['cells'] for r in b['regions']]==[160304,545186] and b['interfaceFaces']==1984
    assert not b['nativeMultiZoneTested'] and not b['directEulerRANSReady'] and b['newQualifiedAreaM2']==0
    for row in s['sources']+b['sourceSHA256']:assert sha(row['path'])==row['sha256']
    r=np.loadtxt(str(BASE/'region01/wall_responsibility.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(r[:,0],np.arange(1,19753)) and set(r[:,2])=={0.,1.}
    assert abs(float(np.sum(r[:,1]))-1.9646292565715553)<1e-12
    z=np.load(str(BASE/'coupling01/viscous_two_zone_index.npz'))
    assert np.array_equal(np.sort(np.r_[z['upstreamCells'],z['aftCells']]),np.arange(705490))
    assert np.array_equal(np.sort(np.r_[z['upstreamOriginalFaces'],z['aftOriginalFaces']]),np.arange(1,19753))
    assert np.array_equal(np.flatnonzero(r[:,2]==1)+1,z['upstreamOriginalFaces'])
    # 模块二：发布小CSV/报告/调用源，完整索引保留本地并登记可重建身份。
    out.mkdir();evidence.mkdir()
    for folder in BASE.glob('process_*'):shutil.copytree(folder,evidence/folder.name)
    for stage in ('source01','region01','coupling01'):
        target=evidence/stage;target.mkdir()
        for p in (BASE/stage).iterdir():
            if p.suffix!='.npz':shutil.copyfile(p,target/p.name)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.glob('*/*.npz')}
    save(evidence/'local_asset_index.json',dict(base='outputs/v233_validation',files=local,
      rebuild='inspectPhaseRegionalClosureV233.py regions/coupling in a fresh version directory; existing stages do not overwrite'))
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    seconds=sum(r['elapsedSeconds'] for r in records.values());elapsed=time.monotonic()-started
    assert seconds+elapsed<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    report=dict(passed=True,goalStatus='active',scientificAcceptance=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
      EulerRegionalCells=[17086,47896],viscousRegionalCells=[160304,545186],originalWallFaces=[4434,15318],
      originalAreaM2=1.9646292565715553,EulerInterfaceFaces=439,EulerInwardCharacteristicFaces=61,
      viscousInterfaceFaces=1984,allRegionsConnected=True,originalCellsAndWallsPreserved=True,
      sourceFiles=len(set(row['path'] for row in s['sources']+b['sourceSHA256'])),
      directEulerRANSReady=False,dualRANSSelectedNotNativeTested=True,nativeBuilds=0,nativeQueries=0,flowUpdates=0,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      processSeconds=seconds,finalizationSeconds=elapsed,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      outputBytesBeforeReport=size,allGroupsEmpty=True,first32BitTestFailurePreserved=True,
      earlierExecutedSourcesBoundBySnapshots=older,unmeteredTargetedReadsDisclosed=True,
      MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,uniqueAssetsDeleted=False,newMeshGenerated=False)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v233/evidence_manifest.json']+['docs/validation/v233/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),processSeconds=seconds,finalizationSeconds=elapsed,bytes=size)))


if __name__=='__main__':main()
