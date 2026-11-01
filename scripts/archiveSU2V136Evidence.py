"""归档v136真实自交检查、两个外域候选和只读UI证据；大网格/运行时保持本地。"""
from pathlib import Path
import json,hashlib,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v136_validation';DEST=ROOT/'docs/validation/v136'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,value):
    assert not p.exists(),str(p);p.write_bytes((json.dumps(value,indent=2,ensure_ascii=False)+'\n').encode('utf8'))
def main():
    assert not (DEST/'delivery_manifest.json').exists();rows=[];unique={}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        s=p.stat();assert s.st_ino;unique[(s.st_dev,s.st_ino)]=s.st_size
        rows.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=s.st_size,fileID=s.st_ino,device=s.st_dev,sha256=sha(p)))
    runtime=ROOT/'tools/runtime/cgal/6.0.1-py311';runtimeRows=[(p,p.stat()) for p in runtime.rglob('*') if p.is_file()]
    runtimeBytes=sum(s.st_size for p,s in runtimeRows)
    write(DEST/'storage_inventory.json',dict(scope='v136 only; runtime separate',files=rows,pathCount=len(rows),distinctFileIDs=len(unique),
        pathLengthBytes=sum(r['bytes'] for r in rows),distinctFileLengthBytes=sum(unique.values()),
        runtimePathLengthBytes=runtimeBytes,runtimeFileCount=len(runtimeRows),isAllocatedDiskUsage=False,deletedFiles=0,historicalOutputsTouched=False))
    copied=[]
    def copy(p,target):
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():assert sha(p)==sha(target)
        else:shutil.copyfile(p,target)
        assert sha(p)==sha(target);copied.append(dict(source=p.relative_to(ROOT).as_posix(),archived=target.relative_to(ROOT).as_posix(),bytes=target.stat().st_size,sha256=sha(target)))
    for row in rows:
        p=ROOT/row['path']
        if row['bytes']<1024**2 and (p.suffix in ['.json','.log','.mat'] or (p.suffix=='.off' and row['bytes']<2048)):
            copy(p,DEST/'evidence'/p.relative_to(BASE))
    runtimeInfo=read(runtime/'runtime.json')
    for name,row in runtimeInfo['files'].items():assert sha(runtime/name)==row['sha256']
    for name in ['runtime.json','import.log']+runtimeInfo['licenseLocations']:
        copy(runtime/name,DEST/'runtime_evidence'/name)
    names=['scripts/installCGALGeometryRuntime.py','scripts/cgalSurfaceIntersections.py','scripts/runCGALSurfaceAudit.py',
        'scripts/reviewSU2TransitionShells.py','scripts/meshSU2ExplicitTransition.py','scripts/reviewSU2ExplicitQuality.py',
        'scripts/meshSU2HybridTransition.py','scripts/meshSU2HybridVolume.py','scripts/su2PyramidGeometry.py',
        'scripts/auditSU2HybridTopology.py','scripts/optimizeSU2FanExterior.py','scripts/prepareSU2CheckpointContinuation.py',
        'scripts/runMatlabResearch.py','scripts/archiveSU2V136Evidence.py','inspectSU2MeshCandidate.m','rocketAeroHeatWorkbench.m',
        'buildSU2MeshQualityPreview.m','tests/testCGALSurfaceIntersections.py','tests/testSU2CandidateInspection.m','tests/verifySU2CandidateInspection.py']
    sources=[]
    for name in names:
        p=ROOT/name;target=DEST/'source_versions'/(p.name+'.'+sha(p)+'.raw');copy(p,target)
        sources.append(dict(path=name,sha256=sha(p),archived=target.relative_to(ROOT).as_posix()))
    def version(path,h):
        if sha(path)==h:return None
        found=[p for p in (DEST/'source_versions').glob('*.'+h+'.raw') if sha(p)==h];assert len(found)==1,(str(path),h)
        return found[0].relative_to(ROOT).as_posix()
    runs=[]
    for name in ['explicit','hxt']:
        p=BASE/name/'mesh_process.json';r=read(p);assert r['status']=='terminal' and r['exitCode']==1 and r['sourceStable'] and r['inputStable']
        changes=[]
        for source,h in r['sourceHashes'].items():
            old=version(Path(source),h)
            if old:changes.append(dict(path=source,sha256=h,archived=old))
        for f,h in r['outputSHA256'].items():assert sha(p.parent/f)==h
        runs.append(dict(case=name,exitCode=1,seconds=r['elapsedSeconds'],changedSinceExecution=changes))
    geometry=[]
    for p in sorted(BASE.rglob('process.json')):
        r=read(p);assert r['status']=='terminal' and r['exitCode']==0 and r['runtimeStableAfter'] and all(r['sourceStableAfter'].values())
        for f,h in r['sourceSHA256'].items():assert sha(Path(f))==h
        for f,h in r['outputSHA256'].items():assert sha(p.parent/f)==h
        geometry.append(dict(receipt=p.relative_to(ROOT).as_posix(),seconds=r['seconds'],exitCode=0))
    task=read(BASE/'v136_inspector01_process.json');assert task['status']=='terminal' and task['exitCode']==0 and all(task['sourceStableAfter'].values())
    for path,h in task['sourceSHA256'].items():assert sha(ROOT/path)==h,path
    existing={r['archived'] for r in copied}
    for p in DEST.rglob('*'):
        if p.is_file() and p.relative_to(ROOT).as_posix() not in existing:copied.append(dict(archived=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    manifest=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),files=copied,sourceBindings=sources,
        nativeMeshRuns=runs,actualCGALProcesses=geometry,matlabTask=dict(exitCode=0,lockedFiles=len(task['sourceSHA256']),sourceStable=True),
        nativeFlowRuns=0,newQualifiedArea=0,meshCandidatesAccepted=0,nativeMouseOrScreenshotTest=False,
        evidenceBytesBeforeManifest=sum(r['bytes'] for r in copied),runtimeBinaryUploaded=False,largeHistoricalResultsBackedUp=False,
        localOnlyCurrentFiles=[row for row in rows if row['bytes']>=1024**2],
        unrelatedChangesNotCommitted=['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py','docs/research/',
            'docs/validation/opencode_feature_retention/','docs/validation/opencode_section_export/','tests/testRocketFeatureRetentionMatrix.m','tests/testSectionExportReasons.js'])
    write(DEST/'delivery_manifest.json',manifest)
    print(json.dumps(dict(files=len(copied),evidenceBytes=manifest['evidenceBytesBeforeManifest'],outputPathBytes=sum(r['bytes'] for r in rows),
        distinctOutputBytes=sum(unique.values()),runtimeBytes=runtimeBytes,CGALTasks=len(geometry),matlabLockedFiles=len(task['sourceSHA256'])),indent=2))
if __name__=='__main__':main()
