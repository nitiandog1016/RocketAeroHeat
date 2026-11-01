"""保存v134实际失败、已测UI和小型证据；核对准确源码，不运行数值或删除文件。"""
from pathlib import Path
import json,hashlib,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v134_validation';DEST=ROOT/'docs/validation/v134'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,d):
    assert not p.exists(),str(p)
    p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))

def main():
    assert not (DEST/'delivery_manifest.json').exists()
    inventory=[];unique={}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        st=p.stat();assert st.st_ino;unique[(st.st_dev,st.st_ino)]=st.st_size
        inventory.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=st.st_size,sha256=sha(p),fileID=st.st_ino,device=st.st_dev))
    write(DEST/'storage_inventory.json',dict(scope='v134 only',files=inventory,pathCount=len(inventory),
        distinctFileIDs=len(unique),pathLengthBytes=sum(x['bytes'] for x in inventory),distinctFileLengthBytes=sum(unique.values()),
        isAllocatedDiskUsage=False,deletedFiles=0,historicalOutputsTouched=False))
    copied=[]
    def copy(p,target):
        target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists();shutil.copyfile(p,target);assert sha(p)==sha(target)
        copied.append(dict(source=p.relative_to(ROOT).as_posix(),archived=target.relative_to(ROOT).as_posix(),bytes=target.stat().st_size,sha256=sha(target)))
    for row in inventory:
        p=ROOT/row['path']
        if row['bytes']<2*1024**2 and p.suffix in ['.json','.log','.png','.mat']:copy(p,DEST/'evidence'/p.relative_to(BASE))
    names=['buildSU2MeshQualityPreview.m','rocketAeroHeatWorkbench.m','web/geometry.js',
        'tests/testSU2MeshQualityWorkbench.m','tests/verifySU2MeshQualityPreview.py',
        'scripts/buildSU2NoseEnvelope.py','scripts/relaxSU2NoseEnvelope.py','scripts/plotSU2MeshQuality.py',
        'scripts/meshSU2HybridVolume.py','scripts/prepareSU2CheckpointContinuation.py','scripts/runMatlabResearch.py',
        'readSU2HybridVolume.m','buildSU2CreaseFanSlice.m','createAerothermalGeometryView.m','updateAerothermalGeometryView.m',
        'web/geometry.html','web/geometry.css','scripts/archiveSU2V134Evidence.py']
    sources=[]
    for name in names:
        p=ROOT/name;target=DEST/'source_versions'/(p.name+'.'+sha(p)+'.raw');copy(p,target)
        sources.append(dict(path=name,sha256=sha(p),archived=target.relative_to(ROOT).as_posix()))
    tasks=[]
    for name in ['mesh_quality_workbench01_process.json','mesh_quality_workbench02_process.json']:
        d=read(BASE/name);assert d['status']=='terminal' and all(d['sourceStableAfter'].values())
        differences=[]
        for path,h in d['sourceSHA256'].items():
            if sha(ROOT/path)==h:continue
            versions=[p for p in (DEST/'failure_versions').glob('*.'+h+'.raw') if sha(p)==h]
            assert len(versions)==1,(path,h)
            differences.append(dict(path=path,sha256=h,archived=versions[0].relative_to(ROOT).as_posix()))
        tasks.append(dict(receipt=name,exitCode=d['exitCode'],lockedFiles=len(d['sourceSHA256']),changedSinceExecution=differences))
    runs=[]
    for name in ['envelope/prepare_process.json','surface/prepare_process.json']:
        p=BASE/name;d=read(p);assert d['status']=='terminal' and d['exitCode']==1 and d['inputStable'] and d['sourceStable']
        for path,h in d['sourceHashes'].items():assert sha(Path(path))==h,path
        for path,h in d['outputSHA256'].items():assert sha(p.parent/path)==h,path
        runs.append(dict(receipt=name,exitCode=1,elapsedSeconds=d['elapsedSeconds']))
    listed={x['archived'] for x in copied}
    for p in sorted(DEST.rglob('*')):
        if p.is_file() and p.relative_to(ROOT).as_posix() not in listed:
            copied.append(dict(archived=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    last=read(BASE/'mesh_quality_workbench02_process.json')
    m=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),files=copied,
        sourceBindings=sources,matlabTasks=tasks,geometryPreparations=runs,nativeSurfaceOptimizationCalls=1,
        native3DMeshGenerations=0,nativeFlowRuns=0,newQualifiedArea=0,defaultChanged=False,
        evidenceBytesBeforeManifest=sum(x['bytes'] for x in copied),postMatlabTools=['tests/verifySU2MeshQualityPreview.py','scripts/plotSU2MeshQuality.py'],
        unmergedResources=[dict(path=n,taskSHA256=last['sourceSHA256'][n],committedByThisVersion=False)
            for n in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py']],
        localOnlyReplayInputs=last['explicitInputs'],largeHistoricalResultsBackedUp=False,
        firstPlotClipped=True,acceptedScientificPlot='evidence/nose_quality_layout02.png',nativeMouseOrScreenshotTest=False)
    write(DEST/'delivery_manifest.json',m)
    print(json.dumps(dict(archivedFiles=len(copied),evidenceBytes=m['evidenceBytesBeforeManifest'],
        outputFileLengthBytes=sum(unique.values()),matlabTasks=tasks),indent=2))

if __name__=='__main__':main()
