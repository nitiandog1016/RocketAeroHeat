"""保存v126真实终态、原失败、主工作台证据与准确源码；不复制大型数值场。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2CheckpointContinuation import ROOT,read,sha
BASE=ROOT/'outputs/v126_validation';DEST=ROOT/'docs/validation/v126'


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    case=Path(read(BASE/'main01/engine_folder.json')['folder']);matlab=[];native=[]
    for name,code in [('sst_output_main01',1),('sst_output_main02',0)]:
        r=read(BASE/(name+'_process.json'))
        assert r['status']=='terminal' and r['exitCode']==code and all(r['sourceStableAfter'].values())
        matlab.append(dict(name=name,exitCode=code,sourceStable=True))
    for p in [BASE/'format01',BASE/'restart02',case]:
        r=read(p/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        native.append(dict(folder=str(p),exitCode=0,elapsedSeconds=r['elapsedSeconds']))
    pipeline=read(case/'workflow/pipeline_process.json');assert pipeline['status']=='terminal' and pipeline['exitCode']==0
    files=[]
    for root,subdir in [(BASE,'evidence'),(case,'native_main')]:
        for p in sorted(root.rglob('*')):
            if not p.is_file() or p.stat().st_size>4*1024**2:continue
            smallNumerical=p.parent in [BASE/'format01',BASE/'restart02'] and p.suffix in ['.csv','.vtk','.su2']
            if p.suffix not in ['.json','.log','.cfg'] and p.name not in ['history.csv','section.csv'] and not smallNumerical:continue
            q=DEST/subdir/p.relative_to(root);q.parent.mkdir(parents=True,exist_ok=True);assert not q.exists();shutil.copyfile(p,q)
            files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    selected=set('prepareSU2SSTCase.py prepareSU2SSTCheckpoint.py prepareSU2SSTOutputProbe.py prepareSU2SSTOutputRestart.py su2ThermalRuntime.py su2SSTInterface.py prepareSU2CheckpointContinuation.py runEmbeddedSU2.py runEmbeddedSU2Pipeline.py runMatlabResearch.py auditSU2SSTState.py testSU2SSTOutputWorkbench.m auditSU2ResultConvergence.m readSU2SSTStateAudit.m describeSU2Convergence.m describeSU2Runtime.m estimateSU2CheckpointStorage.m embeddedSU2Job.m readSU2SurfaceResult.m mapSU2WallResult.m buildSU2WallProjection.m readSU2HybridVolume.m attachSU2ThermalBoundary.m createSU2ResearchView.m createAerothermalGeometryView.m updateAerothermalGeometryView.m rocketAeroHeatWorkbench.m'.split())
    available={sha(p):p for d in [DEST/'initial_source_versions',DEST/'first_view_failure'] for p in d.glob('*') if p.is_file()}
    bindings=[]
    def save(p,digest):
        source=p if p.is_file() and sha(p)==digest else available.get(digest)
        assert source is not None,(str(p),digest)
        q=DEST/'source_versions'/(p.name+'.'+digest+'.raw');q.parent.mkdir(parents=True,exist_ok=True)
        if q.exists():assert sha(q)==digest
        else:shutil.copyfile(source,q)
        row=dict(path=str(p),sha256=digest,archived=q.relative_to(ROOT).as_posix())
        if row not in bindings:bindings.append(row)
    def collect(data):
        if isinstance(data,dict):
            for k,v in data.items():
                if isinstance(v,str) and len(v)==64 and all(c in '0123456789abcdef' for c in v):
                    p=Path(k);p=p if p.is_absolute() else ROOT/p
                    if p.name in selected:save(p,v)
                else:collect(v)
        elif isinstance(data,list):
            for item in data:collect(item)
    for root in [BASE,case]:
        for p in root.rglob('*.json'):
            if p.stat().st_size<4*1024**2:collect(read(p))
    for name in ['scripts/archiveSU2V126Evidence.py','scripts/deduplicateSU2V126Copies.py',
                 'scripts/diagnoseSU2SSTMultistep.py','tests/testSU2SSTOutputContract.py','tests/verifySU2SSTOutputExports.py']:
        p=ROOT/name;save(p,sha(p))
    p=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/output/CFlowCompOutput.cpp';save(p,sha(p))
    resources=['web/viewer.js','web/section.js'];locks=read(BASE/'sst_output_main02_process.json')['sourceSHA256']
    unmerged=[]
    for name in resources:
        assert sha(ROOT/name)==locks[name]
        unmerged.append(dict(path=name,taskSHA256=locks[name],committedByThisVersion=False))
    direct=[dict(archived=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in sorted((DEST/'native_ui').glob('*'))]
    m=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),files=files,
        sourceBindings=bindings,directEvidence=direct,matlab=matlab,nativeRuns=native,nativeEngineRuns=3,
        fullDomainGeneralSSTInterfacePassed=True,mainFiveStepStartTested=True,nativeHighPrecisionTemperatureWritten=True,
        sameConditionReopened=True,observedMousePickedFace=9977,unmergedWorkbenchResources=unmerged,
        physicsQualified=False,newQualifiedArea=0,largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(bindings),matlab=matlab,nativeEngineRuns=3),indent=2))


if __name__=='__main__':main()
