"""归档v132局部离散、官方梯度对照及工作台的小证据，完整大场仅本地。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2SSTGradientStudy import ROOT,BASE,read,sha,write
DEST=ROOT/'docs/validation/v132'

def main():
    assert not (DEST/'delivery_manifest.json').exists()
    matlab=read(BASE/'sst_gradient_workbench01_process.json')
    assert matlab['status']=='terminal' and matlab['exitCode']==0 and all(matlab['sourceStableAfter'].values())
    runs=[]
    for name in ['wls','gg']:
        r=read(BASE/name/'engine_process.json')
        assert r['status']=='terminal' and all(r['inputsStable'].values())
        runs.append(dict(variant=name,exitCode=r['exitCode'],elapsedSeconds=r['elapsedSeconds']))
    files=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or p.stat().st_size>2*1024**2:continue
        if p.suffix not in ['.json','.log','.cfg','.png'] and p.name not in ['history.csv','section.csv']:continue
        q=DEST/'evidence'/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True)
        assert not q.exists();shutil.copyfile(str(p),str(q))
        files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),
                          bytes=p.stat().st_size,sha256=sha(p)))
    names=['scripts/reviewSU2SSTLocalDiffusion.py','scripts/prepareSU2SSTGradientStudy.py',
        'scripts/compareSU2SSTGradientStudy.py','scripts/su2ThermalRuntime.py',
        'scripts/deduplicateSU2V132Copies.py','scripts/archiveSU2V132Evidence.py',
        'scripts/auditSU2SSTState.py','scripts/su2SSTInterface.py','scripts/runEmbeddedSU2.py',
        'scripts/prepareSU2SSTCheckpoint.py','scripts/runMatlabResearch.py',
        'tests/testSU2SSTLocalDiffusion.py','tests/testSU2SSTGradientStudy.py',
        'tests/testSU2SSTGradientWorkbench.m','tests/verifySU2SSTGradientExports.py',
        'readSU2VerifiedSSTStudy.m','auditSU2ResultConvergence.m','describeSU2Runtime.m',
        'readSU2SSTStateAudit.m','describeSU2Convergence.m','rocketAeroHeatWorkbench.m',
        'createSU2ResearchView.m','createAerothermalGeometryView.m','updateAerothermalGeometryView.m']
    bindings=[]
    for rel in names:
        p=ROOT/rel;h=sha(p)
        if rel in matlab['sourceSHA256']:assert h==matlab['sourceSHA256'][rel]
        q=DEST/'source_versions'/(p.name+'.'+h+'.raw');q.parent.mkdir(exist_ok=True)
        assert not q.exists();shutil.copyfile(str(p),str(q))
        bindings.append(dict(path=rel,sha256=h,archived=q.relative_to(ROOT).as_posix()))
    for path,h in read(BASE/'local_diffusion.json')['sourceSHA256'].items():
        p=Path(path)
        if 'su2-thermal-build/source' not in p.as_posix():continue
        assert sha(p)==h;q=DEST/'upstream'/(p.name+'.v86.raw');assert not q.exists()
        shutil.copyfile(str(p),str(q));bindings.append(dict(path=path,sha256=h,archived=q.relative_to(ROOT).as_posix()))
    unmerged=[dict(path=name,taskSHA256=matlab['sourceSHA256'][name],committedByThisVersion=False)
              for name in ['web/viewer.js','web/section.js']]
    for item in unmerged:assert sha(ROOT/item['path'])==item['taskSHA256']
    evidence=sum(p.stat().st_size for p in DEST.rglob('*') if p.is_file())
    logical=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert evidence<10*1024**2 and logical<1500*1024**2
    report=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
        files=files,sourceBindings=bindings,nativeRuns=runs,matlabExitCode=matlab['exitCode'],
        evidenceBytesBeforeManifest=evidence,pathLengthBytes=logical,pathLengthIsNotDiskUsage=True,
        unmergedWorkbenchResources=unmerged,largeResultsBackedUpToGitHub=False,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False)
    write(DEST/'delivery_manifest.json',report)
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(bindings),evidenceBytes=evidence,pathLengthBytes=logical)))

if __name__=='__main__':main()
