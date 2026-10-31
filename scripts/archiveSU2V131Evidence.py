"""归档v131有界输运的小型实际证据；大场/运行包留本地，未合并资源只记散列。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2SSTBoundedStudy import ROOT,BASE,sha,read,write
DEST=ROOT/'docs/validation/v131'
def main():
    assert not (DEST/'delivery_manifest.json').exists()
    matlab=[]
    for p in sorted(BASE.glob('sst_bounded_workbench*_process.json')):
        r=read(p);assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
        matlab.append(dict(path=str(p),exitCode=r['exitCode']))
    runs=[]
    for p in BASE.rglob('engine_process.json'):
        r=read(p);assert r['status']=='terminal' and all(r['inputsStable'].values())
        runs.append(dict(folder=str(p.parent),exitCode=r['exitCode'],elapsedSeconds=r['elapsedSeconds']))
    files=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or p.stat().st_size>3*1024**2:continue
        permitted=p.suffix in ['.json','.log','.cfg','.used','.original','.patched','.cpp','.py','.ninja']
        permitted=permitted or p.name in ['history.csv','derivatives.csv','section.csv']
        permitted=permitted or (p.name=='mesh.su2' and p.parent.name.startswith('small_'))
        if not permitted:continue
        q=DEST/'evidence'/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True)
        assert not q.exists();shutil.copyfile(str(p),str(q))
        files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    selected=set('buildSU2BoundedRuntime.py buildSU2ThermalRuntime.py verifySU2BoundedChainRule.py prepareSU2SSTBoundedStudy.py checkSU2BoundedPackage.py compareSU2SSTBoundedStudy.py reviewSU2BoundedBuild.py su2ThermalRuntime.py su2SSTInterface.py runEmbeddedSU2.py auditSU2SSTState.py testSU2SSTBoundedWorkbench.m readSU2VerifiedSSTStudy.m auditSU2ResultConvergence.m describeSU2Convergence.m describeSU2Runtime.m readSU2SSTStateAudit.m readSU2SurfaceResult.m mapSU2WallResult.m buildSU2WallProjection.m attachSU2ThermalBoundary.m readSU2HybridVolume.m createSU2ResearchView.m createAerothermalGeometryView.m updateAerothermalGeometryView.m rocketAeroHeatWorkbench.m runMatlabResearch.py'.split())
    available={sha(p):p for base in [ROOT/'docs/validation/v130/source_versions',DEST/'failure_versions'] for p in base.glob('*') if p.is_file()}
    available.update({sha(p):p for p in BASE.rglob('*.py') if p.is_file()})
    bindings=[]
    def save(p,h):
        source=p if p.is_file() and sha(p)==h else available.get(h)
        assert source is not None,(str(p),h)
        q=DEST/'source_versions'/(p.name+'.'+h+'.raw');q.parent.mkdir(exist_ok=True)
        if q.exists():assert sha(q)==h
        else:shutil.copyfile(str(source),str(q))
        row=dict(path=str(p),sha256=h,archived=q.relative_to(ROOT).as_posix())
        if row not in bindings:bindings.append(row)
    def collect(v):
        if isinstance(v,dict):
            for key,item in v.items():
                if isinstance(item,str) and len(item)==64 and all(c in '0123456789abcdef' for c in item):
                    p=Path(key);p=p if p.is_absolute() else ROOT/p
                    if p.name in selected:save(p,item)
                else:collect(item)
        elif isinstance(v,list):
            for item in v:collect(item)
    for p in BASE.rglob('*.json'):
        if p.stat().st_size<3*1024**2:collect(read(p))
    for name in ['scripts/archiveSU2V131Evidence.py','scripts/finalizeSU2BoundedRuntime.py','scripts/deduplicateSU2V131Copies.py',
                 'tests/testSU2BoundedContract.py','tests/verifySU2SSTBoundedExports.py']:
        p=ROOT/name;save(p,sha(p))
    unmerged=[]
    if matlab:
        lock=read(Path(matlab[-1]['path']))['sourceSHA256']
        for name in ['web/viewer.js','web/section.js']:
            assert sha(ROOT/name)==lock[name]
            unmerged.append(dict(path=name,taskSHA256=lock[name],committedByThisVersion=False))
    logical=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    evidence=sum(p.stat().st_size for p in DEST.rglob('*') if p.is_file())
    assert logical<=900*1024**2 and evidence<=12*1024**2
    report=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
        files=files,sourceBindings=bindings,nativeRuns=runs,matlab=matlab,
        storage=dict(pathLengthBytes=logical,evidenceBytesBeforeManifest=evidence,notDiskAllocation=True,
                     buildCacheDelta=read(BASE/'build_final/cache_delta.json')),
        unmergedWorkbenchResources=unmerged,largeResultsBackedUpToGitHub=False,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False)
    write(DEST/'delivery_manifest.json',report)
    print(json.dumps(dict(files=len(files),sourceVersions=len(bindings),evidenceBytes=evidence,pathLengthBytes=logical),indent=2))
if __name__=='__main__':main()
