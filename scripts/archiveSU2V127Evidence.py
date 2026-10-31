"""归档完整域同源数值对照、小型科研证据和准确源码；大型体场只留本地。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2CheckpointContinuation import ROOT,read,sha
BASE=ROOT/'outputs/v127_validation';DEST=ROOT/'docs/validation/v127'


def main():
    assert not (DEST/'delivery_manifest.json').exists();matlab=[];native=[]
    for p in sorted(BASE.glob('sst_numerical_workbench*_process.json')):
        r=read(p);assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
        matlab.append(dict(name=p.name,exitCode=r['exitCode'],sourceStable=True))
    assert matlab and matlab[-1]['exitCode']==0
    for name in ['control','cap05']:
        r=read(BASE/name/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        native.append(dict(variant=name,exitCode=0,elapsedSeconds=r['elapsedSeconds']))
    files=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or p.stat().st_size>4*1024**2:continue
        if p.suffix not in ['.json','.log','.cfg'] and p.name not in ['history.csv','section.csv']:continue
        q=DEST/'evidence'/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);assert not q.exists();shutil.copyfile(p,q)
        files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    selected=set('prepareSU2SSTCheckpoint.py prepareSU2SSTNumericalStudy.py su2ThermalRuntime.py su2SSTInterface.py prepareSU2CheckpointContinuation.py runEmbeddedSU2.py runMatlabResearch.py auditSU2SSTState.py testSU2SSTNumericalWorkbench.m auditSU2ResultConvergence.m readSU2SSTStateAudit.m readSU2VerifiedSSTStudy.m describeSU2Convergence.m describeSU2Runtime.m readSU2SurfaceResult.m mapSU2WallResult.m buildSU2WallProjection.m readSU2HybridVolume.m attachSU2ThermalBoundary.m createSU2ResearchView.m createAerothermalGeometryView.m updateAerothermalGeometryView.m rocketAeroHeatWorkbench.m'.split())
    available={sha(p):p for d in [DEST/'initial_source_versions',DEST/'failure_versions'] for p in d.glob('*') if p.is_file()}
    bindings=[]
    def save(p,digest):
        source=p if p.is_file() and sha(p)==digest else available.get(digest);assert source is not None,(str(p),digest)
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
            for x in data:collect(x)
    for p in BASE.rglob('*.json'):
        if p.stat().st_size<4*1024**2:collect(read(p))
    for name in ['scripts/archiveSU2V127Evidence.py','scripts/deduplicateSU2V127Copies.py','scripts/compareSU2SSTNumericalStudy.py','scripts/diagnoseSU2SSTStateExtrema.py',
                 'scripts/plotSU2SSTComparison.py','tests/testSU2SSTNumericalContract.py','tests/testSU2SSTStudyManifest.py','tests/verifySU2SSTNumericalExports.py']:
        p=ROOT/name;save(p,sha(p))
    for rel in ['solvers/CSolver.cpp','output/CFlowCompOutput.cpp','output/CFlowOutput.cpp']:
        p=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src'/rel;save(p,sha(p))
    locks=read(BASE/matlab[-1]['name'])['sourceSHA256'];unmerged=[]
    for name in ['web/viewer.js','web/section.js']:
        assert sha(ROOT/name)==locks[name];unmerged.append(dict(path=name,taskSHA256=locks[name],committedByThisVersion=False))
    direct=[dict(archived=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in [DEST/'comparison.png',DEST/'plot_provenance.json']]
    comparison=read(BASE/'comparison.json')
    m=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),files=files,
        sourceBindings=bindings,directEvidence=direct,matlab=matlab,nativeRuns=native,nativeEngineRuns=2,
        allScreeningGatesPassed=comparison['allScreeningGatesPassed'],mainWorkbenchReadOnlyTested=True,
        newMainSolveLaunchTested=False,newNativeScreenshots=False,unmergedWorkbenchResources=unmerged,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False,largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(bindings),matlab=matlab,nativeEngineRuns=2),indent=2))


if __name__=='__main__':main()
