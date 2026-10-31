"""归档空间诊断、初场对照及主工作台的小型证据；大场保留本地。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2CheckpointContinuation import ROOT,read,sha
BASE=ROOT/'outputs/v128_validation';DEST=ROOT/'docs/validation/v128'


def main():
    assert not (DEST/'delivery_manifest.json').exists();matlab=[];native=[]
    for p in sorted(BASE.glob('sst_pocket_workbench*_process.json')):
        r=read(p);assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
        matlab.append(dict(name=p.name,exitCode=r['exitCode'],sourceStable=True))
    assert matlab and matlab[-1]['exitCode']==0
    for name in ['spatial01','pocket_control','pocket_repaired_roundoff']:
        r=read(BASE/name/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        native.append(dict(variant=name,exitCode=0,elapsedSeconds=r['elapsedSeconds']))
    files=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or p.stat().st_size>4*1024**2:continue
        if p.suffix not in ['.json','.log','.cfg'] and p.name not in ['history.csv','section.csv']:continue
        q=DEST/'evidence'/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);assert not q.exists();shutil.copyfile(p,q)
        files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    selected=set('prepareSU2SSTCheckpoint.py prepareSU2SSTPocketStudy.py prepareSU2SSTSpatialDiagnostic.py su2ThermalRuntime.py su2SSTInterface.py prepareSU2CheckpointContinuation.py runEmbeddedSU2.py runMatlabResearch.py auditSU2SSTState.py testSU2SSTPocketWorkbench.m auditSU2ResultConvergence.m readSU2SSTStateAudit.m readSU2VerifiedSSTStudy.m describeSU2Convergence.m describeSU2Runtime.m readSU2SurfaceResult.m mapSU2WallResult.m buildSU2WallProjection.m readSU2HybridVolume.m attachSU2ThermalBoundary.m createSU2ResearchView.m createAerothermalGeometryView.m updateAerothermalGeometryView.m rocketAeroHeatWorkbench.m diagnoseSU2SSTSpatial.py'.split())
    available={sha(p):p for d in [DEST/'initial_source_versions',DEST/'failure_versions',ROOT/'docs/validation/v127/source_versions'] for p in d.glob('*') if p.is_file()}
    bindings=[]
    def save(p,h):
        src=p if p.is_file() and sha(p)==h else available.get(h);assert src is not None,(str(p),h)
        q=DEST/'source_versions'/(p.name+'.'+h+'.raw');q.parent.mkdir(parents=True,exist_ok=True)
        if q.exists():assert sha(q)==h
        else:shutil.copyfile(src,q)
        row=dict(path=str(p),sha256=h,archived=q.relative_to(ROOT).as_posix())
        if row not in bindings:bindings.append(row)
    def collect(obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                if isinstance(v,str) and len(v)==64 and all(c in '0123456789abcdef' for c in v):
                    p=Path(k);p=p if p.is_absolute() else ROOT/p
                    if p.name in selected:save(p,v)
                else:collect(v)
        elif isinstance(obj,list):
            for x in obj:collect(x)
    for p in BASE.rglob('*.json'):
        if p.stat().st_size<4*1024**2:collect(read(p))
    for name in ['scripts/archiveSU2V128Evidence.py','scripts/deduplicateSU2V128Copies.py','scripts/compareSU2SSTPocketStudy.py',
                 'scripts/plotSU2SSTPocketStudy.py','tests/testSU2SSTSpatialDiagnostic.py','tests/testSU2SSTPocketInitialization.py','tests/verifySU2SSTPocketExports.py']:
        p=ROOT/name;save(p,sha(p))
    nativeRoot=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD'
    for rel in ['src/solvers/CSolver.cpp','src/solvers/CTurbSolver.cpp','src/solvers/CTurbSSTSolver.cpp',
        'src/solvers/CEulerSolver.cpp','src/output/CFlowCompOutput.cpp','src/output/CFlowOutput.cpp','src/iteration/CFluidIteration.cpp',
        'include/solvers/CFVMFlowSolverBase.hpp','include/solvers/CScalarSolver.inl']:
        p=nativeRoot/rel;save(p,sha(p))
    for rel in ['Common/src/CConfig.cpp','config_template.cfg']:
        p=nativeRoot.parent/rel;save(p,sha(p))
    locks=read(BASE/matlab[-1]['name'])['sourceSHA256'];unmerged=[]
    for name in ['web/viewer.js','web/section.js']:
        assert sha(ROOT/name)==locks[name];unmerged.append(dict(path=name,taskSHA256=locks[name],committedByThisVersion=False))
    direct=[dict(archived=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in [DEST/'pocket_comparison.png',DEST/'plot_provenance.json']]
    result=read(BASE/'pocket_comparison.json')
    m=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),files=files,
        sourceBindings=bindings,directEvidence=direct,matlab=matlab,nativeRuns=native,nativeEngineRuns=3,
        unusedPreparedCandidate='pocket_repaired; 111-node exact-equality preflight preparation never run',
        allScreeningGatesPassed=result['allScreeningGatesPassed'],mainWorkbenchReadOnlyTested=True,newMainSolveLaunchTested=False,
        newNativeScreenshots=False,unmergedWorkbenchResources=unmerged,physicsQualified=False,newQualifiedArea=0,
        defaultChanged=False,largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(bindings),matlab=matlab,nativeEngineRuns=3),indent=2))


if __name__=='__main__':main()
