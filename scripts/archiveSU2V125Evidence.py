"""归档v125真实一次更新和只读查看证据；大场不复制，历史失败不改写。"""
from pathlib import Path
import json,hashlib,shutil,subprocess
from prepareSU2CheckpointContinuation import ROOT,read,sha
BASE=ROOT/'outputs/v125_validation';DEST=ROOT/'docs/validation/v125'


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    native=read(BASE/'phase01/engine_process.json');assert native['status']=='terminal' and native['exitCode']==0 and all(native['inputsStable'].values())
    matlab=[]
    for name,code in [('sst_phase_view01',1),('sst_phase_view02',0)]:
        r=read(BASE/(name+'_process.json'));assert r['status']=='terminal' and r['exitCode']==code and all(r['sourceStableAfter'].values())
        matlab.append(dict(name=name,exitCode=code,sourceStable=True))
    files=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or p.stat().st_size>4*1024**2:continue
        if p.suffix not in ['.json','.log','.cfg'] and p.name!='history.csv':continue
        q=DEST/'evidence'/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);assert not q.exists();shutil.copyfile(p,q)
        files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    selected=set('prepareSU2SSTPhaseProbe.py su2ThermalRuntime.py su2SSTInterface.py prepareSU2CheckpointContinuation.py prepareSU2SSTCase.py runEmbeddedSU2.py runMatlabResearch.py auditSU2SSTState.py auditSU2SSTPhaseProbe.py testSU2SSTPhaseView.m auditSU2ResultConvergence.m readSU2SSTStateAudit.m describeSU2Convergence.m describeSU2Runtime.m readSU2SurfaceResult.m mapSU2WallResult.m buildSU2WallProjection.m readSU2HybridVolume.m attachSU2ThermalBoundary.m createSU2ResearchView.m createAerothermalGeometryView.m updateAerothermalGeometryView.m rocketAeroHeatWorkbench.m'.split())
    available={sha(p):p for directory in [ROOT/'docs/validation/v124/source_versions',DEST/'first_view_failure'] for p in directory.glob('*') if p.is_file()}
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
            for x in data:collect(x)
    for p in BASE.rglob('*.json'):
        if p.stat().st_size<4*1024**2:collect(read(p))
    for name in ['scripts/archiveSU2V125Evidence.py','scripts/deduplicateSU2V125Copies.py','tests/testSU2SSTPhaseContract.py',
                 'tests/verifySU2SSTPhaseExports.py']:
        p=ROOT/name;save(p,sha(p))
    sourceRoot=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src'
    for rel in ['iteration/CFluidIteration.cpp','integration/CMultiGridIntegration.cpp','output/COutput.cpp']:
        p=sourceRoot/rel;save(p,sha(p))
    # 本轮不合并OpenCode源码；记录实际科研组件使用的资源身份，不能称远端完整工作区逐字相同。
    resources=['web/viewer.js','web/section.js'];locks=read(BASE/'sst_phase_view02_process.json')['sourceSHA256']
    unmerged=[]
    for name in resources:
        assert sha(ROOT/name)==locks[name]
        unmerged.append(dict(path=name,taskSHA256=locks[name],committedByThisVersion=False))
    m=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),files=files,sourceBindings=bindings,
        matlab=matlab,nativeEngineRuns=1,nativeEngineExitCodes=[0],nativeElapsedSeconds=native['elapsedSeconds'],
        phaseDiagnosticPassed=True,fullDomainGeneralSSTInterfacePassed=False,mainFiveStepStartTested=False,
        nativeHighPrecisionTemperatureWritten=False,oldV124FieldStillRejected=True,unmergedWorkbenchResources=unmerged,
        physicsQualified=False,newQualifiedArea=0,largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(bindings),matlab=matlab,nativeEngineRuns=1),indent=2))

if __name__=='__main__':main()
