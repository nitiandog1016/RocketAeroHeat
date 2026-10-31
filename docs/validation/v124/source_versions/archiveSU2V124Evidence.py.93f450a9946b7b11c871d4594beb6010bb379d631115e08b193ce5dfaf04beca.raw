"""保存v124原终态、后验失败和准确源码，不用原生成功覆盖失败，不复制大场。"""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v124_validation';DEST=ROOT/'docs/validation/v124'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def put(name,raw):
    p=DEST/'source_versions'/(name+'.'+hashlib.sha256(raw).hexdigest()+'.raw');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==raw
    else:p.write_bytes(raw)
    return p

def main():
    assert not (DEST/'delivery_manifest.json').exists()
    case=Path(read(BASE/'workbench03/engine_folder.json')['folder']);case.relative_to(ROOT/'outputs/embedded_workbench')
    matlab=[]
    for name,code in [('sst_plate01',0),('sst_workbench01',1),('sst_workbench02',1),('sst_workbench03',1),('sst_controls04',0)]:
        r=read(BASE/(name+'_process.json'));assert r['status']=='terminal' and r['exitCode']==code and all(r['sourceStableAfter'].values())
        matlab.append(dict(name=name,exitCode=code,sourceStable=True,sourceReceiptSHA256=sha(BASE/(name+'_process.json'))))
    for folder in [BASE/'plate01/plate',case]:
        r=read(folder/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
    assert read(case/'workflow/pipeline_process.json')['exitCode']==1
    roots=[(BASE,DEST/'evidence'),(case,DEST/'evidence'/'actual')];files=[]
    for source,target in roots:
        for p in sorted(source.rglob('*')):
            if not p.is_file() or p.stat().st_size>4*1024**2:continue
            if p.suffix not in ['.json','.log','.cfg'] and p.name!='history.csv':continue
            q=target/p.relative_to(source);q.parent.mkdir(parents=True,exist_ok=True)
            if q.exists():assert sha(q)==sha(p)
            else:shutil.copyfile(p,q)
            files.append(dict(source=p.relative_to(ROOT).as_posix(),archived=q.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    selected={'rocketAeroHeatWorkbench.m','embeddedSU2Job.m','aeroHeatWorkflowState.m','writeSU2RocketCase.m','buildSU2FlatPlateCase.m',
        'readSU2HybridVolume.m','readSU2SurfaceResult.m','mapSU2WallResult.m','buildSU2WallProjection.m','resolveSU2ThermalBoundary.m',
        'attachSU2ThermalBoundary.m','readSU2CoolPropStateAudit.m','auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m',
        'describeSU2Runtime.m','describeSU2Convergence.m','estimateSU2CheckpointStorage.m','createSU2ResearchView.m',
        'runSU2WorkbenchPipeline.py','runEmbeddedSU2.py','su2ThermalRuntime.py','su2SSTInterface.py','prepareSU2SSTCase.py',
        'prepareSU2CheckpointContinuation.py','auditSU2SSTField.py','auditSU2SSTState.py','prepareSU2CoolPropCase.py','auditSU2CoolPropState.py',
        'testSU2SSTInterface.m','testSU2SSTWorkbench.m','testSU2SSTContracts.py','diagnoseSU2SSTPhase.py','testSU2WorkbenchPipeline.py',
        'runMatlabResearch.py','CFluidIteration.cpp','CMultiGridIntegration.cpp'}
    available={sha(p):p for directory in [DEST/'first_audit_failure',ROOT/'docs/validation/v123/source_versions'] for p in directory.glob('*') if p.is_file()}
    # 恢复首次准备脚本（后来只加源历史/谱系核验）；仅以真实任务SHA接受恢复的准确版本。
    p=ROOT/'scripts/prepareSU2SSTCase.py';text=p.read_text(encoding='utf-8')
    for line in ["        hist=np.atleast_1d(np.genfromtxt(str(source/'history.csv'),names=True,delimiter=','))\n",
                 "        assert np.array_equal(hist['Inner_Iter'],np.arange(len(hist))) and step==int(hist['Inner_Iter'][-1])\n",
                 "            assert sha(source/'continuation.json')==r['inputSHA256'][str(source/'continuation.json')]\n",
                 "            assert sha(source/'initial_flow.transfer.json')==r['inputSHA256'][str(source/'initial_flow.transfer.json')]\n"]:
        assert line in text;text=text.replace(line,'',1)
    expected=read(BASE/'sst_plate01_process.json')['sourceSHA256']['scripts/prepareSU2SSTCase.py']
    raw=next(x for x in [text.encode(),text.replace('\n','\r\n').encode()] if hashlib.sha256(x).hexdigest()==expected)
    available[expected]=put(p.name,raw)
    # 平板最终只读审核在移至scripts之前，原表头修正可按已存SHA准确恢复。
    p=DEST/'first_audit_failure/auditSU2SSTField_header_failure.py';text=p.read_text(encoding='utf-8').replace('csv.reader(h)','csv.reader(h,skipinitialspace=True)')
    expected=read(BASE/'plate01/independent_audit.json')['sourceSHA256']['tests\\auditSU2SSTField.py']
    raw=next(x for x in [text.encode(),text.replace('\n','\r\n').encode()] if hashlib.sha256(x).hexdigest()==expected)
    available[expected]=put('auditSU2SSTField.py',raw)
    bindings=[]
    def collect(obj):
        if isinstance(obj,dict):
            for key,value in obj.items():
                if isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value):
                    p=Path(key);p=p if p.is_absolute() else ROOT/p
                    if p.name in selected:
                        source=p if p.exists() and sha(p)==value else available.get(value)
                        assert source is not None,(str(p),value)
                        q=put(p.name,source.read_bytes());bindings.append(dict(path=str(p),sha256=value,archived=q.relative_to(ROOT).as_posix()))
                else:collect(value)
        elif isinstance(obj,list):
            for item in obj:collect(item)
    for root,_ in roots:
        for p in root.rglob('*.json'):
            if p.stat().st_size<4*1024**2:collect(read(p))
    for name in ['tests/testSU2WorkbenchPipeline.py','tests/testSU2SSTContracts.py','scripts/archiveSU2V124Evidence.py']:
        p=ROOT/name;put(p.name,p.read_bytes())
    result=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
        files=files,sourceBindings=bindings,matlab=matlab,nativeEngineRuns=2,nativeEngineExitCodes=[0,0],
        fullDomainPipelineExitCode=1,fullDomainInterfacePassed=False,fullDomainResultAccepted=False,
        physicsQualified=False,newQualifiedArea=0,largeResultsBackedUpToGitHub=False,
        syntheticSupervision=read(BASE/'supervision_contracts/report.json'))
    (DEST/'delivery_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(list((DEST/'source_versions').glob('*.raw'))),matlab=matlab,fullDomainAccepted=False),indent=2))

if __name__=='__main__':main()
