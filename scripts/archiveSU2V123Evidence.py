"""绑定v123真实主UI、阶段终态、独立复核及准确源版本；不复制大场或运行包。"""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v123_validation';DEST=ROOT/'docs/validation/v123'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def putraw(name,raw):
    p=DEST/'source_versions'/(name+'.'+hashlib.sha256(raw).hexdigest()+'.raw');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==raw
    else:p.write_bytes(raw)
def main():
    audit=read(BASE/'independent_audit.json');assert audit['passed']
    # 旧状态模块在首次测试时未修改；精确移除新分支，按原散列选择LF/CRLF。
    p=ROOT/'aeroHeatWorkflowState.m';old=read(BASE/'heos_ui01_process.json')['sourceSHA256'][p.name]
    text=p.read_text(encoding='utf-8');a=text.index("    case 'invalidateboundary'");b=text.index("    case 'setboundary'",a);text=text[:a]+text[b:]
    raw=next(raw for raw in [text.encode(),text.replace('\n','\r\n').encode()] if hashlib.sha256(raw).hexdigest()==old)
    putraw(p.name,raw)
    roots=[(BASE,DEST/'evidence')]+[(Path(p['folder']),DEST/'evidence'/'pipelines'/Path(p['folder']).name) for p in audit['pipelines']]
    files=[]
    for source,target in roots:
        for p in sorted(source.rglob('*')):
            if not p.is_file() or p.stat().st_size>4*1024**2:continue
            if p.suffix not in ['.json','.log','.cfg'] and p.name not in ['fields.csv','short_fields.csv','section.csv','history.csv']:continue
            q=target/p.relative_to(source);q.parent.mkdir(parents=True,exist_ok=True)
            if q.exists():assert sha(q)==sha(p),'归档不能覆盖不同原文'
            else:shutil.copyfile(p,q)
            files.append(dict(source=str(p.relative_to(ROOT)),archived=str(q.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size))
    selected={'aeroHeatWorkflowState.m','rocketAeroHeatWorkbench.m','embeddedSU2Job.m','readSU2VerifiedHEOSCase.m',
        'readSU2SurfaceResult.m','mapSU2WallResult.m','buildSU2WallProjection.m','writeSU2RocketCase.m',
        'readSU2HybridVolume.m','resolveSU2ThermalBoundary.m','attachSU2ThermalBoundary.m','readSU2CoolPropStateAudit.m',
        'auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m','describeSU2Runtime.m','describeSU2Convergence.m',
        'createSU2ResearchView.m','createAerothermalGeometryView.m','updateAerothermalGeometryView.m','estimateSU2CheckpointStorage.m',
        'geometry.html','geometry-su2.html','geometry.js','geometry.css','su2-section-view.js','su2-section.css','ale-section-core.js',
        'runSU2WorkbenchPipeline.py','prepareSU2CoolPropCase.py','auditSU2CoolPropState.py','runEmbeddedSU2.py',
        'su2ThermalRuntime.py','buildSU2ThermalRuntime.py','testSU2HEOSWorkbench.m','runMatlabResearch.py'}
    available={sha(p):p for p in (DEST/'source_versions').glob('*.raw')};bindings=[]
    def collect(obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                if isinstance(v,str) and len(v)==64 and all(c in '0123456789abcdef' for c in v):
                    p=Path(k);p=p if p.is_absolute() else ROOT/p
                    if p.name in selected and p.exists():
                        src=p if sha(p)==v else available.get(v);assert src is not None,(str(p),v)
                        putraw(p.name,src.read_bytes());bindings.append(dict(path=str(p),sha256=v))
                else:collect(v)
        elif isinstance(obj,list):
            for item in obj:collect(item)
    for source,_ in roots:
        for p in source.rglob('*.json'):
            if p.stat().st_size<4*1024**2:collect(read(p))
    # 没在执行时锁定的独立验证与维护程序也保存当前准确版本。
    for name in ['tests/auditSU2HEOSWorkbench.py','tests/testSU2WorkbenchPipeline.py','scripts/deduplicateSU2V123Copies.py']:
        p=ROOT/name;putraw(p.name,p.read_bytes())
    fixed=ROOT/'tools/runtime/su2-thermal-build/source';sources=[]
    for name,start,end in [('Common/include/option_structure.hpp',1176,1234),
        ('SU2_CFD/src/solvers/CTurbSSTSolver.cpp',94,151),('SU2_CFD/src/solvers/CTurbSSTSolver.cpp',507,532),
        ('SU2_CFD/src/solvers/CEulerSolver.cpp',1018,1028),('SU2_CFD/src/variables/CNSVariable.cpp',1,180),
        ('SU2_CFD/src/output/CFlowOutput.cpp',1360,1376),('config_template.cfg',585,598)]:
        p=fixed/name;lines=p.read_text(encoding='utf-8').splitlines()
        sources.append(dict(path=name,sha256=sha(p),startLine=start,endLine=min(end,len(lines)),excerpt='\n'.join(lines[start-1:end])))
    (DEST/'sst_source_review.json').write_text(json.dumps(dict(pinnedSU2Commit='7048fd35f56d4bf52e586429d6746a36a10255f8',
        implementationNotTested=True,sources=sources),indent=2)+'\n')
    result=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
        matlab=audit['matlab'],pipelines=audit['pipelines'],nativeEngineRuns=1,files=files,sourceBindings=bindings,
        syntheticSupervision=read(BASE/'supervision_contracts/report.json'),physicsQualified=False,newQualifiedArea=0,
        largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(evidenceFiles=len(files),sourceVersions=len(list((DEST/'source_versions').glob('*.raw'))),
                         matlab=result['matlab'],nativeEngineRuns=1),indent=2))
if __name__=='__main__':main()
