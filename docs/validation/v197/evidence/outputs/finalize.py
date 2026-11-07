"""v197只收集本轮小证据和保护文件回执，不打包工程或删除独有数据。"""
from pathlib import Path
import json,shutil,sys,time,py_compile
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'outputs/v197_validation'
sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import sha
from prepareHybridRocketEulerNative import save
start=time.monotonic()
owned=['scripts/buildEilmerNKGeometryRuntime.py','scripts/hybridNKTetraControl.py',
       'scripts/auditHybridNKControl.py','scripts/hybridNKRocket.py','scripts/auditHybridNKRocket.py',
       'scripts/plotHybridNKDevelopment.py','tests/testHybridNKAdmission.py','tests/testHybridNKClosed.py',
       'docs/HYBRID_NEWTON_KRYLOV_V197.md','docs/PROJECT_GOALS.md','AGENTS.md']
for name in owned:
    if name.endswith('.py'):py_compile.compile(str(ROOT/name),doraise=True)
protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
stable={name:sha(ROOT/Path(name.replace('\\','/')))==digest for name,digest in protected.items()}
assert len(stable)==29 and all(stable.values())
processes=[]
for path in BASE.rglob('report.json'):
    r=json.loads(path.read_text(encoding='utf-8'))
    if r.get('schema')=='RAH.EilmerProcess.v1':
        assert r['status']=='process-completed' and not r['liveGroupAfterCleanup']
        processes.append(dict(path=path.relative_to(BASE).as_posix(),seconds=r['elapsedSeconds'],
                              nativeExitCode=r['nativeExitCode'],passed=r['passed'],inputsStable=r['inputsStable']))
save(BASE/'final_receipt.json',dict(passed=True,protectedFilesStable=stable,physicalQualification=False,
    originalWallFaces=19752,newSolvedHeatAreaM2=0,uniqueFilesDeleted=0,nativeProcesses=processes,
    nativeSupervisorSeconds=sum(x['seconds'] for x in processes),
    timingScope='distinct build/preparation/CFD supervisor processes; wrappers, development and tool coordination excluded',
    localVersionBytesBeforeEvidence=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),
    elapsedSeconds=time.monotonic()-start,sourceSHA256={name:sha(ROOT/name) for name in owned}))
save(BASE/'manual_tool_failures.json',dict(scope='observed read-only tool failures; no fabricated native receipt',failures=[
    dict(stage='source-preview',reason='Windows default GBK could not decode source; explicit UTF-8 recovery'),
    dict(stage='path-discovery',reason='guessed legacy/source paths absent; located actual paths without CFD rerun'),
    dict(stage='NK-metadata-discovery',reason='explicit snapshot-times-metadata absent; native NK uses restart-metadata, verified against exact source and diagnostics')]))
dest=ROOT/'docs/validation/v197/evidence';dest.mkdir(exist_ok=False);copied={}
def copy(src,relative):
    dst=dest/relative;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    assert sha(src)==sha(dst)
    copied[relative.as_posix()]=dict(sha256=sha(src),bytes=src.stat().st_size)
small_native={'config','nk-config','commands-to-steady-mode','fluid.metadata','restart-metadata','reference-residuals','nk-diagnostics'}
for path in BASE.rglob('*'):
    if not path.is_file() or '__pycache__' in path.parts or path.name=='surface.json':continue
    rel=path.relative_to(BASE)
    if 'lmrsim' in rel.parts and path.name not in small_native:continue
    if path.suffix in {'.json','.txt','.py','.d','.lua','.gas','.csv','.md'} or path.name in small_native or path.name=='nk_development.png':
        copy(path,Path('outputs')/rel)
for name in owned:copy(ROOT/name,Path('final_source')/name)
copy(ROOT/'tools/runtime/eilmer-nk-geometric-v1/runtime.json',Path('runtime.json'))
save(ROOT/'docs/validation/v197/evidence_manifest.json',dict(files=copied,count=len(copied),
    bytes=sum(x['bytes'] for x in copied.values()),largeNativeFieldsPublished=False,
    protectedStableCount=29,uniqueFilesDeleted=0,physicalQualification=False))
print(json.dumps(dict(passed=True,evidenceCount=len(copied),evidenceBytes=sum(x['bytes'] for x in copied.values()),
    nativeSupervisorSeconds=sum(x['seconds'] for x in processes),elapsedSeconds=time.monotonic()-start)))
