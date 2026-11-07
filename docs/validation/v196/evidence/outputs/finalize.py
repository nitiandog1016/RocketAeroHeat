"""收尾仅核对本轮小证据/保护资产；不清理独有数据或复制工程。"""
from pathlib import Path
import hashlib,json,shutil,sys,time,py_compile
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'outputs/v196_validation'
sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import sha
from hybridRocketEuler import read_wall
from geometricRocketEuler import check_admission
from prepareHybridRocketEulerNative import save

start=time.monotonic()
owned=['scripts/hybridRocketEuler.py','scripts/checkEilmerCellTimeScale.py','scripts/observeEilmerCellABI.py',
       'scripts/buildEilmerGeometryRuntime.py','scripts/geometricRocketEuler.py','scripts/auditGeometricRocketEuler.py',
       'scripts/plotGeometricRocketFailure.py','tests/testEilmerGeometricLength.py','tests/testGeometricRocketAdmission.py',
       'tests/testGeometricRocketFailure.py','patches/eilmer/geometric-length-v1.patch','docs/PROJECT_GOALS.md',
       'docs/HYBRID_GEOMETRIC_TIME_SCALE_V196.md','AGENTS.md']
for name in owned:
    if name.endswith('.py'):py_compile.compile(str(ROOT/name),doraise=True)
# 最终共享读取器的旧/新实际原面初场相容；只读，没有CFD或假造外缘。
old=ROOT/'outputs/v195_validation';a=read_wall(old/'rocket01',0,old/'preparation01/report.json')
b=read_wall(BASE/'rocket01',0,BASE/'preparation01/report.json')
assert a['fields']==b['fields'] and a['faceIDs']==b['faceIDs']
assert all(x is None for x in b['wallHeatFluxWm2']) and not b['engineeringHeatReady']
check_admission(json.loads((BASE/'admission.json').read_text()))
protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
stable={name:sha(ROOT/Path(name.replace('\\','/')))==d for name,d in protected.items()}
assert len(stable)==29 and all(stable.values())
processes=[]
for p in BASE.rglob('report.json'):
    r=json.loads(p.read_text())
    if r.get('schema')=='RAH.EilmerProcess.v1':
        processes.append(dict(path=str(p.relative_to(BASE)).replace('\\','/'),elapsedSeconds=r['elapsedSeconds'],
                              exitCode=r['nativeExitCode'],passed=r['passed'],inputsStable=r['inputsStable']))
save(BASE/'final_receipt.json',dict(passed=True,protectedFilesStable=stable,originalFacesCompared=19752,
    exactInitialFieldCompatibility=True,newWallHeatAreaM2=0,nativeProcesses=processes,
    sumNativeSupervisorSeconds=sum(r['elapsedSeconds'] for r in processes),
    timingScope='sum of distinct recorded native/ABI/build/prep supervisor processes; excludes Python wrappers, installation, tool parents and development',
    elapsedSeconds=time.monotonic()-start,physicalQualification=False,
    sourceSHA256={name:sha(ROOT/name) for name in owned}))
save(BASE/'manual_tool_failures.json',dict(scope='Observed tool failures, not reconstructed process receipts',failures=[
    dict(stage='initial-symbol-observation',reason='No DWARF local symbol table; original native report preserved'),
    dict(stage='inline-offset-query',reason='PowerShell inline quoting caused Python SyntaxError before execution; no native process'),
    dict(stage='offsets01',reason='import path missed src; actual compiler failure/source preserved'),
    dict(stage='first-small-tetra-audit',reason='three mutable runtime outputs incorrectly locked as inputs; native normal, read-only recovery'),
    dict(stage='first-rocket-preparation',reason='Windows/WSL source-key KeyError before case creation; exact source preserved',toolParentSeconds=2.196914)
],dependencyInstall=dict(toolParentSeconds=6.2232852,packages=7,downloadApproxMB=4.646,installedApproxMB=12.4)))
dest=ROOT/'docs/validation/v196/evidence';dest.mkdir(exist_ok=False)
copied={}
def copy(src,relative):
    target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,target)
    assert sha(src)==sha(target)
    copied[str(relative).replace('\\','/')]=dict(sha256=sha(src),bytes=src.stat().st_size)
for p in BASE.rglob('*'):
    if not p.is_file() or '__pycache__' in p.parts:continue
    rel=p.relative_to(BASE)
    if p.suffix in {'.json','.txt','.py','.d','.gdb','.csv'} or p.name=='rejected_development_final.png':
        copy(p,Path('outputs')/rel)
for name in owned:copy(ROOT/name,Path('final_source')/name)
# 包的清单可公开；运行库和大场不上传。
copy(ROOT/'tools/runtime/eilmer-geometric-length-v1/runtime.json',Path('runtime.json'))
save(ROOT/'docs/validation/v196/evidence_manifest.json',dict(files=copied,count=len(copied),
    bytes=sum(x['bytes'] for x in copied.values()),largeNativeFieldsPublished=False,
    protectedStableCount=29,uniqueFilesDeleted=0,physicalQualification=False))
print(json.dumps(dict(passed=True,evidenceCount=len(copied),evidenceBytes=sum(x['bytes'] for x in copied.values()),
    nativeSupervisorSumSeconds=sum(r['elapsedSeconds'] for r in processes),elapsedSeconds=time.monotonic()-start)))
