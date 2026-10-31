"""归档本轮小型证据与实际源码；大体场、运行时与失败二进制留本地。"""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v122_validation';DEST=ROOT/'docs/validation/v122'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def putraw(name,raw):
    p=DEST/'source_versions'/(name+'.'+hashlib.sha256(raw).hexdigest()+'.raw');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==raw
    else:p.write_bytes(raw)
def recover_previous_sources():
    # 准确回退两次CLI参数扩展，仅用于保留已执行旧版本；散列必须匹配原回执。
    p=ROOT/'tests/diagnoseCoolPropActualException.py';t=p.read_text(encoding='utf-8')
    t=t.replace('import hashlib,json,os,sys,argparse','import hashlib,json,os,sys')
    t=t.replace("    parser=argparse.ArgumentParser();parser.add_argument('--guarded',action='store_true');args=parser.parse_args()\n    old=ROOT/'outputs/v122_validation'/('actual02' if args.guarded else 'actual01')\n    out=ROOT/'outputs/v122_validation'/('exception_diagnostic03' if args.guarded else 'exception_diagnostic02')", "    old=ROOT/'outputs/v122_validation/actual01';out=ROOT/'outputs/v122_validation/exception_diagnostic02'")
    t=t.replace("    exe=ROOT/'tools/runtime/su2'/('8.6.0-rah-coolprop-recon-v1' if args.guarded else '8.6.0-rah-coolprop-v1')/'bin/SU2_CFD.exe'\n    gdb=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/gdb.exe'", "    exe=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-v1/bin/SU2_CFD.exe';gdb=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/gdb.exe'")
    expected=json.loads((OUT/'exception_diagnostic02/input_hashes.json').read_text())[str(p)]
    raw=t.encode();assert hashlib.sha256(raw).hexdigest()==expected;putraw(p.name,raw)
    p=ROOT/'tests/auditCoolPropGuardEquivalence.py';t=p.read_text(encoding='utf-8').replace('import hashlib,json,argparse','import hashlib,json')
    start=t.index('    parser=argparse.ArgumentParser()');end=t.index("    r=json.loads",start)
    t=t[:start]+"    old=ROOT/'outputs/v121_validation/thermal20';new=ROOT/'outputs/v122_validation/guard_control'\n"+t[end:]
    expected=json.loads((OUT/'guard_control/equivalence.json').read_text())['sources'][str(p)]
    raw=t.encode();assert hashlib.sha256(raw).hexdigest()==expected;putraw(p.name,raw)
def main():
    for p in OUT.rglob('engine_process.json'):
        r=json.loads(p.read_text());assert r['status']=='terminal','实际引擎仍运行，不能归档瞬时回执'
    recover_previous_sources();index=[]
    for p in sorted(OUT.rglob('*')):
        if not p.is_file():continue
        if p.suffix.lower() not in ['.json','.log','.gdb','.diff','.original','.used','.cpp','.cfg'] and p.name not in ['history.csv','section.csv','fields.csv']:continue
        if p.stat().st_size>6*1024**2:continue
        q=DEST/'evidence'/p.relative_to(OUT);q.parent.mkdir(parents=True,exist_ok=True)
        if q.exists():assert sha(q)==sha(p),'不可变归档冲突 '+str(q)
        else:shutil.copyfile(str(p),str(q))
        index.append(dict(source=str(p.relative_to(ROOT)),archived=str(q.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size))
    # 运行中/待写文件不能作为交付终态。
    runs={}
    for p in OUT.rglob('engine_process.json'):
        r=json.loads(p.read_text());assert r['status']=='terminal' and all(r['inputsStable'].values())
        runs[str(p.parent.relative_to(OUT))]={k:r[k] for k in ['status','exitCode','elapsedSeconds','runtimeFlavor']}
    # 仅冻结本轮实际绑定的脚本版本。历史源或大输入另由原SHA指向，不复制整工程。
    versions={sha(p):p for p in (DEST/'source_versions').glob('*.raw')}
    for p in ROOT.glob('docs/validation/v121/source_versions/*.raw'):versions.setdefault(sha(p),p)
    for p in OUT.rglob('*'):
        if p.is_file() and p.suffix in ['.used','.original']:versions.setdefault(sha(p),p)
    bindings=[]
    selected={'runEmbeddedSU2.py','su2ThermalRuntime.py','prepareSU2CoolPropCase.py','prepareActualCoolPropRun.py',
        'auditCoolPropFullFarfield.py','diagnoseCoolPropActualException.py','buildSU2ThermalRuntime.py',
        'buildSU2CoolPropReconstructionRuntime.py','prepareCoolPropGuardRuns.py','auditCoolPropGuardEquivalence.py',
        'fixCoolPropDiagnosticRuntime.py','linkCoolPropDiagnosticCorrection.py','auditCoolPropFormatProperties.py',
        'auditSU2CoolPropState.py','testCoolPropInitializationContracts.py','auditActualCoolPropField.py',
        'reviewActualCoolPropField.m','readSU2CoolPropStateAudit.m','auditSU2ResultConvergence.m','describeSU2Runtime.m',
        'describeSU2Convergence.m','writeSU2CoolPropInitialState.m','auditActualCoolPropExports.py',
        'CEulerSolver.cpp','CCoolProp.cpp','FlashRoutines.cpp','RAHCoolPropState.cpp'}
    def visit(value):
        if isinstance(value,dict):
            for k,v in value.items():
                if isinstance(v,str) and len(v)==64 and all(c in '0123456789abcdef' for c in v):
                    p=Path(k)
                    if not p.is_absolute():p=ROOT/p
                    if (p.name in selected or p.suffix=='.gdb') and p.exists():
                        if sha(p)==v:source=p
                        elif v in versions:source=versions[v]
                        else:continue
                        putraw(p.name,source.read_bytes());bindings.append(dict(path=str(p),sha256=v))
                else:visit(v)
        elif isinstance(value,list):
            for q in value:visit(q)
    for p in OUT.rglob('*.json'):
        if p.stat().st_size<6*1024**2:visit(json.loads(p.read_text(encoding='utf-8-sig')))
    # 每个实际引擎的宿主脚本都须可找到准确版本，不以当前文件冒充历史。
    allhash={sha(p) for p in (DEST/'source_versions').glob('*.raw')}
    for p in OUT.rglob('engine_process.json'):
        for name,h in json.loads(p.read_text())['inputSHA256'].items():
            f=Path(name)
            if f.suffix=='.py':assert h in allhash,(name,h)
    matlab={}
    for p in OUT.glob('review*_process.json'):
        r=json.loads(p.read_text());assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
        matlab[p.stem]=dict(status=r['status'],exitCode=r['exitCode'],sourceStable=True)
    result=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),engines=runs,matlab=matlab,
        files=index,sourceBindings=bindings,physicsQualified=False,newQualifiedArea=0,largeResultsBackedUpToGitHub=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(files=len(index),sourceVersions=len(allhash),engines=runs),indent=2))
if __name__=='__main__':main()
