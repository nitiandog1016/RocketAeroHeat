"""复核已完成构建：区分生成器依赖缓存与真实编译源码，不重新编译。"""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v121_validation'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source';BASE=ROOT/'tools/runtime/su2-coolprop-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not (OUT/'compiled.json').exists()
    process=json.loads((OUT/'compile_process.json').read_text());assert process['status']=='terminal' and process['exitCode']==0 and not process['timedOut']
    before=json.loads((OUT/'compile_source_tree.json').read_text());restore=json.loads((OUT/'source_restore.json').read_text());assert all(restore.values())
    restored={str(SOURCE/'meson.build'):OUT/'meson.build.used',str(SOURCE/'subprojects/CoolProp/dev/generate_headers.py'):OUT/'generate_headers.py.used'}
    caches={str(SOURCE/'subprojects/CoolProp/dev'/name) for name in ['.fluiddepcache','.incompdepcache']}
    cacheChanges=[];checked=0
    for path,h in before.items():
        p=Path(path);assert p.is_file()
        if path in restored:
            assert sha(restored[path])==h;assert sha(p)==sha(OUT/(p.name+'.original'));checked+=1
        elif path in caches:
            cacheChanges.append(dict(path=path,beforeSHA256=h,afterSHA256=sha(p),bytes=p.stat().st_size))
        else:assert sha(p)==h,path;checked+=1
    assert len(cacheChanges)==2 and all(c['beforeSHA256']!=c['afterSHA256'] for c in cacheChanges)
    inputs=json.loads((OUT/'build_inputs.json').read_text())
    for path,h in inputs.items():assert sha(restored.get(path,Path(path)))==h,path
    exe=BASE/'cfd-build/SU2_CFD/src/SU2_CFD.exe';assert exe.is_file()
    buildBytes=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert buildBytes<3*1024**3
    note=dict(firstPostBuildAudit='failed because two generated dependency-list caches changed; compiler completed normally',
        failedInstallAttempt='compiled.json absent, stopped before compiling helper or installing runtime',
        rerunCompiler=False,actualCompileSourcesAndGeneratedPropertyTablesStable=True,checkedFiles=checked,
        exclusions=cacheChanges,sourceRestored=True,buildBytes=buildBytes,actualNTFSAllocationMeasured=False,
        scriptSHA256=sha(Path(__file__).resolve()))
    (OUT/'build_finalization.json').write_text(json.dumps(note,indent=2)+'\n')
    state=dict(status='compiled-not-yet-runtime-tested',executable=str(exe),executableSHA256=sha(exe),
        coolPropCommit='98b3523d5daa98454618d381d2ae53f7471d216b',sourceStable=True,
        sourceStableScope='all frozen source/property files except two explicitly recorded generated dependency caches; restored build inputs checked against used copies',
        sourceTreeSHA256=sha(OUT/'compile_source_tree.json'),finalizationSHA256=sha(OUT/'build_finalization.json'),physicsQualified=False)
    (OUT/'compiled.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(note,indent=2))
if __name__=='__main__':main()
