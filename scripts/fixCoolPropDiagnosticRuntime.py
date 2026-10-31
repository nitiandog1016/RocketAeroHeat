"""修复固定CoolProp的一处异常消息缺参；不改EOS公式，保留前两项真实失败。"""
from pathlib import Path
import hashlib,json,os,subprocess,re,difflib,shutil
import buildSU2ThermalRuntime as supervisor
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tools/runtime/su2-coolprop-build'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source';CP=SOURCE/'subprojects/CoolProp'
OUT=ROOT/'outputs/v122_validation/diagnostic_format_build';FLAVOR='rah-coolprop-recon-v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    assert not OUT.exists();OUT.mkdir();supervisor.OUT=OUT
    diagnostic=ROOT/'outputs/v122_validation/exception_diagnostic03'
    assert read(diagnostic/'debugger_process.json')['status']=='terminal'
    assert 'fmt::v12::report_error' in (diagnostic/'debugger.log').read_text(errors='replace')
    path=CP/'src/Backends/Helmholtz/FlashRoutines.cpp';raw=path.read_bytes();text=raw.decode()
    old='format("DP yielded T > 1.5Tmax w/ T (%g) K")'
    assert text.count(old)==1;fixed=text.replace(old,'format("DP yielded T > 1.5Tmax w/ T (%g) K", HEOS._T)')
    (OUT/'FlashRoutines.cpp.original').write_bytes(raw);(OUT/'FlashRoutines.cpp.used').write_bytes(fixed.encode())
    (OUT/'diagnostic_format.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),fixed.splitlines(True),fromfile='a/src/Backends/Helmholtz/FlashRoutines.cpp',tofile='b/src/Backends/Helmholtz/FlashRoutines.cpp')))
    build=BASE/'cfd-build';cc=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'
    ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    env=dict(os.environ,PATH=str(cc)+os.pathsep+os.environ['PATH'])
    obj='subprojects/CoolProp/libCoolProp.a.p/src_Backends_Helmholtz_FlashRoutines.cpp.obj';library=build/'subprojects/CoolProp/libCoolProp.a'
    prior=read(ROOT/'outputs/v122_validation/reconstruction_direct/link_inputs_after.json')
    assert all(sha(Path(p))==h for p,h in prior.items()),'增量链接输入已改变'
    (OUT/'link_inputs_before.json').write_text(json.dumps(prior,indent=2)+'\n')
    def command(target):
        line=subprocess.check_output([str(ninja),'-C',str(build),'-t','commands',target],env=env).decode().splitlines()[-1]
        assert '\\"' not in line
        return [v[1:-1] if v.startswith('"') else v for v in re.findall(r'"[^"]*"|\S+',line)]
    path.write_bytes(fixed.encode())
    try:
        supervisor.run('compile',command(obj),build,env,300)
        supervisor.run('archive_index',[str(cc/'ar.exe'),'s','subprojects/CoolProp/libCoolProp.a'],build,env,30)
        supervisor.run('link',command('SU2_CFD/src/SU2_CFD.exe'),build,env,120)
        helper=OUT/'RAH_CoolProp_State.exe';query=ROOT/'third_party/su2_thermal/RAHCoolPropState.cpp'
        supervisor.run('helper',[str(cc/'g++.exe'),'-std=c++17','-O2','-I'+str(CP/'include'),'-I'+str(CP/'externals/fmtlib/include'),str(query),str(library),'-o',str(helper)],ROOT,env,120)
        changed={str(build/obj),str(library)};assert all(p in changed or sha(Path(p))==h for p,h in prior.items())
        assert path.read_bytes()==fixed.encode()
        (OUT/'link_inputs_after.json').write_text(json.dumps({p:sha(Path(p)) for p in prior},indent=2)+'\n')
        parent=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-recon-v1';manifest=read(parent/'runtime.json')
        dest=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR);assert not dest.exists()
        for e in manifest['files']:
            p=parent/e['path'];assert sha(p)==e['sha256'];q=dest/e['path'];q.parent.mkdir(parents=True,exist_ok=True)
            if e['path']=='bin/SU2_CFD.exe':shutil.copyfile(str(build/'SU2_CFD/src/SU2_CFD.exe'),str(q))
            elif e['path']=='bin/RAH_CoolProp_State.exe':shutil.copyfile(str(helper),str(q))
            else:os.link(str(p),str(q))
        shutil.copyfile(str(OUT/'diagnostic_format.diff'),str(dest/'diagnostic_format.diff'))
        manifest.update(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(dest/'bin/SU2_CFD.exe'),propertyExecutable=str(dest/'bin/RAH_CoolProp_State.exe'),
            parentRuntimeSHA256=sha(parent/'runtime.json'),originalRuntimeSHA256=manifest['parentRuntimeSHA256'],
            diagnosticCorrection='DP_flash OutOfRangeError message supplies its missing temperature argument; numerical EOS unchanged',
            buildSourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),library,query,OUT/'link_inputs_after.json',OUT/'diagnostic_format.diff',build/'SU2_CFD/src/SU2_CFD.exe',helper]})
        manifest['status']='compiled-diagnostic-correction-tests-pending'
        manifest['files']=[dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.rglob('*')) if p.is_file()]
        (dest/'runtime.json').write_text(json.dumps(manifest,indent=2)+'\n');(OUT/'runtime_installed.json').write_bytes((dest/'runtime.json').read_bytes())
        print(json.dumps(dict(installed=str(dest),physicsQualified=False)))
    finally:
        path.write_bytes(raw);(OUT/'source_restore.json').write_text(json.dumps({str(path):sha(path)==hashlib.sha256(raw).hexdigest()})+'\n')
if __name__=='__main__':main()
