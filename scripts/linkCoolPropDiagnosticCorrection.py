"""纠正普通静态库误用薄归档索引：真实替换已编译成员并验证二进制确实改变。"""
from pathlib import Path
import hashlib,json,os,subprocess,shutil,argparse
import buildSU2ThermalRuntime as supervisor
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v122_validation/diagnostic_archive_fix'
BASE=ROOT/'tools/runtime/su2-coolprop-build';SOURCE=ROOT/'tools/runtime/su2-thermal-build/source';CP=SOURCE/'subprojects/CoolProp'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--resume-after-archive',action='store_true');args=parser.parse_args()
    if args.resume_after_archive:
        assert read(OUT/'replace_archive_member_process.json')['status']=='terminal'
        assert read(OUT/'replace_archive_member_process.json')['exitCode']==0
    else:assert not OUT.exists();OUT.mkdir()
    supervisor.OUT=OUT
    prior=ROOT/'outputs/v122_validation/diagnostic_format_build';build=BASE/'cfd-build'
    inputs=read(prior/'link_inputs_after.json')
    cc=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin';env=dict(os.environ,PATH=str(cc)+os.pathsep+os.environ['PATH'])
    lib=build/'subprojects/CoolProp/libCoolProp.a';obj=build/'subprojects/CoolProp/libCoolProp.a.p/src_Backends_Helmholtz_FlashRoutines.cpp.obj'
    assert lib.open('rb').read(8)==b'!<arch>\n';oldhash=inputs[str(lib)]
    assert all((args.resume_after_archive and p==str(lib)) or sha(Path(p))==h for p,h in inputs.items())
    probe=OUT/'diagnostic_type.cpp';probe.write_text('''#include "AbstractState.h"
#include "CoolProp.h"
#include "Exceptions.h"
#include <iostream>
#include <memory>
int main() {
  std::unique_ptr<CoolProp::AbstractState> s(CoolProp::AbstractState::factory("HEOS", "Air"));
  try { s->update(CoolProp::DmassP_INPUTS,5000.0/(287.0*4000.0),5000.0); }
  catch (const CoolProp::OutOfRangeError& e) { std::cout << "typed-OutOfRangeError " << e.what(); return 0; }
  catch (const std::exception& e) { std::cout << "other-exception " << e.what(); return 0; }
  std::cout << "unexpected-accept"; return 1;
}
''')
    def compile_source(name,source,dest):
        supervisor.run(name,[str(cc/'g++.exe'),'-std=c++17','-O2','-I'+str(CP/'include'),'-I'+str(CP/'externals/fmtlib/include'),str(source),str(lib),'-o',str(dest)],ROOT,env,120)
    if not args.resume_after_archive:
        compile_source('old_probe_compile',probe,OUT/'old_probe.exe')
        supervisor.run('old_probe',[str(OUT/'old_probe.exe')],ROOT,env,20)
    assert 'other-exception' in (OUT/'old_probe.log').read_text()
    if not args.resume_after_archive:supervisor.run('replace_archive_member',[str(cc/'ar.exe'),'r',str(lib),str(obj)],build,env,30)
    assert sha(lib)!=oldhash
    # 独立提取归档成员，与真实编译对象原字节比较；不能只凭ar退出0。
    # Windows ar p会把二进制LF变CRLF；用独立新目录原生提取，不能对管道字节声称对象不等。
    extracted=OUT/'extracted_member';assert not extracted.exists();extracted.mkdir()
    supervisor.run('extract_member',[str(cc/'ar.exe'),'x',str(lib),obj.name],extracted,env,30)
    assert sha(extracted/obj.name)==sha(obj)
    command=read(prior/'link_process.json')['command'];supervisor.run('link',command,build,env,120)
    helper=OUT/'RAH_CoolProp_State.exe';compile_source('helper',ROOT/'third_party/su2_thermal/RAHCoolPropState.cpp',helper)
    compile_source('new_probe_compile',probe,OUT/'new_probe.exe');supervisor.run('new_probe',[str(OUT/'new_probe.exe')],ROOT,env,20)
    assert 'typed-OutOfRangeError DP yielded T > 1.5Tmax' in (OUT/'new_probe.log').read_text()
    parent=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-recon-v2';manifest=read(parent/'runtime.json')
    assert sha(helper)!=sha(Path(manifest['propertyExecutable'])) and sha(build/'SU2_CFD/src/SU2_CFD.exe')!=sha(Path(manifest['executable']))
    flavor='rah-coolprop-recon-v3';dest=ROOT/'tools/runtime/su2'/('8.6.0-'+flavor);assert not dest.exists()
    for e in manifest['files']:
        p=parent/e['path'];assert sha(p)==e['sha256'];q=dest/e['path'];q.parent.mkdir(parents=True,exist_ok=True)
        if e['path']=='bin/SU2_CFD.exe':shutil.copyfile(str(build/'SU2_CFD/src/SU2_CFD.exe'),str(q))
        elif e['path']=='bin/RAH_CoolProp_State.exe':shutil.copyfile(str(helper),str(q))
        else:os.link(str(p),str(q))
    manifest.update(flavor=flavor,version='8.6.0+'+flavor,executable=str(dest/'bin/SU2_CFD.exe'),propertyExecutable=str(dest/'bin/RAH_CoolProp_State.exe'),
        parentRuntimeSHA256=sha(parent/'runtime.json'),status='diagnostic-type-probed-real-flow-pending',
        buildSourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),lib,obj,prior/'FlashRoutines.cpp.used',OUT/'old_probe.log',OUT/'new_probe.log',helper,build/'SU2_CFD/src/SU2_CFD.exe']})
    manifest['files']=[dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.rglob('*')) if p.is_file()]
    (dest/'runtime.json').write_text(json.dumps(manifest,indent=2)+'\n');(OUT/'runtime_installed.json').write_bytes((dest/'runtime.json').read_bytes())
    report=dict(previousV2BinaryIdenticalToV1=True,previousCorrectionNotLinked=True,archiveMemberVerified=True,typedExceptionVerified=True,newBinariesChanged=True,
        oldLibrarySHA256=oldhash,newLibrarySHA256=sha(lib),originalSourceRestored=sha(CP/'src/Backends/Helmholtz/FlashRoutines.cpp')==sha(prior/'FlashRoutines.cpp.original'),physicsQualified=False)
    assert report['originalSourceRestored'];(OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
