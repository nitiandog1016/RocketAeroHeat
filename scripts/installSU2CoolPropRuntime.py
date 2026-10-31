"""编译同库状态查询器并登记独立运行包；仅标为研究能力，不授予物理资格。"""
from pathlib import Path
import hashlib,json,subprocess,os,shutil,re
import buildSU2ThermalRuntime as supervisor
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v121_validation'
OLD=ROOT/'tools/runtime/su2-thermal-build';BASE=ROOT/'tools/runtime/su2-coolprop-build';CP=OLD/'source/subprojects/CoolProp'
FLAVOR='rah-coolprop-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    compiled=json.loads((OUT/'compiled.json').read_text());assert compiled['sourceStable']
    exe=Path(compiled['executable']);assert sha(exe)==compiled['executableSHA256']
    cc=OLD/'toolchain/w64devkit/bin';library=BASE/'cfd-build/subprojects/CoolProp/libCoolProp.a'
    helper=OUT/'RAH_CoolProp_State.exe';source=ROOT/'third_party/su2_thermal/RAHCoolPropState.cpp'
    supervisor.OUT=OUT;env=dict(os.environ,PATH=str(cc)+os.pathsep+os.environ['PATH'])
    locks={str(p):sha(p) for p in [exe,library,source,Path(__file__).resolve()]}
    supervisor.run('state_helper_compile',[str(cc/'g++.exe'),'-std=c++17','-O2','-I'+str(CP/'include'),'-I'+str(CP/'externals/fmtlib/include'),
        str(source),str(library),'-o',str(helper)],ROOT,env,120)
    assert all(sha(Path(p))==h for p,h in locks.items())
    dest=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR);assert not dest.exists();(dest/'bin').mkdir(parents=True)
    for p in [exe,helper]:shutil.copyfile(str(p),str(dest/'bin'/p.name))
    old=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1';base=json.loads((old/'runtime.json').read_text())
    for p in old.iterdir():
        if p.is_file() and p.name not in ['runtime.json','source_patch.diff']:shutil.copyfile(str(p),str(dest/p.name))
    shutil.copyfile(str(old/'source_patch.diff'),str(dest/'axial_wall_source_patch.diff'))
    shutil.copyfile(str(OUT/'build_compatibility.diff'),str(dest/'build_compatibility.diff'))
    # 同时保留CoolProp及实际静态依赖许可；未调用私有REFPROP后端。
    licenseFiles=[CP/'LICENSE',CP/'CoolPropBibTeXLibrary.bib',CP/'dev/fluids/Air.json']
    for folder in [CP/'externals'/n for n in ['Eigen','IF97','REFPROP-headers','fmtlib','msgpack-c','rapidjson','incbin','miniz-3.0.2','nlohmann-json']]+[CP/'boost_CoolProp']:
        licenseFiles += [p for p in folder.iterdir() if p.is_file() and (p.name.lower().startswith(('license','copying','copyright')))]
    for p in licenseFiles:
        if p.exists():q=dest/'licenses'/p.relative_to(CP);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(q))
    dll={}
    for p in (dest/'bin').iterdir():
        data=subprocess.check_output([str(cc/'objdump.exe'),'-p',str(p)]).decode();dll[p.name]=re.findall(r'DLL Name:\s*(\S+)',data)
        assert not any(q.lower().startswith(('libstdc++','libgcc','libgomp','libwinpthread')) for q in dll[p.name]),dll
    state=dict(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(dest/'bin/SU2_CFD.exe'),propertyExecutable=str(dest/'bin/RAH_CoolProp_State.exe'),
       archiveSHA256=base['archiveSHA256'],sourceCommit=base['sourceCommit'],coolPropCommit=compiled['coolPropCommit'],coolPropVersion='7.2.0',
       status='compiled-research-runtime-tests-pending',physicsQualified=False,systemDLLs=dll,buildSourceSHA256=locks,
       restrictions=['single-zone steady compressible NS','CoolProp HEOS Air','SI dimensional','ROE or HLLC','Riemann open boundary','constant or adiabatic wall only initially'],
       files=[dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.rglob('*')) if p.is_file()])
    (dest/'runtime.json').write_text(json.dumps(state,indent=2)+'\n');(OUT/'runtime_installed.json').write_bytes((dest/'runtime.json').read_bytes())
    print(json.dumps(dict(installed=str(dest),runtimeFiles=len(state['files']),physicsQualified=False)))
if __name__=='__main__':main()
