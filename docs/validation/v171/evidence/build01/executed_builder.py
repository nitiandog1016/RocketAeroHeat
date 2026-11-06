"""模块用途：准确observed基线复用空间Tw适配，有限编译并恢复全部共享状态。"""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from buildSU2IdealAirRuntime import ROOT,SOURCE,BUILD,BASE
from buildSU2IdealAirFirstSample import arguments
from su2ThermalRuntime import read,sha
from validateIdealAir import run,save

FLAVOR='rah-ideal-air-spatial-wall-v1'
PARENT=ROOT/'outputs/v163_validation/build03'
OLD=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-observed-v1'
PLAN=ROOT/'docs/validation/v171/PLAN.md'
OUT=ROOT/'outputs/v171_validation/build01'
FILES=('Common/src/geometry/CGeometry.cpp','SU2_CFD/src/python_wrapper_structure.cpp',
       'SU2_CFD/src/variables/CNSVariable.cpp')
HEADER='Common/src/geometry/RAHAxialWallTemperature.hpp'
TURB_OBJECT=Path('SU2_CFD/src/libSU2core.a.p/solvers_CTurbSSTSolver.cpp.obj')
CACHED_TURB=ROOT/'tools/runtime/su2-coolprop-build/cfd-build'/TURB_OBJECT
SPATIAL_WALL=dict(mode='axial-node-temperature-table-v1',schema='RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K',
    environmentVariable='RAH_WALL_TEMPERATURE_FILE',marker='wall',appliedReceipt='rah_wall_temperature_applied.csv',
    interpolation='linear-in-X-no-extrapolation',boundaryFormulaChanged=False,supportsNoTableBaseline=True)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def prepared_sources():
    """模块一：由准确已安装包的执行源恢复基线，只更换两个热输入适配函数。"""
    tree=read(PARENT/'source_tree.json')
    previous=read(OLD/'runtime.json')
    assert sha(PARENT/'source_tree.json')==previous['sourceTreeSHA256']
    for entry in previous['files']:assert sha(OLD/entry['path'])==entry['sha256']
    snapshots={SOURCE/p.relative_to(PARENT/'executed_sources'):p.read_bytes()
               for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    assert all(sha_bytes(raw)==tree[str(path)] for path,raw in snapshots.items())
    for name in FILES:
        path=SOURCE/name
        if path not in snapshots:snapshots[path]=path.read_bytes()
        assert sha_bytes(snapshots[path])==tree[str(path)]
    inherited=dict(snapshots)
    header=ROOT/'third_party/su2_thermal/RAHAxialWallTemperature.hpp'
    body=ROOT/'third_party/su2_thermal/SetCustomBoundary.cpp.inc'
    assert sha(header)=='9c09b738ea9327b3571ccbb716bafac2c5c82f723dcb6205d85cd3389f889920'
    assert sha(body)=='e2c0c754a668dba7b9667146d74d71f40c72664978e3d1015c0d2bde15e76dc1'
    text=snapshots[SOURCE/FILES[0]].decode('utf8')
    begin=text.index('void CGeometry::SetCustomBoundary(CConfig* config) {')
    end=text.index('\nvoid CGeometry::UpdateCustomBoundaryConditions',begin)
    text=text[:begin]+body.read_text('utf8')+'\n'+text[end:]
    assert text.count('#include <unordered_set>')==1
    text=text.replace('#include <unordered_set>','#include <unordered_set>\n#include "RAHAxialWallTemperature.hpp"',1)
    snapshots[SOURCE/FILES[0]]=text.encode('utf8')
    marker='void CDriver::PreprocessPythonInterface(CConfig** config, CGeometry**** geometry, CSolver***** solver) {'
    guard='''
  const char* rah_file = std::getenv("RAH_WALL_TEMPERATURE_FILE");
  if (rah_file && *rah_file && (nZone != 1 || config[0]->GetnMarker_PyCustom() != 1))
    SU2_MPI::Error("RAH wall table requires exactly one custom wall marker in one zone", CURRENT_FUNCTION);
'''
    text=snapshots[SOURCE/FILES[1]].decode('utf8');assert text.count(marker)==1
    snapshots[SOURCE/FILES[1]]=('#include <cstdlib>\n'+text.replace(marker,marker+guard,1)).encode('utf8')
    snapshots[SOURCE/HEADER]=header.read_bytes()
    return previous,tree,inherited,snapshots


def main():
    """模块二：三编译/一个逐字缓存恢复/一次链接；不是全量构建，不混入v167数值。"""
    assert PLAN.is_file() and 'v171' in PLAN.read_text('utf8')
    OUT.mkdir(parents=True,exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'plan_at_build.md').write_bytes(PLAN.read_bytes())
    previous,parent_tree,inherited,snapshots=prepared_sources()
    parent=read(PARENT/'objects_after.json')
    before={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    assert len(before)==len(parent)==298 and set(before)==set(parent)
    variable=Path('SU2_CFD/src/libSU2core.a.p/variables_CNSVariable.cpp.obj')
    assert {n for n in before if before[n]!=parent[n]}=={str(TURB_OBJECT),str(variable)}
    assert before==read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    assert sha(CACHED_TURB)==parent[str(TURB_OBJECT)]
    save(OUT/'objects_before.json',before)
    entries=read(BUILD/'compile_commands.json')
    selected=[next(row for row in entries if row['file'].replace('\\','/').endswith('/'+name)) for name in FILES]
    compiled={str(Path(row['output'])) for row in selected}
    changed=compiled|{str(TURB_OBJECT)};assert len(changed)==4
    commands={'compile%d'%i:arguments(row['command']) for i,row in enumerate(selected)}
    commands['link']=read(PARENT/'commands.json')['link'];save(OUT/'commands.json',commands)
    # 上游库是thin archive，编译对象实际被引用；两个库都须审核，不能只核对目录中有obj。
    ar=BASE/'toolchain/w64devkit/bin/ar.exe';archives={}
    for name in ('Common/src/libSU2Common.a','SU2_CFD/src/libSU2core.a'):
        path=BUILD/name;assert path.read_bytes()[:8]==b'!<thin>\n'
        archives[name]=subprocess.check_output([str(ar),'t',str(path)],cwd=str(BUILD)).decode().splitlines()
    assert all(any(Path(member).name==Path(n).name for members in archives.values() for member in members) for n in changed)
    save(OUT/'archive_members.json',archives)
    inputs=[Path(__file__).resolve(),PLAN,ROOT/'third_party/su2_thermal/RAHAxialWallTemperature.hpp',
        ROOT/'third_party/su2_thermal/SetCustomBoundary.cpp.inc',PARENT/'objects_after.json',PARENT/'source_tree.json',
        OLD/'runtime.json',CACHED_TURB,BUILD/'compile_commands.json',ar,BASE/'toolchain/w64devkit/bin/g++.exe']
    input_locks={str(p):sha(p) for p in inputs};save(OUT/'input_manifest.json',input_locks)
    original_sources={p:p.read_bytes() if p.exists() else None for p in snapshots}
    binary_paths={BUILD/n for n in changed}|{BUILD/'SU2_CFD/src/SU2_CFD.exe'}
    for command in commands.values():
        if '-MF' in command:binary_paths.add(BUILD/command[command.index('-MF')+1])
    original_binaries={p:p.read_bytes() if p.exists() else None for p in binary_paths}
    save(OUT/'shared_binary_before.json',{str(p):dict(sha256=sha_bytes(b),bytes=len(b)) if b is not None else None
                                         for p,b in original_binaries.items()})
    target=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    assert not target.exists()
    worker=[sys.executable,str(Path(__file__).resolve()),'--worker']
    try:
        for path,raw in snapshots.items():path.write_bytes(raw)
        altered={str(SOURCE/name) for name in FILES[:2]}
        assert all(sha(Path(n))==h for n,h in parent_tree.items() if n not in altered)
        tree={n:sha(Path(n)) for n in parent_tree};tree[str(SOURCE/HEADER)]=sha(SOURCE/HEADER)
        save(OUT/'source_tree.json',tree)
        (BUILD/TURB_OBJECT).write_bytes(CACHED_TURB.read_bytes())
        save(OUT/'cached_object_restore.json',dict(source=str(CACHED_TURB),target=str(BUILD/TURB_OBJECT),
            sha256=sha(CACHED_TURB),matchesObserved=True))
        started=time.monotonic()
        for i in range(3):
            remaining=120-(time.monotonic()-started);assert remaining>0
            run(OUT,'compile%d'%i,worker+['compile%d'%i],min(60,remaining))
        assert sha(BUILD/variable)==parent[str(variable)],'纯恢复CNSVariable对象未逐字等于observed'
        unchanged={n:sha(BUILD/n)==h==parent[n] for n,h in before.items() if n not in changed}
        assert len(unchanged)==294 and all(unchanged.values())
        inherited_objects={n:sha(BUILD/n)==h for n,h in parent.items() if n not in compiled-{str(variable)}}
        assert len(inherited_objects)==296 and all(inherited_objects.values())
        save(OUT/'unchanged_objects.json',unchanged);save(OUT/'observed_inherited_objects.json',inherited_objects)
        run(OUT,'link',worker+['link'],30)
        after={n:sha(BUILD/n) for n in before};save(OUT/'objects_after.json',after)
        assert all(sha(Path(n))==h for n,h in tree.items())
        assert all(sha(Path(n))==h for n,h in input_locks.items())
        # 保存本包必要对象，可在未来逐字恢复，避免再为相同版本重编；只归档四个，不复制整库。
        for n in changed:
            p=OUT/'produced_objects'/Path(n);p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes((BUILD/n).read_bytes())
        (target/'bin').mkdir(parents=True,exist_ok=False)
        shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'),str(target/'bin/SU2_CFD.exe'))
    finally:
        restored={}
        for path,raw in original_sources.items():
            snap=OUT/'executed_sources'/path.relative_to(SOURCE);snap.parent.mkdir(parents=True,exist_ok=True)
            if path.exists():snap.write_bytes(path.read_bytes())
            if raw is None:
                path.resolve().relative_to(SOURCE.resolve())
                if path.exists():path.unlink()
                restored[str(path)]=not path.exists()
            else:path.write_bytes(raw);restored[str(path)]=path.read_bytes()==raw
        save(OUT/'source_restore.json',restored)
        binary_restore={}
        for path,raw in original_binaries.items():
            if raw is None:
                path.resolve().relative_to(BUILD.resolve())
                if path.exists():path.unlink()
                binary_restore[str(path)]=not path.exists()
            else:path.write_bytes(raw);binary_restore[str(path)]=path.read_bytes()==raw
        save(OUT/'shared_binary_restore.json',binary_restore)
        assert all(restored.values()) and all(binary_restore.values())
        assert {str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}==before
    # 模块三：独立包保留旧物性/观察身份；空间输入能力不自动授予物理或组合数值资格。
    for entry in previous['files']:
        if entry['path'].endswith(('.exe','.diff')):continue
        p=target/entry['path'];p.parent.mkdir(parents=True,exist_ok=True)
        os.link(str(OLD/entry['path']),str(p))
    delta=''
    for name in (*FILES[:2],HEADER):
        path=SOURCE/name
        delta+=''.join(difflib.unified_diff(inherited.get(path,b'').decode('utf8').splitlines(True),
            snapshots[path].decode('utf8').splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (OUT/'spatial_wall.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((OLD/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime=dict(previous)
    runtime.update(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-combination-tests-pending',parentRuntimeSHA256=sha(OLD/'runtime.json'),
        numericalBaselineFlavor=previous['numericalBaselineFlavor'],numericalParentFlavor=previous['flavor'],
        samePhysicsAndNumerics=False,samePhysicsAndNumericsWithoutSpatialTable=True,
        spatialWall=dict(SPATIAL_WALL),incrementalObjectsRebuilt=3,adaptationObjectsRebuilt=2,
        restoredParentObjectsRecompiled=1,restoredParentObjectsFromExactCache=1,
        inheritedVerifiedObjects=296,unchangedSharedObjects=294,allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'),buildReceipt=read(OUT/'link_process.json'),
        buildSourceSHA256={str(p):sha(p) for p in (Path(__file__).resolve(),OUT/'executed_builder.py',
            PLAN,ROOT/'third_party/su2_thermal/RAHAxialWallTemperature.hpp',ROOT/'third_party/su2_thermal/SetCustomBoundary.cpp.inc',
            OUT/'commands.json',OUT/'objects_before.json',OUT/'source_tree.json')},
        files=[dict(path=p.relative_to(target).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json',runtime);save(OUT/'runtime_installed.json',runtime)
    print(json.dumps(dict(installed=str(target),objectsRebuilt=3,cachedParentRestored=1,sharedStateRestored=True)))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]],cwd=str(BUILD)))
    main()
