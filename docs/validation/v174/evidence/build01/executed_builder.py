"""模块用途：唯一输出对象有限构建，恢复全部共享状态后发布独立原生定比热包。"""
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from prepareSU2NativeGasOutputPatch import (ROOT, SOURCE, BUILD, BASE, PARENT, OLD,
    FLAVOR, FILE, OBJECT, CACHE_CANDIDATES, OUTPUT_CONTRACT, MODEL, GAMMA, R,
    prepared_sources, object_plan)
from buildSU2IdealAirFirstSample import arguments
from su2ThermalRuntime import read, sha
from validateIdealAir import run, save

OUT = ROOT/'outputs/v174_validation/build01'
PLAN = ROOT/'docs/validation/v174/PLAN.md'
REVIEW = ROOT/'docs/validation/v174/BUILD_REVIEW.md'
BUDGET = dict(compileMaximumSeconds=60, linkMaximumSeconds=30,
    auxiliaryMaximumSeconds=20, parentMaximumSeconds=110,
    maximumMemoryBytes=4*1024**3, maximumObjects=1, maximumLinks=1, maximumCFDRuns=0)


def registration_specification():
    """模块一：返回唯一登记内容；普通import不会构建或覆盖共享文件。"""
    return dict(flavor=FLAVOR, parentFlavor='rah-ideal-air-spatial-wall-v1',
        compileSources=[FILE], cachedObjectCandidates=list(CACHE_CANDIDATES),
        exactSourceOverlayCount=19, inheritedParentObjects=297,
        buildReview=REVIEW.relative_to(ROOT).as_posix(), **BUDGET)


def registered_plan():
    blocks = re.findall(r'```v174-output-build-registration\s*\n(.*?)\n```',
                        PLAN.read_text(encoding='utf8'), flags=re.S)
    assert len(blocks) == 1, 'PLAN缺唯一输出构建登记；禁止启动编译'
    assert json.loads(blocks[0]) == registration_specification(), '输出构建预算或对象登记不符'
    assert REVIEW.is_file()
    return registration_specification()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def restore_paths(originals, boundary, executed=None):
    """模块二：独立恢复每个原路径；保存证据失败也不能跳过后续恢复。"""
    restored, errors = {}, {}
    for path, raw in originals.items():
        if executed is not None:
            try:
                destination = executed/path.relative_to(boundary)
                if path.exists():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(path.read_bytes())
            except Exception as error:
                errors[str(path)+'::executed-source'] = repr(error)
        try:
            path.resolve().relative_to(boundary.resolve())
            if raw is None:
                if path.exists(): path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)] = path.read_bytes() == raw
        except Exception as error:
            restored[str(path)] = False
            errors[str(path)+'::restore'] = repr(error)
    return restored, errors


def main():
    """模块三：单对象/单链接；准确继承297父对象，不递归复制任何运行包。"""
    started = time.monotonic(); process_seconds = {}
    registration = registered_plan()
    previous, parent_tree, inherited, snapshots = prepared_sources()
    parent, before, cached = object_plan()
    target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    target.resolve().relative_to((ROOT/'tools/runtime/su2').resolve())
    assert not target.exists() and not OUT.exists(), '独立输出包/构建证据已存在；拒绝覆盖或重编'
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'plan_at_build.md').write_bytes(PLAN.read_bytes())
    save(OUT/'registered_build.json', registration)
    save(OUT/'objects_before.json', before)

    def remaining(reserve=2):
        elapsed = time.monotonic()-started
        assert elapsed-sum(process_seconds.values()) <= BUDGET['auxiliaryMaximumSeconds'], '辅助阶段预算耗尽'
        value = BUDGET['parentMaximumSeconds']-elapsed-reserve
        assert value > 0, '总预算耗尽；停止后续工作并恢复共享状态'
        return value

    entries = read(BUILD/'compile_commands.json')
    selected = [row for row in entries if row['file'].replace('\\','/').endswith('/'+FILE)]
    assert len(selected) == 1 and str(Path(selected[0]['output'])) == str(Path(OBJECT))
    compiled = str(Path(OBJECT)); touched = set(cached) | {compiled}
    assert all(before[name] == parent[name] for name in parent if name not in touched)
    commands = dict(compile0=arguments(selected[0]['command']), link=read(PARENT/'commands.json')['link'])
    save(OUT/'commands.json', commands)
    ar = BASE/'toolchain/w64devkit/bin/ar.exe'
    archives = [BUILD/name for name in ('Common/src/libSU2Common.a','SU2_CFD/src/libSU2core.a')]
    archive_members = {}
    for path in archives:
        assert path.read_bytes()[:8] == b'!<thin>\n'
        archive_members[str(path.relative_to(BUILD))] = subprocess.check_output(
            [str(ar),'t',str(path)], cwd=str(BUILD), timeout=min(5,remaining())).decode().splitlines()
    assert all(any(Path(member).name == Path(name).name for members in archive_members.values()
                   for member in members) for name in touched)
    save(OUT/'archive_members.json', archive_members)

    project_inputs = [Path(__file__).resolve(), PLAN, REVIEW] + [ROOT/'scripts'/name for name in (
        'prepareSU2NativeGasOutputPatch.py','nacaConstantCpModel.py','naca3344Model.py',
        'buildSU2IdealAirRuntime.py','buildSU2IdealAirFirstSample.py','su2ThermalRuntime.py','validateIdealAir.py')]
    inputs = project_inputs + [PARENT/'objects_after.json',PARENT/'commands.json',PARENT/'source_tree.json',
        OLD/'runtime.json',BUILD/'compile_commands.json',ar,BASE/'toolchain/w64devkit/bin/g++.exe']
    inputs += list(cached.values()) + [path for path in (PARENT/'executed_sources').rglob('*') if path.is_file()]
    inputs += [OLD/item['path'] for item in previous['files']]
    locks = {str(path): sha(path) for path in inputs}
    save(OUT/'input_manifest.json', locks)
    for path in project_inputs:
        copy = OUT/'project_sources'/path.relative_to(ROOT)
        copy.parent.mkdir(parents=True, exist_ok=True); copy.write_bytes(path.read_bytes())
    for path in snapshots: path.resolve().relative_to(SOURCE.resolve())
    original_sources = {p: p.read_bytes() if p.exists() else None for p in snapshots}
    binary_paths = {BUILD/name for name in touched} | {BUILD/'SU2_CFD/src/SU2_CFD.exe'}
    command = commands['compile0']; assert command.count('-MF') == 1
    binary_paths.add(BUILD/command[command.index('-MF')+1])
    for path in binary_paths: path.resolve().relative_to(BUILD.resolve())
    original_binaries = {p: p.read_bytes() if p.exists() else None for p in binary_paths}
    for name, values in [('shared_sources_before',original_sources),('shared_binary_before',original_binaries)]:
        save(OUT/(name+'.json'), {str(p): dict(sha256=digest(raw),bytes=len(raw)) if raw is not None else None
                                  for p, raw in values.items()})
    unchanged_build = archives + [BUILD/name for name in ('build.ninja','.ninja_deps','.ninja_log')]
    unchanged_locks = {str(p): sha(p) for p in unchanged_build}
    save(OUT/'unchanged_build_inputs.json', unchanged_locks)
    worker = [sys.executable,str(Path(__file__).resolve()),'--worker']
    completed = False

    def invoke(name, maximum):
        try:
            run(OUT,name,worker+[name],min(maximum,remaining()))
        finally:
            path = OUT/(name+'_process.json')
            if path.exists(): process_seconds[name] = read(path).get('elapsedSeconds',0.)

    try:
        remaining()
        for path, raw in snapshots.items(): path.write_bytes(raw)
        tree = dict(parent_tree); tree[str(SOURCE/FILE)] = sha(SOURCE/FILE)
        assert all(sha(Path(name)) == value for name,value in tree.items()), '准确父源覆盖之外发生变化'
        save(OUT/'source_tree.json', tree)
        cached_receipt = []
        for name, path in cached.items():
            (BUILD/name).write_bytes(path.read_bytes()); assert sha(BUILD/name) == parent[name]
            cached_receipt.append(dict(source=str(path),target=str(BUILD/name),sha256=parent[name]))
        save(OUT/'cached_object_restore.json', cached_receipt)
        invoke('compile0',BUDGET['compileMaximumSeconds'])
        parent_objects = {name: sha(BUILD/name) == value for name,value in parent.items() if name != compiled}
        assert len(parent_objects) == 297 and all(parent_objects.values())
        save(OUT/'inherited_parent_objects.json', parent_objects)
        invoke('link',BUDGET['linkMaximumSeconds'])
        after = {name: sha(BUILD/name) for name in before}; save(OUT/'objects_after.json',after)
        assert all(sha(Path(name)) == value for name,value in tree.items())
        assert all(sha(Path(name)) == value for name,value in locks.items())
        assert all(sha(Path(name)) == value for name,value in unchanged_locks.items())
        produced = OUT/'produced_objects'/OBJECT; produced.parent.mkdir(parents=True,exist_ok=True)
        produced.write_bytes((BUILD/OBJECT).read_bytes())
        shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'),str(OUT/'SU2_CFD.exe'))
        completed = True
    finally:
        sources_restored, source_errors = restore_paths(original_sources,SOURCE,OUT/'executed_sources')
        binaries_restored, binary_errors = restore_paths(original_binaries,BUILD)
        restored_objects = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
        save(OUT/'source_restore.json',sources_restored)
        save(OUT/'shared_binary_restore.json',binaries_restored)
        save(OUT/'restoration_errors.json',dict(source=source_errors,binary=binary_errors))
        save(OUT/'objects_restored.json',restored_objects)
        save(OUT/'build_terminal.json',dict(buildCompleted=completed,elapsedSeconds=time.monotonic()-started,
            processSeconds=process_seconds,auxiliarySeconds=time.monotonic()-started-sum(process_seconds.values()),
            sourceRestorePassed=all(sources_restored.values()) and not source_errors,
            binaryRestorePassed=all(binaries_restored.values()) and not binary_errors,
            all298ObjectsRestored=restored_objects == before))
        assert len(sources_restored) == 19 and all(sources_restored.values()) and not source_errors
        assert all(binaries_restored.values()) and not binary_errors and restored_objects == before
        assert all(sha(Path(name)) == value for name,value in unchanged_locks.items())

    # 模块四：仅发布新二进制、明确许可文件和增量补丁；不继承递归目录副本。
    remaining()
    assert all(sha(Path(name)) == value for name,value in locks.items())
    (target/'bin').mkdir(parents=True,exist_ok=False)
    shutil.copyfile(str(OUT/'SU2_CFD.exe'),str(target/'bin/SU2_CFD.exe'))
    for entry in previous['files']:
        if entry['path'].endswith(('.exe','.diff')): continue
        destination = target/entry['path']; destination.parent.mkdir(parents=True,exist_ok=True)
        os.link(str(OLD/entry['path']),str(destination))
    delta = ''.join(difflib.unified_diff(inherited[SOURCE/FILE].decode('utf8').splitlines(True),
        snapshots[SOURCE/FILE].decode('utf8').splitlines(True),fromfile='a/'+FILE,tofile='b/'+FILE))
    (OUT/'native_ideal_gas_output.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((OLD/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime = dict(previous)
    for key in ('temperatureRangeK','samePhysicsAndNumericsWithoutSpatialTable','mainCompileAttempts',
                'restoredEulerObjectSHA256','allCFDObjectsRebuilt'):
        runtime.pop(key,None)
    runtime.update(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-output-extension-tests-pending',parentRuntimeSHA256=sha(OLD/'runtime.json'),
        numericalParentFlavor=previous['flavor'],numericalBaselineFlavor=previous['numericalBaselineFlavor'],
        inheritedPackageThermodynamicModel=previous['thermodynamicModel'],
        inheritedPackageTemperatureRangeK=previous['temperatureRangeK'],
        thermodynamicModel=MODEL,supportedModel='IDEAL_GAS',gasConstant=R,gamma=GAMMA,
        modelMetadataMeaning='supported selected native IDEAL_GAS case; inherited variable-cp implementation retained',
        nativeIdealGasOutputContract=dict(OUTPUT_CONTRACT),firstFrameOutputModelExtension='IDEAL_GAS',
        samePhysicsAndNumerics=False,numericalUpdateSourceUnchanged=True,sameNumericalUpdateSources=True,
        outputScheduleChanged=True,nativeSharedFrameEqualityVerified=False,
        incrementalObjectsRebuilt=1,adaptationObjectsRebuilt=1,restoredParentObjectsRecompiled=0,
        restoredParentObjectsFromExactCache=len(cached),inheritedVerifiedObjects=297,
        unchangedSharedObjects=298-len(touched),allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'),buildReceipt=read(OUT/'link_process.json'),
        buildBudget=dict(BUDGET),sharedStateRestored=True,
        buildSourceSHA256={str(p):sha(p) for p in project_inputs+[OUT/'executed_builder.py',OUT/'commands.json',
            OUT/'objects_before.json',OUT/'source_tree.json',OUT/'input_manifest.json',OUT/'source_restore.json',
            OUT/'shared_binary_restore.json',OUT/'objects_restored.json',OUT/'registered_build.json']},
        files=[dict(path=p.relative_to(target).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    remaining(reserve=0)
    save(target/'runtime.json',runtime); save(OUT/'runtime_installed.json',runtime)
    elapsed = time.monotonic()-started
    save(OUT/'parent_receipt.json',dict(status='terminal',elapsedSeconds=elapsed,
        auxiliarySeconds=elapsed-sum(process_seconds.values()),processSeconds=process_seconds,
        maximumSeconds=BUDGET['parentMaximumSeconds'],maximumMemoryBytes=BUDGET['maximumMemoryBytes'],
        objectsRebuilt=1,cachedObjectsRestored=len(cached),links=1,sharedStateRestored=True))
    remaining(reserve=0)
    print(json.dumps(dict(installed=str(target),objectsRebuilt=1,cachedObjectsRestored=len(cached),sharedStateRestored=True)))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        assert sys.argv[2] in ('compile0','link')
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]],cwd=str(BUILD)))
    assert len(sys.argv) == 1
    main()
