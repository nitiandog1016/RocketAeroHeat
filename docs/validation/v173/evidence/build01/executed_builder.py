"""模块用途：准确空间Tw父包的三对象有限构建，所有终态恢复共享源及二进制。"""
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

from prepareSU2SpatialStepPatch import (
    ROOT, SOURCE, BUILD, BASE, PARENT, OLD, FLAVOR, FILE, SOLVER, EULER,
    CHANGED, prepared_sources,
)
from buildSU2IdealAirFirstSample import arguments
from prepareNACASpatialStep import QUERY_CONTEXT
from su2IdealAirInterface import CONSERVATIVE_STEP
from su2ThermalRuntime import read, sha
from validateIdealAir import run, save

OUT = ROOT/'outputs/v173_validation/build01'
PLAN = ROOT/'docs/validation/v173/PLAN.md'
FIXTURE = ROOT/'outputs/v173_validation/fixture01'
FILES = (FILE, SOLVER, EULER)
CACHED = (
    'Common/src/libSU2Common.a.p/geometry_CGeometry.cpp.obj',
    'SU2_CFD/src/libSU2core.a.p/python_wrapper_structure.cpp.obj',
    'SU2_CFD/src/libSU2core.a.p/solvers_CTurbSSTSolver.cpp.obj',
)
BUDGET = dict(compileMaximumSeconds=120, perObjectMaximumSeconds=60,
              linkMaximumSeconds=30, parentMaximumSeconds=210,
              maximumMemoryBytes=4*1024**3, maximumLinks=1, maximumCFDRuns=0)
QUERY_CONTEXT_DETAILS = dict(
    mode='ideal-air-domain-error-query-context-v1', model='IDEAL_AIR_VARIABLE_CP',
    queries=['rhoe', 'PT', 'Prho', 'rhoT', 'hs', 'Ps'],
    delegatedQueries=['SetEnergy_Prho', 'ComputeDerivativeNRBC_Prho'],
    caughtException='std::domain_error', otherExceptions='propagate-unchanged',
    payload=['query-mode', 'input1', 'input2', 'original-reason', 'phase',
             'localKnown', 'local', 'globalKnown', 'global', 'markerKnown',
             'marker', 'innerKnown', 'inner'],
    numericPrecision=17, numericLocale='classic',
    storage='function-local-static-thread_local-with-nested-scope-restore',
    nodePhases=['node.trial', 'node.old', 'node.backtrack', 'node.zero-step',
                'node.legacy-recovery'],
    nodeIdentity='local-RCM-only; global context separately in conservative_steps.csv',
    riemannPhases=['riemann.internal', 'riemann.exterior', 'riemann.mixed'],
    riemannIdentity='local-RCM-and-global-node-marker-inner-iteration',
    otherCalls='unscoped-with-unknown-identity-flags',
    logging='exception-text-only; existing step trace and termination decide persistence',
    successfulQueriesUnchanged=True, boundaryNumericsChanged=False,
    boundaryAdmissibilityProtection=False, reconstructionAdmissibilityProtection=False,
    changesQueryNumerics=False, physicsQualified=False,
)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def registered_plan():
    """模块一：必须有独立、明确的构建登记；第一阶段只读计划不构成编译许可。"""
    blocks = re.findall(r'```v173-build-registration\s*\n(.*?)\n```',
                        PLAN.read_text(encoding='utf8'), flags=re.S)
    assert len(blocks) == 1, 'PLAN缺唯一v173-build-registration构建登记'
    expected = dict(flavor=FLAVOR, compileSources=list(FILES), cachedObjects=list(CACHED),
                    fixtureReport=FIXTURE.relative_to(ROOT).as_posix()+'/report.json', **BUDGET)
    assert json.loads(blocks[0]) == expected, '构建登记与准确对象/资源合同不同'
    return expected


def restore_paths(originals, boundary, executed=None):
    """模块二：逐路径独立恢复；一次留证失败不能阻止其余原路径恢复。"""
    restored, errors = {}, {}
    for path, raw in originals.items():
        if executed is not None:
            try:
                target = executed/path.relative_to(boundary)
                if path.exists():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(path.read_bytes())
            except Exception as error:
                errors[str(path)+'::executed-source'] = repr(error)
        try:
            path.resolve().relative_to(boundary.resolve())
            if raw is None:
                if path.exists():
                    path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)] = path.read_bytes() == raw
        except Exception as error:
            restored[str(path)] = False
            errors[str(path)+'::restore'] = repr(error)
    return restored, errors


def main():
    """模块三：三缓存恢复、三对象重建、一次链接；独立发布在全部共享恢复之后。"""
    started = time.monotonic()
    registration = registered_plan()
    fixture_report = read(FIXTURE/'report.json')
    assert fixture_report['passed'] is True and fixture_report['budgetPassed'] is True, '新组合夹具未通过'
    assert fixture_report['engineObjectsUnchanged'] == 298
    assert len(fixture_report['sourcesRestored']) == 21 and all(fixture_report['sourcesRestored'].values())
    assert all(fixture_report['inputsStable'].values()) and all(fixture_report['sharedBinariesStable'].values())
    fixture_locks = read(FIXTURE/'input_locks.json')
    assert fixture_locks == fixture_report['inputSHA256']
    fixture_bindings = {}
    for name, digest in fixture_locks.items():
        original = Path(name)
        # 夹具按第一阶段原计划执行；后续构建登记只改变当前计划，不追改旧回执。
        bound = FIXTURE/'project_sources'/original.relative_to(ROOT) if original == PLAN else original
        assert sha(bound) == digest, '夹具准确输入已改变: '+str(original)
        fixture_bindings[name] = dict(boundPath=str(bound), sha256=digest,
            historicalPlanSnapshot=original == PLAN)
    previous, parent_tree, inherited, snapshots = prepared_sources()
    fixture_tree = read(FIXTURE/'executed_source_tree.json')
    assert set(fixture_tree) == {str(p) for p in snapshots}
    assert all(fixture_tree.get(str(p)) == sha_bytes(raw) for p, raw in snapshots.items()), \
        '新夹具执行源与本次准备源不一致'
    assert all(sha(FIXTURE/'executed_sources'/p.relative_to(SOURCE)) == fixture_tree[str(p)]
               for p in snapshots), '夹具准确执行源文件身份不符'
    assert len(inherited) == 19 and len(snapshots) == 21
    assert previous['flavor'] == 'rah-ideal-air-spatial-wall-v1'
    assert previous['asciiSignificantDigits'] == 17
    assert previous['terminationObservation'] == 'cpp-current-exception-windows-stack-delegate-v1'
    target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    target.resolve().relative_to((ROOT/'tools/runtime/su2').resolve())
    assert not target.exists() and not OUT.exists()
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'plan_at_build.md').write_bytes(PLAN.read_bytes())
    save(OUT/'registered_build.json', registration)
    save(OUT/'fixture_binding.json', fixture_bindings)

    def remaining(reserve=15):
        # 为finally中的原状态恢复预留时间；不以终止恢复进程来伪造父预算通过。
        value = BUDGET['parentMaximumSeconds']-(time.monotonic()-started)-reserve
        assert value > 0, '父构建预算耗尽；停止后续阶段并恢复共享状态'
        return value

    parent = read(PARENT/'objects_after.json')
    before = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
    assert len(before) == len(parent) == 298 and set(before) == set(parent)
    assert before == read(PARENT/'objects_before.json'), '共享对象缓存已改变；拒绝混包'
    save(OUT/'objects_before.json', before)
    entries = read(BUILD/'compile_commands.json')
    selected = []
    for name in FILES:
        rows = [row for row in entries if row['file'].replace('\\', '/').endswith('/'+name)]
        assert len(rows) == 1
        selected.append(rows[0])
    compiled = {str(Path(row['output'])) for row in selected}
    cached = {str(Path(name)): PARENT/'produced_objects'/name for name in CACHED}
    touched = compiled | set(cached)
    assert len(compiled) == 3 and len(touched) == 6 and touched.issubset(before)
    assert all(before[n] == parent[n] for n in before if n not in touched)
    assert all(sha(path) == parent[name] for name, path in cached.items())
    commands = {'compile%d'%i: arguments(row['command']) for i, row in enumerate(selected)}
    commands['link'] = read(PARENT/'commands.json')['link']
    save(OUT/'commands.json', commands)

    ar = BASE/'toolchain/w64devkit/bin/ar.exe'
    archives = {}
    archive_paths = [BUILD/name for name in ('Common/src/libSU2Common.a', 'SU2_CFD/src/libSU2core.a')]
    for path in archive_paths:
        assert path.read_bytes()[:8] == b'!<thin>\n'
        archives[str(path.relative_to(BUILD))] = subprocess.check_output(
            [str(ar), 't', str(path)], cwd=str(BUILD), timeout=min(10, remaining())).decode().splitlines()
    assert all(any(Path(member).name == Path(name).name for members in archives.values()
                   for member in members) for name in touched)
    save(OUT/'archive_members.json', archives)

    project_inputs = [Path(__file__).resolve(), PLAN, ROOT/'docs/validation/v173/BUILD_REVIEW.md',
        ROOT/'scripts/prepareSU2SpatialStepPatch.py', ROOT/'scripts/prepareSU2ConservativeStepPatch.py',
        ROOT/'scripts/prepareNACASpatialStep.py',
        ROOT/'scripts/buildSU2IdealAirRuntime.py', ROOT/'scripts/buildSU2IdealAirFirstSample.py',
        ROOT/'scripts/su2IdealAirInterface.py', ROOT/'scripts/su2ThermalRuntime.py',
        ROOT/'scripts/validateIdealAir.py', ROOT/'third_party/ideal_air/ConservativeStep.inc',
        ROOT/'third_party/ideal_air/ConservativeStepTrace.hpp', ROOT/'third_party/ideal_air/PropertyQueryContext.hpp']
    inputs = project_inputs + [PARENT/'objects_after.json', PARENT/'objects_before.json',
        PARENT/'source_tree.json', PARENT/'commands.json', OLD/'runtime.json',
        FIXTURE/'report.json', FIXTURE/'executed_source_tree.json', FIXTURE/'input_locks.json',
        BUILD/'compile_commands.json',
        ar, BASE/'toolchain/w64devkit/bin/g++.exe'] + list(cached.values())
    inputs += [Path(value['boundPath']) for value in fixture_bindings.values()]
    inputs += [FIXTURE/'executed_sources'/p.relative_to(SOURCE) for p in snapshots]
    inputs += [p for p in (PARENT/'executed_sources').rglob('*') if p.is_file()]
    inputs += [OLD/row['path'] for row in previous['files']]
    input_locks = {str(p): sha(p) for p in inputs}
    save(OUT/'input_manifest.json', input_locks)
    for path in project_inputs:
        copy = OUT/'project_sources'/path.relative_to(ROOT)
        copy.parent.mkdir(parents=True, exist_ok=True)
        copy.write_bytes(path.read_bytes())
    for path in snapshots:
        path.resolve().relative_to(SOURCE.resolve())
    original_sources = {p: p.read_bytes() if p.exists() else None for p in snapshots}
    binary_paths = {BUILD/name for name in touched} | {BUILD/'SU2_CFD/src/SU2_CFD.exe'}
    for name, command in commands.items():
        if name.startswith('compile'):
            assert command.count('-MF') == 1
            binary_paths.add(BUILD/command[command.index('-MF')+1])
    assert len(binary_paths) == 10
    for path in binary_paths:
        path.resolve().relative_to(BUILD.resolve())
    original_binaries = {p: p.read_bytes() if p.exists() else None for p in binary_paths}
    save(OUT/'shared_binary_before.json', {str(p): dict(sha256=sha_bytes(raw), bytes=len(raw))
        if raw is not None else None for p, raw in original_binaries.items()})
    save(OUT/'shared_sources_before.json', {str(p): dict(sha256=sha_bytes(raw), bytes=len(raw))
        if raw is not None else None for p, raw in original_sources.items()})
    unchanged_build = archive_paths + [BUILD/name for name in ('build.ninja', '.ninja_deps', '.ninja_log')]
    unchanged_build_locks = {str(p): sha(p) for p in unchanged_build}
    save(OUT/'unchanged_build_inputs.json', unchanged_build_locks)
    worker = [sys.executable, str(Path(__file__).resolve()), '--worker']
    build_complete = False
    try:
        remaining()
        for path, raw in snapshots.items():
            path.write_bytes(raw)
        altered = {str(SOURCE/name) for name in CHANGED}
        assert all(sha(Path(name)) == value for name, value in parent_tree.items() if name not in altered)
        tree = dict(parent_tree)
        tree.update({str(path): sha(path) for path in snapshots})
        save(OUT/'source_tree.json', tree)
        cached_receipt = []
        for name, source in cached.items():
            (BUILD/name).write_bytes(source.read_bytes())
            assert sha(BUILD/name) == parent[name]
            cached_receipt.append(dict(source=str(source), target=str(BUILD/name),
                                       sha256=parent[name], matchesSpatialParent=True))
        save(OUT/'cached_object_restore.json', cached_receipt)
        compilation_started = time.monotonic()
        for i in range(3):
            compile_remaining = BUDGET['compileMaximumSeconds']-(time.monotonic()-compilation_started)
            assert compile_remaining > 0
            run(OUT, 'compile%d'%i, worker+['compile%d'%i],
                min(BUDGET['perObjectMaximumSeconds'], compile_remaining, remaining()))
        assert sum(read(OUT/('compile%d_process.json'%i))['elapsedSeconds'] for i in range(3)) <= 120
        unchanged = {name: sha(BUILD/name) == value for name, value in before.items() if name not in touched}
        inherited_objects = {name: sha(BUILD/name) == value for name, value in parent.items() if name not in compiled}
        assert len(unchanged) == 292 and all(unchanged.values())
        assert len(inherited_objects) == 295 and all(inherited_objects.values())
        save(OUT/'unchanged_objects.json', unchanged)
        save(OUT/'spatial_inherited_objects.json', inherited_objects)
        run(OUT, 'link', worker+['link'], min(BUDGET['linkMaximumSeconds'], remaining()))
        after = {name: sha(BUILD/name) for name in before}
        save(OUT/'objects_after.json', after)
        assert all(sha(Path(name)) == value for name, value in tree.items())
        assert all(sha(Path(name)) == value for name, value in input_locks.items())
        assert all(sha(Path(name)) == value for name, value in unchanged_build_locks.items())
        for name in touched:
            produced = OUT/'produced_objects'/name
            produced.parent.mkdir(parents=True, exist_ok=True)
            produced.write_bytes((BUILD/name).read_bytes())
        shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'), str(OUT/'SU2_CFD.exe'))
        build_complete = True
    finally:
        # 先恢复所有路径，随后才写回执/断言，保证一类失败不会跳过另一类恢复。
        source_restore, source_errors = restore_paths(original_sources, SOURCE, OUT/'executed_sources')
        binary_restore, binary_errors = restore_paths(original_binaries, BUILD)
        restored_objects = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
        save(OUT/'source_restore.json', source_restore)
        save(OUT/'shared_binary_restore.json', binary_restore)
        save(OUT/'restoration_errors.json', dict(source=source_errors, binary=binary_errors))
        save(OUT/'objects_restored.json', restored_objects)
        save(OUT/'build_terminal.json', dict(buildCompleted=build_complete,
            elapsedSeconds=time.monotonic()-started, parentMaximumSeconds=BUDGET['parentMaximumSeconds'],
            sourceRestorePassed=all(source_restore.values()) and not source_errors,
            binaryRestorePassed=all(binary_restore.values()) and not binary_errors,
            all298ObjectsRestored=restored_objects == before))
        assert all(source_restore.values()) and not source_errors
        assert all(binary_restore.values()) and not binary_errors
        assert restored_objects == before
        assert all(sha(Path(name)) == value for name, value in unchanged_build_locks.items())

    # 模块四：恢复完成才发布独立包；不继承父包“无空间表时数值相同”的旧声明。
    remaining(reserve=0)
    assert all(sha(Path(name)) == value for name, value in input_locks.items())
    (target/'bin').mkdir(parents=True, exist_ok=False)
    shutil.copyfile(str(OUT/'SU2_CFD.exe'), str(target/'bin/SU2_CFD.exe'))
    for entry in previous['files']:
        if entry['path'].endswith(('.exe', '.diff')):
            continue
        destination = target/entry['path']
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(str(OLD/entry['path']), str(destination))
    delta = ''
    for name in CHANGED:
        path = SOURCE/name
        delta += ''.join(difflib.unified_diff(inherited.get(path, b'').decode('utf8').splitlines(True),
            snapshots[path].decode('utf8').splitlines(True), fromfile='a/'+name, tofile='b/'+name))
    (OUT/'spatial_step.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((OLD/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime = dict(previous)
    for stale in ('admissibility', 'restoredEulerObjectSHA256', 'samePhysicsAndNumericsWithoutSpatialTable',
                  'mainCompileAttempts'):
        runtime.pop(stale, None)
    runtime.update(flavor=FLAVOR, version='8.6.0+'+FLAVOR, executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-combination-tests-pending', parentRuntimeSHA256=sha(OLD/'runtime.json'),
        numericalBaselineFlavor=FLAVOR, inheritedNumericalBaselineFlavor=previous['numericalBaselineFlavor'],
        numericalParentFlavor=previous['flavor'], samePhysicsAndNumerics=False,
        samePhysicsAndNumericsWithoutSpatialTable=False,
        conservativeStep=dict(CONSERVATIVE_STEP), queryContext=dict(QUERY_CONTEXT),
        queryContextDetails=dict(QUERY_CONTEXT_DETAILS),
        spatialWall=dict(previous['spatialWall']), incrementalObjectsRebuilt=3,
        adaptationObjectsRebuilt=3, restoredParentObjectsRecompiled=0,
        restoredParentObjectsFromExactCache=3, inheritedVerifiedObjects=295,
        unchangedSharedObjects=292, allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'), buildReceipt=read(OUT/'link_process.json'),
        buildBudget=dict(BUDGET), sharedStateRestored=True,
        buildSourceSHA256={str(p): sha(p) for p in project_inputs + [OUT/'executed_builder.py',
            OUT/'commands.json', OUT/'objects_before.json', OUT/'source_tree.json',
            OUT/'input_manifest.json', OUT/'source_restore.json', OUT/'shared_binary_restore.json',
            OUT/'objects_restored.json', OUT/'registered_build.json']},
        files=[dict(path=p.relative_to(target).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    remaining(reserve=0)
    save(target/'runtime.json', runtime)
    save(OUT/'runtime_installed.json', runtime)
    save(OUT/'parent_receipt.json', dict(status='terminal', elapsedSeconds=time.monotonic()-started,
        maximumSeconds=BUDGET['parentMaximumSeconds'], maximumMemoryBytes=BUDGET['maximumMemoryBytes'],
        objectsRebuilt=3, cachedObjectsRestored=3, links=1, sharedStateRestored=True))
    assert time.monotonic()-started <= BUDGET['parentMaximumSeconds']
    print(json.dumps(dict(installed=str(target), objectsRebuilt=3, cachedParentRestored=3,
                          sharedStateRestored=True)))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        assert sys.argv[2] in ('compile0', 'compile1', 'compile2', 'link')
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]], cwd=str(BUILD)))
    assert len(sys.argv) == 1
    main()
