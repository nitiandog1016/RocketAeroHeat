"""v173 小证据归档：只在任务冻结后执行；新构建及准确失败保留，大运行产物不复制。"""

from pathlib import Path
import argparse
import hashlib
import json
import os
import time


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/v173_validation'
DOC = ROOT / 'docs/validation/v173'
DEST = DOC / 'evidence'
MANIFEST = DOC / 'evidence_manifest.json'
RUNTIME = ROOT / 'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-step-v1'
PARENT_RUNTIME = ROOT / 'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
LIMIT = 8 * 1024**2
SMALL_SUFFIXES = {'.py', '.json', '.csv', '.cfg', '.log', '.txt', '.md', '.diff',
                  '.cpp', '.hpp', '.h', '.c', '.inc', '.inl', '.ps1', '.out', '.err', '.sha256', '.png', '.svg'}
NEVER_SUFFIXES = {'.obj', '.o', '.a', '.lib', '.exe', '.dll', '.pdb', '.pdf',
                  '.npz', '.npy', '.su2', '.vtk', '.vtu', '.vtp', '.vtm', '.szplt', '.plt'}
PROOFS = (
    dict(version=172, path=ROOT/'docs/validation/v172/evidence_manifest.json',
         sha256='8947384621aadc08a08f054f4965bf9025466d06709384637441ba3dbc746d08'),
    dict(version=171, path=ROOT/'docs/validation/v171/evidence_manifest.json',
         sha256='b20ee03659b14bbaff1e2b1db3ad9e05cb932f43cf13c45e29054d29a9155435'),
)
SOURCE_NAMES = [
    'docs/PROJECT_GOALS.md', 'AGENTS.md', 'README.md',
    'scripts/prepareSU2SpatialStepPatch.py', 'scripts/buildSU2SpatialStepRuntime.py',
    'scripts/verifySU2SpatialStep.py', 'scripts/reviewNACAAdmissibilityCoverage.py',
    'scripts/prepareNACASpatialStep.py', 'scripts/su2NACASpatialStepInterface.py',
    'scripts/prepareSU2ConservativeStepPatch.py',
    'scripts/prepareNACAAdaptiveContinuation.py', 'scripts/su2NACAAdaptiveInterface.py',
    'scripts/reviewSU2AdaptiveSource.py', 'scripts/reviewNACAColumnResolution.py',
    'scripts/analyzeNACASpatialHeat.py', 'scripts/nacaSavedFrameStudy.py',
    'scripts/su2NACASpatialInterface.py', 'scripts/su2NACALaminarInterface.py',
    'scripts/su2IdealAirCheckpoints.py', 'scripts/su2IdealAirInterface.py',
    'scripts/su2ThermalRuntime.py', 'scripts/runEmbeddedSU2.py',
    'scripts/prepareNACASpatialCase.py', 'scripts/prepareSU2NACALaminar.py',
    'scripts/nacaSpatialWallFlux.py', 'scripts/nacaWallFlux.py',
    'scripts/nacaMeasurementSegments.py', 'scripts/nacaBoundaryLayerEdge.py',
    'scripts/naca3344Model.py', 'scripts/naca3344Normalization.py',
    'scripts/reviewNACASpatialDevelopment.py', 'scripts/plotNACASpatialReview.py',
    'scripts/plotNACAAdaptiveHistory.py', 'scripts/reviewNACAAdaptiveFailure.py',
    'scripts/validateIdealAir.py', 'scripts/idealAirProperties.py',
    'tests/testSU2SpatialStep.cpp',
    'third_party/ideal_air/ConservativeStep.inc',
    'third_party/ideal_air/ConservativeStepTrace.hpp',
    'third_party/ideal_air/PropertyQueryContext.hpp',
    'third_party/ideal_air/CSU2IdealAir.hpp', 'third_party/ideal_air/IdealAir.hpp',
]




def relative(path):
    """所有读写路径必须解析到本工作区；不跟随工作区外符号链接。"""
    path = Path(path)
    path.resolve().relative_to(ROOT.resolve())
    return path.relative_to(ROOT).as_posix()


def metadata(path):
    stat = path.stat()
    return dict(bytes=stat.st_size, mtimeNs=stat.st_mtime_ns)


def sha_file(path):
    """仅最终阶段流式核验本地大文件，不生成副本。"""
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def category(path):
    """小型派生壁面CSV及新旧PNG均保留；原生网格、初场、重启及体表场留本地。"""
    name = path.name.lower()
    suffix = path.suffix.lower()
    if suffix in NEVER_SUFFIXES:
        return 'excluded-binary-mesh-field-or-pdf'
    if name.startswith(('initial_', 'restart_', 'volume_', 'surface_flow')):
        return 'excluded-native-field'
    if '__pycache__' in path.parts or suffix in {'.pyc', '.pyo'}:
        return 'excluded-python-cache'
    if suffix not in SMALL_SUFFIXES and name != '.gitattributes':
        return 'excluded-unlisted-format'
    if path.stat().st_size > LIMIT:
        return 'excluded-above-small-file-limit'
    return 'small-evidence'


def collect(extra_sources=()):
    """只枚举路径和文件元数据；这一阶段不读取输入或大场内容。"""
    found = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file():
            found[relative(path)] = (path, 'files', path.relative_to(OUT))
    for name in list(SOURCE_NAMES)+list(extra_sources):
        path = ROOT/name
        source_name = relative(path)
        if not path.is_file():
            raise RuntimeError('缺少已列明最终源码：'+name)
        found[source_name] = (path, 'finalSources', Path('final_source')/source_name)
    for path in sorted(DOC.rglob('*')):
        if not path.is_file() or DEST in path.parents or path == MANIFEST:
            continue
        found[relative(path)] = (path, 'finalSources', Path('final_source')/relative(path))
    # 新构建只显式登记三个运行包文件；不递归遍历运行包或复制obj/exe。
    for runtime, label, names in (
            (RUNTIME, 'new', ('runtime.json', 'ideal_air.diff', 'bin/SU2_CFD.exe')),
            (PARENT_RUNTIME, 'parent', ('runtime.json', 'bin/SU2_CFD.exe'))):
        for name in names:
            path = runtime/name
            if not path.is_file():
                raise RuntimeError('缺少准确运行时身份：'+relative(path))
            found[relative(path)] = (path, 'finalSources', Path('runtime_identity')/label/name)
    # 后续本轮源码/总结及测试在冻结时自动收录；其他文件可用--final-source明确加入。
    for base, pattern in ((ROOT/'tests', '*V173.py'), (ROOT/'scripts', '*SpatialStep*.py'),
                          (ROOT/'docs', '*V173.md')):
        for path in sorted(base.glob(pattern)):
            if path.is_file():
                found[relative(path)] = (path, 'finalSources', Path('final_source')/relative(path))
    return [dict(source=name, path=p, group=group, target=target,
                 category=category(p), **metadata(p))
            for name, (p, group, target) in sorted(found.items())]


def snapshot(rows):
    return {r['source']: dict(bytes=r['bytes'], mtimeNs=r['mtimeNs']) for r in rows}


def read_json(path):
    relative(path)
    if path.stat().st_size > LIMIT:
        raise RuntimeError('回执超小文件上限：'+relative(path))
    return json.loads(path.read_text(encoding='utf-8-sig'))


def terminal_guards(report_names):
    """不能用命令行确认替代真实终态回执；后处理报告路径须由根任务指定。"""
    required = [OUT/'native_real01_process.json', OUT/'real01/engine_process.json',
                OUT/'native_small01_process.json', OUT/'small01/engine_process.json',
                OUT/'build01_supervision_process.json', OUT/'build01/parent_receipt.json']
    for path in required:
        record = read_json(path)
        if record.get('status') != 'terminal':
            raise RuntimeError('真实任务尚非terminal，拒绝归档：'+relative(path))
    # 只检查实际阶段回执；合同拒绝夹具内复制/构造的process不代表活进程。
    process_paths = set(OUT.glob('*process.json')) | set(OUT.glob('*/engine_process.json'))
    for path in sorted(process_paths):
        record = read_json(path)
        if isinstance(record, dict) and record.get('status') in {'running', 'starting'}:
            raise RuntimeError('仍有运行回执，拒绝归档：'+relative(path))
    reports = []
    for name in report_names:
        path = (ROOT/name).resolve()
        path.relative_to(OUT.resolve())
        if path.suffix.lower() != '.json' or not path.is_file():
            raise RuntimeError('后处理报告必须是v173输出内已存在的JSON：'+str(path))
        read_json(path)
        reports.append(relative(path))
    if not reports:
        raise RuntimeError('未指定已结束的后处理报告，拒绝归档')
    return reports


def build_runtime_guards():
    """新构建只核对准确包身份与保存回执，不遍历/重编共享源树或运行库。"""
    build = OUT/'build01'
    terminal = read_json(build/'build_terminal.json')
    required = ('buildCompleted', 'sourceRestorePassed', 'binaryRestorePassed', 'all298ObjectsRestored')
    if not all(terminal.get(key) is True for key in required):
        raise RuntimeError('新构建或共享恢复未全部通过')
    parent_receipt = read_json(build/'parent_receipt.json')
    supervision = read_json(OUT/'build01_supervision_process.json')
    if (parent_receipt.get('status') != 'terminal' or supervision.get('status') != 'terminal'
            or supervision.get('exitCode') != 0 or parent_receipt.get('objectsRebuilt') != 3
            or parent_receipt.get('cachedObjectsRestored') != 3 or parent_receipt.get('links') != 1
            or parent_receipt.get('sharedStateRestored') is not True):
        raise RuntimeError('准确构建父回执身份不符')
    source_restore = read_json(build/'source_restore.json')
    binary_restore = read_json(build/'shared_binary_restore.json')
    before = read_json(build/'objects_before.json')
    restored = read_json(build/'objects_restored.json')
    after = read_json(build/'objects_after.json')
    if (len(source_restore) != 21 or not all(value is True for value in source_restore.values())
            or len(binary_restore) != 10 or not all(value is True for value in binary_restore.values())
            or len(before) != 298 or restored != before or set(after) != set(before)):
        raise RuntimeError('21源/10二进制/298对象的准确恢复回执不相容')

    runtime_path = RUNTIME/'runtime.json'
    runtime = read_json(runtime_path)
    runtime_digest = sha_file(runtime_path)
    if runtime_digest != sha_file(build/'runtime_installed.json'):
        raise RuntimeError('安装运行包与准确构建快照不同')
    if (runtime.get('flavor') != 'rah-ideal-air-spatial-step-v1'
            or runtime.get('numericalParentFlavor') != 'rah-ideal-air-spatial-wall-v1'
            or runtime.get('samePhysicsAndNumerics') is not False
            or runtime.get('samePhysicsAndNumericsWithoutSpatialTable') is not False
            or runtime.get('sharedStateRestored') is not True):
        raise RuntimeError('新数值包身份或等价声明不符')
    if runtime['sourceTreeSHA256'] != sha_file(build/'source_tree.json'):
        raise RuntimeError('准确构建源树SHA不符')
    if runtime['parentRuntimeSHA256'] != sha_file(PARENT_RUNTIME/'runtime.json'):
        raise RuntimeError('准确spatial-wall父运行包身份变化')
    if runtime['buildReceipt'] != read_json(build/'link_process.json'):
        raise RuntimeError('运行包与唯一准确链接回执不同')
    identities = []
    for package, manifest in ((RUNTIME, runtime), (PARENT_RUNTIME, read_json(PARENT_RUNTIME/'runtime.json'))):
        for name in ('bin/SU2_CFD.exe', 'ideal_air.diff') if package == RUNTIME else ('bin/SU2_CFD.exe',):
            records = [row for row in manifest['files'] if row['path'] == name]
            if len(records) != 1:
                raise RuntimeError('运行包目标身份不唯一：'+name)
            path = package/name
            digest = sha_file(path)
            if digest != records[0]['sha256'] or path.stat().st_size != records[0]['bytes']:
                raise RuntimeError('准确运行包文件SHA/长度不符：'+relative(path))
            identities.append(dict(path=relative(path), bytes=path.stat().st_size,
                sha256=digest, archiveClass='small-evidence' if name.endswith('.diff') else 'local-only'))
    if sha_file(build/'SU2_CFD.exe') != next(row['sha256'] for row in runtime['files']
                                          if row['path'] == 'bin/SU2_CFD.exe'):
        raise RuntimeError('构建保留exe与安装exe不相等')
    produced_names = (
        'Common/src/libSU2Common.a.p/geometry_CGeometry.cpp.obj',
        'SU2_CFD/src/libSU2core.a.p/python_wrapper_structure.cpp.obj',
        'SU2_CFD/src/libSU2core.a.p/solvers_CTurbSSTSolver.cpp.obj',
        'SU2_CFD/src/libSU2core.a.p/variables_CNSVariable.cpp.obj',
        'SU2_CFD/src/libSU2core.a.p/solvers_CNSSolver.cpp.obj',
        'SU2_CFD/src/libSU2core.a.p/solvers_CEulerSolver.cpp.obj',
    )
    produced = []
    for name in produced_names:
        path = build/'produced_objects'/name
        digest = sha_file(path)
        if digest != after[str(Path(name))]:
            raise RuntimeError('准确保存对象与新构建记录不同：'+name)
        produced.append(dict(path=relative(path), bytes=path.stat().st_size,
                             sha256=digest, localOnly=True, copied=False))
    return dict(runtimePath=relative(runtime_path), runtimeSHA256=runtime_digest,
        sourceTreeSHA256=runtime['sourceTreeSHA256'], parentRuntimeSHA256=runtime['parentRuntimeSHA256'],
        files=identities, producedObjects=produced, objectsRebuilt=3, cachedObjectsRestored=3, links=1,
        savedRestorationReceiptsVerified=True, currentSharedRuntimeTreeReaudited=False,
        runtimeRecursiveInventory=False, newBuildLaunchedByArchive=False)


def write_new(path, value):
    """排他创建，不覆盖已有失败、报告或归档文件。"""
    relative(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(value)


def execute(args):
    started = time.monotonic()
    if not args.final:
        rows = collect(args.final_source)
        selected = [r for r in rows if r['category'] == 'small-evidence']
        excluded = [r for r in rows if r['category'] != 'small-evidence']
        print(json.dumps(dict(mode='metadata-preview-only', eligibleSmallFiles=len(selected),
              eligibleLogicalBytes=sum(r['bytes'] for r in selected),
              excludedLocalFiles=len(excluded), excludedLogicalBytes=sum(r['bytes'] for r in excluded),
              maximumSmallFileBytes=LIMIT, contentsRead=False, finalArchiveCreated=False,
              inheritedProofs=[dict(version=p['version'], path=relative(p['path']),
                  bytes=p['path'].stat().st_size, contentsRead=False) for p in PROOFS],
              runtimeRecursiveInventory=False,
              localOnly=[dict(path=r['source'], bytes=r['bytes'], reason=r['category']) for r in excluded]),
              ensure_ascii=False, indent=2))
        return
    if not args.confirmed_postprocessing_complete:
        raise RuntimeError('须收到根任务终态及后处理完成通知后才能执行最终归档')
    if DEST.exists() or MANIFEST.exists():
        raise RuntimeError('归档目录或清单已存在；禁止覆盖或自动清理')
    reports = terminal_guards(args.postprocessing_report)
    build_identity = build_runtime_guards()
    inherited = {}
    for entry in PROOFS:
        if sha_file(entry['path']) != entry['sha256']:
            raise RuntimeError('既有证据清单SHA已改变：'+relative(entry['path']))
        proof = read_json(entry['path'])
        for record in proof['files']+proof['finalSources']:
            inherited.setdefault(record['sha256'], record['archive'])
        inherited.setdefault(entry['sha256'], relative(entry['path']))
    inherited_used = {}
    rows = collect(args.final_source)
    before = snapshot(rows)
    seen = {}
    small_records = []
    final_records = []
    local_only = []
    pending = []
    # 在创建目标前完成来源hash和元数据稳定检查，保留所有失败目录原路径。
    for row in rows:
        path = row['path']
        if row['category'] != 'small-evidence':
            local_only.append(dict(path=row['source'], bytes=row['bytes'], sha256=sha_file(path),
                                   localOnly=True, reason=row['category'], copied=False))
            continue
        value = path.read_bytes()
        digest = hashlib.sha256(value).hexdigest()
        if digest not in seen:
            if digest in inherited:
                # 逐字核验被引用的已有证据；不复制旧字段、源码或原构建证据树。
                previous = ROOT/inherited[digest]
                relative(previous)
                if sha_file(previous) != digest:
                    raise RuntimeError('既有证据SHA不符：'+relative(previous))
                seen[digest] = relative(previous)
                inherited_used[digest] = relative(previous)
            else:
                target = DEST/row['target']
                seen[digest] = relative(target)
                pending.append((target, value))
        record = dict(bytes=len(value), sha256=digest, archive=seen[digest])
        if row['group'] == 'files':
            small_records.append(dict(source=row['source'], **record))
        else:
            final_records.append(dict(path=row['source'], **record))
    if snapshot(collect(args.final_source)) != before:
        raise RuntimeError('归档读取期间源文件集合或元数据变化，拒绝创建最终归档')
    terminal_guards(args.postprocessing_report)
    DEST.mkdir(parents=True, exist_ok=False)
    for target, value in pending:
        write_new(target, value)
    # 新归档逐字hash核验；读取不改源目录，未知失败均不清理。
    for record in small_records+final_records:
        if sha_file(ROOT/record['archive']) != record['sha256']:
            raise RuntimeError('归档逐字核验失败：'+record['archive'])
    if snapshot(collect(args.final_source)) != before or any(sha_file(p['path']) != p['sha256'] for p in PROOFS):
        raise RuntimeError('归档写入期间源文件集合/元数据/既有proof变化，保留已写证据并拒绝最终清单')
    report = dict(files=small_records, finalSources=final_records, localOnly=local_only,
        distinctArchivedFiles=len(pending), archivedBytes=sum(len(v) for _,v in pending),
        distinctReferencedSHA256=len(seen), inheritedProofDistinctReferences=len(inherited_used),
        inheritedProofs=[dict(version=p['version'], path=relative(p['path']),
            sha256=p['sha256'], copied=False) for p in PROOFS],
        inheritedArchiveReferences=[dict(sha256=h,archive=name) for h,name in sorted(inherited_used.items())],
        runtimeRecursiveInventory=False, runtimeBuiltInVersion=173, runtimeParentVersion=171,
        buildRuntimeIdentity=build_identity, binaryObjectOrExecutableCopies=0,
        logicalSmallEvidenceFiles=len(small_records)+len(final_records),
        logicalSmallEvidenceBytes=sum(r['bytes'] for r in small_records+final_records),
        sameSHA256ReferenceCount=len(small_records)+len(final_records)-len(seen),
        localOnlyFiles=len(local_only), localOnlyLogicalBytes=sum(r['bytes'] for r in local_only),
        maximumSmallFileBytes=LIMIT, sourcesStable=True, archiveSHA256Verified=True,
        postprocessingReports=reports, uniqueFilesDeleted=0, inputFilesModified=0,
        newCFDTasks=0, newRuntimeBuilds=0, newPhysicalQualification=False,
        operationScope='archive-only; task CFD/build counts remain in their original receipts',
        elapsedSeconds=time.monotonic()-started)
    write_new(MANIFEST, (json.dumps(report, ensure_ascii=False, indent=2)+'\n').encode('utf-8'))
    print(json.dumps({key: value for key,value in report.items()
                      if key not in {'files','finalSources','localOnly','inheritedArchiveReferences'}}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--final', action='store_true')
    parser.add_argument('--confirmed-postprocessing-complete', action='store_true')
    parser.add_argument('--postprocessing-report', action='append', default=[])
    parser.add_argument('--final-source', action='append', default=[])
    execute(parser.parse_args())
