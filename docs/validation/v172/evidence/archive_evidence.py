"""v172 临时小证据归档：预览不读大场，最终归档须等真实任务和后处理结束。"""

from pathlib import Path
import argparse
import hashlib
import json
import os
import time


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/v172_validation'
DOC = ROOT / 'docs/validation/v172'
DEST = DOC / 'evidence'
MANIFEST = DOC / 'evidence_manifest.json'
RUNTIME = ROOT / 'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
LIMIT = 8 * 1024**2
SMALL_SUFFIXES = {'.py', '.json', '.csv', '.cfg', '.log', '.txt', '.md', '.diff',
                  '.cpp', '.hpp', '.h', '.c', '.inc', '.ps1', '.out', '.err', '.sha256', '.png'}
NEVER_SUFFIXES = {'.obj', '.o', '.a', '.lib', '.exe', '.dll', '.pdb', '.pdf',
                  '.npz', '.npy', '.su2', '.vtk', '.vtu', '.vtp', '.vtm', '.szplt', '.plt'}
PROOF = ROOT / 'docs/validation/v171/evidence_manifest.json'
PROOF_SHA256 = 'b20ee03659b14bbaff1e2b1db3ad9e05cb932f43cf13c45e29054d29a9155435'
SOURCE_NAMES = [
    'docs/NACA_ADAPTIVE_DEVELOPMENT_V172.md', 'docs/PROJECT_GOALS.md',
    'AGENTS.md', 'README.md',
    'scripts/prepareNACAAdaptiveContinuation.py', 'scripts/su2NACAAdaptiveInterface.py',
    'scripts/reviewSU2AdaptiveSource.py', 'scripts/reviewNACAColumnResolution.py',
    'tests/verifyNACAAdaptiveContractV172.py', 'tests/verifyNACAAdaptiveRunV172.py',
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
    'tests/verifyNACASavedFrameStudyV171.py',
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
    if suffix not in SMALL_SUFFIXES:
        return 'excluded-unlisted-format'
    if path.stat().st_size > LIMIT:
        return 'excluded-above-small-file-limit'
    return 'small-evidence'


def collect():
    """只枚举路径和文件元数据；这一阶段不读取输入或大场内容。"""
    found = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file():
            found[relative(path)] = (path, 'files', path.relative_to(OUT))
    for name in SOURCE_NAMES:
        path = ROOT/name
        if not path.is_file():
            raise RuntimeError('缺少已列明最终源码：'+name)
        found[relative(path)] = (path, 'finalSources', Path('final_source')/name)
    for path in sorted(DOC.rglob('*')):
        if not path.is_file() or DEST in path.parents or path == MANIFEST:
            continue
        found[relative(path)] = (path, 'finalSources', Path('final_source')/relative(path))
    # 本轮没有新构建；只登记两个准确运行时身份，不递归遍历或复制运行包。
    for path in (RUNTIME/'runtime.json', RUNTIME/'bin/SU2_CFD.exe'):
        if not path.is_file():
            raise RuntimeError('缺少继承运行时身份：'+relative(path))
        found[relative(path)] = (path, 'finalSources', Path('runtime_identity')/path.relative_to(RUNTIME))
    # 后续新增的本轮测试自动收录；源码新增项由根任务冻结前补入SOURCE_NAMES。
    for path in sorted((ROOT/'tests').glob('*V172.py')):
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
    required = [OUT/'native_real01_process.json', OUT/'real01/engine_process.json']
    for path in required:
        record = read_json(path)
        if record.get('status') != 'terminal':
            raise RuntimeError('真实任务尚非terminal，拒绝归档：'+relative(path))
    for path in sorted(OUT.rglob('*process.json')):
        record = read_json(path)
        if isinstance(record, dict) and record.get('status') in {'running', 'starting'}:
            raise RuntimeError('仍有运行回执，拒绝归档：'+relative(path))
    reports = []
    for name in report_names:
        path = (ROOT/name).resolve()
        path.relative_to(OUT.resolve())
        if path.suffix.lower() != '.json' or not path.is_file():
            raise RuntimeError('后处理报告必须是v172输出内已存在的JSON：'+str(path))
        read_json(path)
        reports.append(relative(path))
    if not reports:
        raise RuntimeError('未指定已结束的后处理报告，拒绝归档')
    return reports


def write_new(path, value):
    """排他创建，不覆盖已有失败、报告或归档文件。"""
    relative(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(value)


def execute(args):
    started = time.monotonic()
    if not args.final:
        rows = collect()
        selected = [r for r in rows if r['category'] == 'small-evidence']
        excluded = [r for r in rows if r['category'] != 'small-evidence']
        print(json.dumps(dict(mode='metadata-preview-only', eligibleSmallFiles=len(selected),
              eligibleLogicalBytes=sum(r['bytes'] for r in selected),
              excludedLocalFiles=len(excluded), excludedLogicalBytes=sum(r['bytes'] for r in excluded),
              maximumSmallFileBytes=LIMIT, contentsRead=False, finalArchiveCreated=False,
              inheritedProof=dict(path=relative(PROOF),bytes=PROOF.stat().st_size,contentsRead=False),
              runtimeRecursiveInventory=False,
              localOnly=[dict(path=r['source'], bytes=r['bytes'], reason=r['category']) for r in excluded]),
              ensure_ascii=False, indent=2))
        return
    if not args.confirmed_postprocessing_complete:
        raise RuntimeError('须收到根任务终态及后处理完成通知后才能执行最终归档')
    if DEST.exists() or MANIFEST.exists():
        raise RuntimeError('归档目录或清单已存在；禁止覆盖或自动清理')
    reports = terminal_guards(args.postprocessing_report)
    if sha_file(PROOF) != PROOF_SHA256:
        raise RuntimeError('既有v171证据清单偏离2865300基线身份')
    proof = read_json(PROOF)
    inherited = {r['sha256']:r['archive'] for r in proof['files']+proof['finalSources']}
    inherited[PROOF_SHA256] = relative(PROOF)
    inherited_used = {}
    rows = collect()
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
    if snapshot(collect()) != before:
        raise RuntimeError('归档读取期间源文件集合或元数据变化，拒绝创建最终归档')
    terminal_guards(args.postprocessing_report)
    DEST.mkdir(parents=True, exist_ok=False)
    for target, value in pending:
        write_new(target, value)
    # 新归档逐字hash核验；读取不改源目录，未知失败均不清理。
    for record in small_records+final_records:
        if sha_file(ROOT/record['archive']) != record['sha256']:
            raise RuntimeError('归档逐字核验失败：'+record['archive'])
    if snapshot(collect()) != before or sha_file(PROOF) != PROOF_SHA256:
        raise RuntimeError('归档写入期间源文件集合/元数据/既有proof变化，保留已写证据并拒绝最终清单')
    report = dict(files=small_records, finalSources=final_records, localOnly=local_only,
        distinctArchivedFiles=len(pending), archivedBytes=sum(len(v) for _,v in pending),
        distinctReferencedSHA256=len(seen), inheritedProofDistinctReferences=len(inherited_used),
        inheritedProof=dict(path=relative(PROOF),sha256=PROOF_SHA256,baselineCommit='2865300',copied=False),
        inheritedArchiveReferences=[dict(sha256=h,archive=name) for h,name in sorted(inherited_used.items())],
        runtimeRecursiveInventory=False, runtimeReusedFromVersion=171,
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
    execute(parser.parse_args())
