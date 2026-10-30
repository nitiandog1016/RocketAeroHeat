"""清理十轮已结束研究内嵌的预条件矩阵；保留全部科学字段原字节。"""

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

from scipy.io import loadmat, whosmat

from compactRetiredOutputCaches import no_matlab, read, save, snapshot
from mat5NestedCache import prune_result, same_science, scan_result
from planLegacyOutputCleanup import ROOT, OUT, rec, safe, sha

DEST = ROOT / 'docs/validation/storage_20261030/nested_operators'
RETAINED = OUT / 'retained_operator_results_20261030'
VERSIONS = (1, 2, 3, 4, 5, 6, 7, 8, 18, 19)
CACHE_FIELDS = {'diagonal', 'lower', 'upper', 'minusCoupling', 'plusCoupling',
                'secondLower', 'secondUpper', 'wlsNormalFactor'}


def records(path):
    values = []
    with safe(path).open('rb') as file:
        header = file.read(128)
        assert header.startswith(b'MATLAB 5.0 MAT-file,') and header[-2:] == b'IM'
        for name, shape, kind in whosmat(path):
            offset = file.tell()
            tag, size = struct.unpack('<II', file.read(8))
            assert tag == 15 and size > 0
            values.append(dict(name=name, shape=shape, kind=kind, offset=offset, bytes=size + 8))
            file.seek(size, 1)
        assert file.tell() == path.stat().st_size
    return header, values


def raw_hash(file, position, count, destination=None):
    file.seek(position)
    digest = hashlib.sha256()
    while count:
        block = file.read(min(count, 1024 * 1024))
        assert block
        digest.update(block)
        if destination is not None:
            destination.write(block)
        count -= len(block)
    return digest.hexdigest()


def migrate(original, retained):
    """模块一：逐记录复制；仅重写result.firstOperator的白名单构造字段。"""
    header, original_records = records(original)
    result = next(value for value in original_records if value['name'] == 'result')
    with original.open('rb') as file:
        scientific = scan_result(file, result['offset'] + 8, result['bytes'] - 8, CACHE_FIELDS)
    retained.parent.mkdir(parents=True, exist_ok=True)
    copied = []
    with original.open('rb') as old, safe(retained).open('xb') as new:
        new.write(header)
        for item in original_records:
            if item['name'] == 'result':
                changed = prune_result(old, new, item['offset'] + 8, item['bytes'] - 8, scientific)
            else:
                copied.append(dict(name=item['name'], originalOffset=item['offset'],
                                   retainedOffset=new.tell(), bytes=item['bytes'],
                                   sha256=raw_hash(old, item['offset'], item['bytes'], new)))
    # 新开文件再读，并由SciPy确认仍是合法的MAT5变量表。
    retained_header, retained_records = records(retained)
    assert header == retained_header
    assert [(r['name'], r['shape'], r['kind']) for r in original_records] == [
        (r['name'], r['shape'], r['kind']) for r in retained_records]
    result_after = next(value for value in retained_records if value['name'] == 'result')
    with retained.open('rb') as file:
        actual = scan_result(file, result_after['offset'] + 8, result_after['bytes'] - 8,
                             CACHE_FIELDS, require_removed=False)
    same_science(scientific, actual)
    with original.open('rb') as old, retained.open('rb') as new:
        for item in copied:
            assert raw_hash(old, item['originalOffset'], item['bytes']) == item['sha256']
            assert raw_hash(new, item['retainedOffset'], item['bytes']) == item['sha256']
    return dict(originalResultRecord=result, retainedResultRecord=result_after,
                resultScan=scientific, copiedRecords=copied,
                retainedScientificFields=len(scientific['kept']),
                removedRawMatrixBytes=sum(v['bytes'] for v in scientific['removed'].values()))


def plan():
    """模块二：固定十份原文件，绑定实际终态、几何和准确历史构造源码。"""
    no_matlab()
    assert not DEST.exists() and not RETAINED.exists(), '不覆盖旧清理回执或数值结果。'
    before = snapshot()
    DEST.mkdir(parents=True)
    with (DEST / 'inventory_before.json.gz').open('wb') as file:
        with gzip.GzipFile(filename='', mode='wb', fileobj=file, mtime=0) as compressed:
            compressed.write(json.dumps(before).encode('utf-8'))
    historical = {}
    for file in (ROOT / 'docs/validation').glob('storage_*/**/plan.json'):
        for item in read(file).get('sourceFiles', []):
            if (ROOT / item['path']).is_file():
                historical.setdefault(item['sha256'], ROOT / item['path'])
    output_sources = {}
    for name in before['files']:
        if name.endswith('.m'):
            output_sources.setdefault(Path(name).name, []).append(ROOT / name)
    sources, protected = {}, {}

    def protect(path, expected=None):
        item = rec(path, expected)
        protected[item['path']] = item
        return item

    def source(relative, expected):
        relative = relative.replace('\\', '/')
        key = (relative, expected)
        if key in sources:
            return
        current = safe(ROOT / relative)
        candidate = current if current.is_file() and sha(current) == expected else historical.get(expected)
        if candidate is None:
            candidate = next((p for p in output_sources.get(current.name, []) if sha(p) == expected), None)
        if candidate is None:
            for revision in subprocess.check_output(['git', 'log', '--all', '--format=%H', '--', relative],
                                                     cwd=ROOT, text=True).splitlines():
                old = subprocess.run(['git', 'show', revision + ':' + relative], cwd=ROOT, capture_output=True)
                if old.returncode == 0 and hashlib.sha256(old.stdout).hexdigest() == expected:
                    folder = DEST / 'historical_sources'
                    folder.mkdir(exist_ok=True)
                    candidate = folder / (current.name + '.' + expected[:16] + '.raw')
                    if not candidate.exists():
                        candidate.write_bytes(old.stdout)
                    break
        assert candidate is not None, '缺少准确历史源码：' + relative
        sources[key] = dict(rec(candidate, expected), restorePath=relative)

    for name in ('v83_validation/krylov41/wall_velocity_krylov_evidence.mat',
                 'v93_validation/continuous01/step01/active_direction.mat',
                 'v93_validation/continuous01/step01/direction_context.mat',
                 'v80_validation/domain02/shell_domain_evidence.mat',
                 'v101_validation/kernel02/shock_fitted_kernel_evidence.mat'):
        protect(OUT / name)
    # 真实几何、当前Q41、当前v101研究以及独有失败文件完全不改。
    for name, value in before['files'].items():
        if ('failure' in Path(name).name.lower() or 'checkpoint' in Path(name).name.lower()) and name.endswith('.mat'):
            protect(ROOT / name)
    groups = []
    for version in VERSIONS:
        tail = 'v83_validation/krylov%02d/wall_velocity_krylov_evidence.mat' % version
        original = safe(OUT / tail)
        original_stat = original.stat()
        assert original_stat.st_nlink == 1
        process_file = OUT / ('v83_validation/krylov%02d_process.json' % version)
        process = read(process_file)
        assert process['status'] == 'terminal' and process['exitCode'] == 0
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        protect(process_file)
        for relative, digest in process['sourcesBefore'].items():
            normalized = relative.replace('\\', '/')
            if normalized.endswith('.m') and (normalized.startswith('internalFlow/') or
                    Path(normalized).name in {'setupInternalFlow.m', 'studyRocketWallVelocityKrylov.m',
                                              'continueRocketWallVelocityKrylov.m'}):
                source(normalized, digest)
        assert any(Path(r).name == 'buildRocketShellLineOperator.m' for r in process['sourcesBefore'])
        report = loadmat(original, variable_names=['report'], simplify_cells=True)['report']
        assert report['engineProcessID'] in process['engineProcessIDs']
        retained = safe(RETAINED / tail)
        original_record = rec(original)
        migration = migrate(original, retained)
        assert original_stat.st_mtime_ns == original.stat().st_mtime_ns
        assert rec(original)['sha256'] == original_record['sha256']
        item = dict(original=original_record, retained=rec(retained), process=rec(process_file),
                    markerPath=(original.parent / 'OPERATOR_CACHE_PRUNED_20261030.json').relative_to(ROOT).as_posix(),
                    **migration)
        groups.append(item)
        print('VERIFIED krylov%02d science=%d removedRawMB=%.3f netMB=%.3f' % (
            version, item['retainedScientificFields'], item['removedRawMatrixBytes'] / 1e6,
            (item['original']['bytes'] - item['retained']['bytes']) / 1e6), flush=True)
    result = dict(kind='retired-nested-operator-cache-cleanup', status='planned-originals-not-deleted',
                  workspace=str(ROOT), createdUTC=datetime.now(timezone.utc).isoformat(),
                  gitHead=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  groups=groups, protectedFiles=list(protected.values()), sourceFiles=list(sources.values()),
                  originalBytes=sum(g['original']['bytes'] for g in groups),
                  retainedBytes=sum(g['retained']['bytes'] for g in groups),
                  scientificFields=sum(g['retainedScientificFields'] for g in groups),
                  note='只删预条件矩阵；其余每个科学矩阵的原MAT数据元素逐字节保留。新文件不是原容器，历史SHA不套用；需在新目录按准确源码重建缓存与来源。')
    result['netContentBytesReclaimed'] = result['originalBytes'] - result['retainedBytes']
    assert result['netContentBytesReclaimed'] > 0
    save(DEST / 'plan.json', result)
    print(json.dumps({k: result[k] for k in ('originalBytes', 'retainedBytes', 'netContentBytesReclaimed', 'scientificFields')}), flush=True)


def audit():
    """模块三：删除后独立读取保留MAT，复核全部科学字段、保留依赖和原目录身份。"""
    plan = read(DEST / 'plan.json')
    receipt = read(DEST / 'execution.json')
    assert receipt['status'] == 'completed-nested-operator-cache-cleanup'
    assert receipt['planSHA256'] == sha(DEST / 'plan.json')
    with gzip.open(DEST / 'inventory_before.json.gz', 'rt', encoding='utf-8') as file:
        before = json.load(file)
    deleted = set()
    for group in plan['groups']:
        assert not (ROOT / group['original']['path']).exists()
        deleted.add(group['original']['path'])
        file = ROOT / group['retained']['path']
        rec(file, group['retained']['sha256'])
        header, items = records(file)
        record = next(r for r in items if r['name'] == 'result')
        with file.open('rb') as stream:
            actual = scan_result(stream, record['offset'] + 8, record['bytes'] - 8,
                                 CACHE_FIELDS, require_removed=False)
            same_science(group['resultScan'], actual)
            for item in group['copiedRecords']:
                assert raw_hash(stream, item['retainedOffset'], item['bytes']) == item['sha256']
        marker = read(ROOT / group['markerPath'])
        assert marker['original'] == group['original'] and marker['retained'] == group['retained']
        print('AUDITED ' + group['original']['path'], flush=True)
    for item in plan['protectedFiles'] + plan['sourceFiles']:
        rec(ROOT / item['path'], item['sha256'])
    after = snapshot()
    unchanged = 0
    for name, value in before['files'].items():
        if name in deleted:
            assert name not in after['files']
        else:
            assert value == after['files'][name], name
            unchanged += 1
    result = dict(status='passed', auditedUTC=datetime.now(timezone.utc).isoformat(),
                  scientificFields=plan['scientificFields'], deletedOriginals=len(deleted),
                  otherOriginalFilesUnchanged=unchanged, protectedFiles=len(plan['protectedFiles']),
                  exactSourceVersions=len(plan['sourceFiles']),
                  netContentBytesReclaimed=plan['netContentBytesReclaimed'],
                  freeBytesBefore=before['freeBytes'], freeBytesAfter=after['freeBytes'],
                  freeIncreaseSinceBeforeMigration=after['freeBytes'] - before['freeBytes'],
                  compressedDataBytesAfter=after['compressedDataBytes'],
                  numericalSolvePerformed=False, largeOutputsUploaded=False)
    save(DEST / 'postcheck.json', result)
    print(json.dumps(result), flush=True)


def final_audit():
    """模块四：用途说明更新后只复核身份变化、归档源码与既有数值审计绑定。"""
    data = read(DEST / 'plan.json')
    receipt = read(DEST / 'execution.json')
    previous = read(DEST / 'postcheck.json')
    update = read(DEST / 'readme_update.json')
    assert previous['status'] == 'passed' and receipt['planSHA256'] == sha(DEST / 'plan.json')
    assert update['path'] == 'outputs/README.md'
    rec(ROOT / update['rawOriginal'], update['before']['sha256'])
    rec(ROOT / update['path'], update['after']['sha256'])
    with gzip.open(DEST / 'inventory_before.json.gz', 'rt', encoding='utf-8') as file:
        before = json.load(file)
    after = snapshot()
    deleted = {group['original']['path'] for group in data['groups']}
    unchanged = 0
    for name, value in before['files'].items():
        if name in deleted:
            assert name not in after['files']
        elif name == update['path']:
            assert value['fileID'] == update['before']['fileID']
            assert value['mtimeNS'] == update['before']['mtimeNS']
            assert value['bytes'] == update['before']['bytes']
        else:
            assert value == after['files'][name], name
            unchanged += 1
    for item in data['protectedFiles'] + data['sourceFiles'] + [g['retained'] for g in data['groups']]:
        rec(ROOT / item['path'], item['sha256'])
    archived = read(DEST / 'source_archive.json')
    for item in archived['files']:
        rec(ROOT / item['path'], item['sha256'])
    tested = read(DEST / 'tests.json')
    assert tested['status'] == 'passed' and tested['exitCode'] == 0
    for item in tested['sources']:
        rec(ROOT / item['path'], item['sha256'])
    result = dict(status='passed-after-documented-README-update',
                  auditedUTC=datetime.now(timezone.utc).isoformat(),
                  otherOriginalFilesUnchanged=unchanged, documentedChangedFile=update['path'],
                  scientificAudit=rec(DEST / 'postcheck.json'),
                  planSHA256=sha(DEST / 'plan.json'), sourceArchiveFiles=len(archived['files']),
                  currentSources=[rec(ROOT / name) for name in ('scripts/mat5NestedCache.py',
                                  'scripts/compactNestedOutputOperators.py',
                                  'scripts/Invoke-NestedOperatorCleanup.ps1', 'tests/testMat5NestedCache.py')],
                  compressedDataBytesAfter=after['compressedDataBytes'],
                  netContentBytesReclaimed=data['netContentBytesReclaimed'], largeOutputsUploaded=False)
    save(DEST / 'final_postcheck.json', result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    {'plan': plan, 'audit': audit, 'final-audit': final_audit}[sys.argv[1]]()
