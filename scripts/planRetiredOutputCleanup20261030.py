"""2026-10-30输出维护：核对重复中间文件和两份已退役QR构造，不删除独有结果。"""

from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import json
import subprocess
import sys

from scipy.io import loadmat, whosmat

from planLegacyOutputCleanup import ROOT, OUT, inventory, rec, safe, sha
from planRetiredQRCacheCleanup import equal
from planOutputStorageMaintenance import allocated

DEST = ROOT / 'docs/validation/storage_20261030'
CONTEXTS = [('v91_validation', 'interfacefull91'), ('v92_validation', 'repreviewfull92')]
RAW_FOLDERS = {
    'v70_validation/newton_step_comparison': ('newton_step_trial_raw.mat', 'newton_step_comparison.mat'),
    'v71_validation/rocket_failure_final': ('newton_diagnostic_raw.mat', 'newton_diagnostic_evidence.mat'),
    'v71_validation/rocket_failure_replay': ('newton_diagnostic_raw.mat', 'newton_diagnostic_evidence.mat'),
    'v72_validation/continuation_trial': ('continuation_raw.mat', 'continuation_evidence.mat'),
    'v72_validation/warm_pressure_trial': ('warm_pressure_raw.mat', 'warm_pressure_evidence.mat'),
    'v73_validation/pseudo_cfl_001': ('pseudo_cfl_raw.mat', 'pseudo_cfl_evidence.mat'),
    'v73_validation/pseudo_cfl_010': ('pseudo_cfl_raw.mat', 'pseudo_cfl_evidence.mat'),
    'v73_validation/pseudo_cfl_100': ('pseudo_cfl_raw.mat', 'pseudo_cfl_evidence.mat'),
    'v73_validation/fine_difference_1em09': ('fine_difference_raw.mat', 'fine_difference_evidence.mat'),
    'v73_validation/fine_difference_7em08': ('fine_difference_raw.mat', 'fine_difference_evidence.mat'),
    'v73_validation/warm_hlle_davis': ('warm_flux_raw.mat', 'warm_flux_evidence.mat'),
    'v73_validation/warm_hlle_einfeldt': ('warm_flux_raw.mat', 'warm_flux_evidence.mat'),
    'v74_validation/local_difference_2em06': ('local_difference_raw.mat', 'local_difference_evidence.mat'),
    'v74_validation/local_difference_2em07': ('local_difference_raw.mat', 'local_difference_evidence.mat'),
    'v75_validation/spatial_root_double': ('spatial_root_raw.mat', 'spatial_root_evidence.mat'),
    'v75_validation/spatial_root_two_double_euler': ('spatial_root_raw.mat', 'spatial_root_evidence.mat'),
    'v83_validation/localrocket02': ('local_coarse_construction_evidence.mat', 'local_coarse_evidence.mat'),
    'v83_validation/localunit02': ('local_coarse_construction_evidence.mat', 'local_coarse_evidence.mat'),
}


def read(path):
    return json.loads(safe(path).read_text(encoding='utf-8-sig'))


def plan():
    """模块一：只读逐变量比较，固定删除范围；所有差异文件保留。"""
    assert not DEST.exists(), '已有维护记录，不覆盖。'
    before = inventory()
    aliases = {}
    for relative, item in before['files'].items():
        aliases.setdefault((item['device'], item['fileID']), []).append(relative)
    protected, versions, archives, groups, rejected = {}, {}, {}, [], []
    historical = {}
    by_name = {}
    for relative in before['files']:
        if relative.endswith('.m'):
            by_name.setdefault(Path(relative).name, []).append(ROOT / relative)
    for base in (ROOT / 'docs/validation',):
        for file in base.rglob('*.raw'):
            historical.setdefault(sha(file), file)

    def protect(path, expected=None):
        item = rec(path, expected)
        protected[item['path']] = item
        return item

    def source(relative, expected):
        relative = relative.replace('\\', '/')
        key = (relative, expected)
        if key in versions:
            return
        current = safe(ROOT / relative)
        original = current if current.is_file() and sha(current) == expected else historical.get(expected)
        if original is not None:
            versions[key] = dict(rec(original, expected), restorePath=relative)
            return
        content = subprocess.run(['git', 'show', 'HEAD:' + relative], cwd=ROOT, capture_output=True)
        if content.returncode != 0 or hashlib.sha256(content.stdout).hexdigest() != expected:
            candidate = next((p for p in by_name.get(current.name, []) if sha(p) == expected), None)
            original_bytes = candidate.read_bytes() if candidate is not None else None
            if original_bytes is None:
                revisions = subprocess.check_output(['git', 'log', '--all', '--format=%H', '--', relative], cwd=ROOT, text=True).splitlines()
                for revision in revisions:
                    old = subprocess.run(['git', 'show', revision + ':' + relative], cwd=ROOT, capture_output=True)
                    if old.returncode == 0 and hashlib.sha256(old.stdout).hexdigest() == expected:
                        original_bytes = old.stdout
                        break
            assert original_bytes is not None, '缺少历史源码：' + relative
        else:
            original_bytes = content.stdout
        file = DEST / 'historical_sources' / (current.name + '.' + expected[:12] + '.raw')
        archives[file] = original_bytes
        versions[key] = dict(path=file.relative_to(ROOT).as_posix(), bytes=len(original_bytes), sha256=expected, restorePath=relative)

    def add(file, kind, **extra):
        s = file.stat()
        links = aliases[(s.st_dev, s.st_ino)]
        assert len(links) == s.st_nlink, '发现outputs外硬链接。'
        markers = sorted({(Path(p).parent / 'STORAGE_PRUNED_20261030.json').as_posix() for p in links})
        assert all(not (ROOT / p).exists() for p in markers)
        groups.append(dict(kind=kind, file=rec(file), aliases=[rec(ROOT / p, sha(file)) for p in links],
                           allocatedBytesBefore=allocated(file), markerPaths=markers, **extra))

    for folder, (old, final) in RAW_FOLDERS.items():
        old_file, final_file = OUT / folder / old, OUT / folder / final
        old_names = [v[0] for v in whosmat(old_file)]
        final_names = {v[0] for v in whosmat(final_file)}
        missing = set(old_names) - final_names
        different = []
        if not missing:
            # 不比较文件名或近似数值：旧MAT每个变量都必须在最终MAT中完全相等。
            a = loadmat(old_file, simplify_cells=True)
            b = loadmat(final_file, variable_names=old_names, simplify_cells=True)
            different = [name for name in old_names if not equal(a[name], b[name])]
            del a, b
        if missing or different:
            rejected.append(dict(path=old_file.relative_to(ROOT).as_posix(), missingVariables=sorted(missing), differentVariables=different))
            print('RETAIN_DIFFERENT ' + folder, flush=True)
            continue
        replacement = protect(final_file)
        add(old_file, 'all-old-MAT-variables-exactly-retained-in-final', replacement=replacement, variables=old_names)
        print('DUPLICATE_VERIFIED ' + folder, flush=True)

    # 模块二：仅两份旧QR因子可重建；温度行、方向、失败/前缀和原状态逐项保留。
    for version, job in CONTEXTS:
        file = OUT / version / 'continuous01/step01/direction_context.mat'
        assert {v[0] for v in whosmat(file)} == {'preparationFile', 'preparationSHA256', 'p', 'factor', 'initialRows'}
        m = loadmat(file, variable_names=['preparationFile', 'preparationSHA256', 'p', 'initialRows'], simplify_cells=True)
        prep = safe(m['preparationFile'])
        protect(prep, m['preparationSHA256'])
        p = read(prep)
        assert equal(p, m['p']), '保存准备清单与上下文不一致。'
        active = file.parent / 'active_direction.mat'
        d = loadmat(active, variable_names=['initialRows', 'preparationSHA256'], simplify_cells=True)
        assert d['preparationSHA256'] == m['preparationSHA256'] and equal(d['initialRows'], m['initialRows'])
        protect(active)
        protect(p['sourceBinding']['sourceFile'], p['sourceBinding']['sourceFileSHA256'])
        protect(p['geometryBinding']['file'], p['geometryBinding']['sha256'])
        protect(prep.parent / 'CACHE_PRUNED.json')
        process_file = OUT / version / (job + '_process.json')
        process = read(process_file)
        assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
        protect(process_file)
        driver_tests = {'tests/studyRocketTemperatureLookahead.m', 'tests/prepareRocketCurrentPeakSystem.m',
                        'tests/loadRocketPeakSystem.m', 'tests/createRocketNSDirectionProbe.m',
                        'tests/continueRocketWallVelocityKrylov.m'}
        for relative, digest in process['sourceSHA256'].items():
            # 全工程进程快照还列出未调用的测试；重建绑定核验计算模块及本入口的驱动依赖。
            if relative.endswith('.m') and (not relative.startswith('tests/') or relative in driver_tests):
                source(relative, digest)
        for path in file.parent.rglob('*'):
            if path.is_file() and path != file:
                protect(path)
        add(file, 'retired-rebuildable-QR-factor', preparationFile=rec(prep), rowsExactlyRetainedIn=rec(active),
            rebuildEntry='tests/studyRocketTemperatureLookahead.m', rebuildRequired=True)
        print('RETIRED_QR_VERIFIED ' + version, flush=True)

    # 模块三：锁定当前Q41和最新研究；计算源码、当前输入和独有结果均不修改。
    for relative in ('v93_validation/continuous01/step01/direction_context.mat',
                     'v93_validation/continuous01/step01/active_direction.mat',
                     'v93_validation/continuous01/step01/krylov01/wall_velocity_krylov_evidence.mat',
                     'v80_validation/domain02/shell_domain_evidence.mat'):
        protect(OUT / relative)
    for version in range(94, 99):
        for path in (OUT / ('v%d_validation' % version)).rglob('*'):
            if path.is_file():
                protect(path)
    for path in (ROOT / 'regionalFlow/shockLayer').glob('*.m'):
        protect(path)
    delete = [item for g in groups for item in g['aliases']]
    assert not set(protected).intersection(item['path'] for item in delete)
    DEST.mkdir(parents=True)
    for path, content in archives.items():
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(content)
    inv = DEST / 'inventory_before.json.gz'
    with inv.open('wb') as stream:
        with gzip.GzipFile(filename='', mode='wb', fileobj=stream, mtime=0) as zipped:
            zipped.write(json.dumps(before).encode('utf8'))
    data = dict(kind='verified-redundant-intermediates-and-two-retired-QR-factors', workspace=str(ROOT),
                createdUTC=datetime.now(timezone.utc).isoformat(), groups=groups, deleteFiles=delete,
                protectedFiles=list(protected.values()), sourceFiles=list(versions.values()),
                rejectedCandidates=rejected, inventoryBefore=rec(inv),
                plannedUniqueBytes=sum(g['file']['bytes'] for g in groups),
                plannedAllocatedBytes=sum(g['allocatedBytesBefore'] for g in groups),
                numericalRebuildExecuted=False)
    (DEST / 'plan.json').write_bytes((json.dumps(data, ensure_ascii=False, indent=2) + '\n').encode('utf8'))
    print(json.dumps(dict(objects=len(groups), paths=len(delete), bytes=data['plannedUniqueBytes'],
                          rejected=len(rejected), protected=len(protected), sources=len(versions))))


def audit():
    """模块四：重新遍历全outputs，核对只有计划路径消失，所有保留SHA及文件身份不变。"""
    data, receipt = read(DEST / 'plan.json'), read(DEST / 'execution.json')
    assert receipt['status'] == 'completed-verified-output-maintenance'
    assert receipt['planSHA256'] == sha(DEST / 'plan.json')
    rec(ROOT / data['inventoryBefore']['path'], data['inventoryBefore']['sha256'])
    with gzip.open(DEST / 'inventory_before.json.gz', 'rb') as stream:
        before = json.load(stream)
    after = inventory()
    removed = {v['path'] for v in data['deleteFiles']}
    markers = {p for g in data['groups'] for p in g['markerPaths']}
    assert set(before['files']) - set(after['files']) == removed
    assert set(after['files']) - set(before['files']) == markers
    assert all(before['files'][p] == after['files'][p] for p in before['files'] if p not in removed)
    for item in data['protectedFiles'] + data['sourceFiles']:
        rec(ROOT / item['path'], item['sha256'])
    report = dict(status='passed-full-output-inventory-and-preserved-SHA-audit',
                  deletedObjects=len(data['groups']), deletedPaths=len(removed),
                  unchangedOriginalFiles=len(before['files'])-len(removed),
                  protectedFiles=len(data['protectedFiles']), sourceVersions=len(data['sourceFiles']),
                  deletedUniqueBytes=data['plannedUniqueBytes'],
                  actualFreeIncreaseBytes=receipt['actualFreeIncreaseBytes'],
                  uniqueContentBytesAfter=after['uniqueFileLengthBytes'], currentQ41AndV94ThroughV98Retained=True,
                  numericalComputationExecuted=False)
    assert not (DEST / 'postcheck.json').exists()
    (DEST / 'postcheck.json').write_bytes((json.dumps(report, indent=2) + '\n').encode('utf8'))
    print(json.dumps(report))


def final_audit():
    """模块五：说明更新后再核对全目录，只允许已记录的README文档变化。"""
    data, update = read(DEST / 'plan.json'), read(DEST / 'readme_update.json')
    assert read(DEST / 'postcheck.json')['status'] == 'passed-full-output-inventory-and-preserved-SHA-audit'
    with gzip.open(DEST / 'inventory_before.json.gz', 'rb') as stream:
        before = json.load(stream)
    after = inventory()
    removed = {v['path'] for v in data['deleteFiles']}
    markers = {p for g in data['groups'] for p in g['markerPaths']}
    assert set(before['files']) - set(after['files']) == removed
    assert set(after['files']) - set(before['files']) == markers
    changed = {p for p in before['files'] if p not in removed and before['files'][p] != after['files'][p]}
    assert changed == {'outputs/README.md'}
    rec(ROOT / update['originalBytes']['path'], update['originalBytes']['sha256'])
    rec(ROOT / update['after']['path'], update['after']['sha256'])
    assert after['files']['outputs/README.md'] == {key: update['after'][key] for key in ('bytes', 'mtimeNS', 'fileID', 'device')}
    for item in data['protectedFiles'] + data['sourceFiles']:
        rec(ROOT / item['path'], item['sha256'])
    report = dict(status='passed-final-output-inventory-after-documented-README-update',
                  deletedPaths=len(removed), unchangedOriginalFiles=before['fileCount']-len(removed)-1,
                  onlyAdditionalDocumentationChange=sorted(changed), protectedFiles=len(data['protectedFiles']),
                  sourceVersions=len(data['sourceFiles']), uniqueContentBytesAfter=after['uniqueFileLengthBytes'],
                  currentQ41AndV94ThroughV98Retained=True, numericalComputationExecuted=False)
    assert not (DEST / 'final_postcheck.json').exists()
    (DEST / 'final_postcheck.json').write_bytes((json.dumps(report, indent=2) + '\n').encode('utf8'))
    print(json.dumps(report))


if __name__ == '__main__':
    {'plan': plan, 'audit': audit, 'final-audit': final_audit}[sys.argv[1]]()
