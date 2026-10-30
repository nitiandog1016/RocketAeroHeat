"""回收八份旧研究的图预条件矩阵；沿用已核验的流式MAT迁存。"""

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import gzip
import json
import subprocess
import sys

import compactNestedOutputOperators as core
from compactRetiredOutputCaches import no_matlab, read, save
from planLegacyOutputCleanup import ROOT, OUT, rec, safe, sha

DEST = ROOT / 'docs/validation/storage_20261031/graph_operators'
RETAINED = OUT / 'retained_graph_results_20261031'
VERSIONS = (9, 10, 11, 12, 13, 14, 15, 17)
CACHE_FIELDS = core.CACHE_FIELDS | {'lineInverse', 'lineUpperFactor'}
PROCESS_NAMES = dict(zip(VERSIONS, ('coarsetrial01', 'layertrial01', 'targettrial01',
    'mixedtrial01', 'localtrial01', 'guesstrial01', 'recycledtrial01', 'wlsnormaltrial02')))


def configure():
    """模块一：固定本轮目录与八份白名单，禁止覆盖上一轮回执。"""
    core.DEST = DEST
    core.RETAINED = RETAINED
    core.VERSIONS = VERSIONS
    core.CACHE_FIELDS = CACHE_FIELDS


def plan():
    """模块二：先完整迁存及逐字段核对；此步骤不删除原文件。"""
    configure()
    no_matlab()
    assert not (DEST / 'plan.json').exists() and not RETAINED.exists()
    DEST.mkdir(parents=True, exist_ok=True)
    previous = DEST / 'inventory_before.json.gz'
    if previous.exists():
        # 首稿误按krylov编号寻找进程回执，在迁存/删除之前拒绝；原清单保留。
        folder = DEST / 'preflight'; folder.mkdir(exist_ok=True)
        safe(previous).rename(safe(folder / 'inventory_attempt01.json.gz'))
    before = core.snapshot()
    with previous.open('wb') as file:
        with gzip.GzipFile(filename='', mode='wb', fileobj=file, mtime=0) as compressed:
            compressed.write(json.dumps(before).encode('utf-8'))
    historical = {}
    for file in (ROOT / 'docs/validation').glob('storage_*/**/plan.json'):
        for item in read(file).get('sourceFiles', []):
            if (ROOT / item['path']).is_file():
                historical.setdefault(item['sha256'], ROOT / item['path'])
    output_sources = {}
    for relative in before['files']:
        if relative.endswith('.m'):
            output_sources.setdefault(Path(relative).name, []).append(ROOT / relative)
    sources, protected = {}, {}

    def protect(path, expected=None):
        item = rec(path, expected); protected[item['path']] = item; return item

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
                old = subprocess.run(['git', 'show', revision + ':' + relative], cwd=ROOT,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if old.returncode == 0 and hashlib.sha256(old.stdout).hexdigest() == expected:
                    folder = DEST / 'historical_sources'; folder.mkdir(exist_ok=True)
                    candidate = folder / (current.name + '.' + expected[:16] + '.raw')
                    if not candidate.exists():
                        candidate.write_bytes(old.stdout)
                    break
        assert candidate is not None, 'Missing accurate source: ' + relative
        sources[key] = dict(rec(candidate, expected), restorePath=relative)

    for relative in ('v83_validation/krylov41/wall_velocity_krylov_evidence.mat',
                     'v93_validation/continuous01/step01/active_direction.mat',
                     'v93_validation/continuous01/step01/direction_context.mat',
                     'v80_validation/domain02/shell_domain_evidence.mat',
                     'v101_validation/kernel02/shock_fitted_kernel_evidence.mat'):
        protect(OUT / relative)
    for relative in before['files']:
        if relative.endswith('.mat') and ('failure' in Path(relative).name.lower() or
                                         'checkpoint' in Path(relative).name.lower()):
            protect(ROOT / relative)
    groups = []
    for version in VERSIONS:
        tail = 'v83_validation/krylov%02d/wall_velocity_krylov_evidence.mat' % version
        original = safe(OUT / tail)
        original_stat = original.stat(); assert original_stat.st_nlink == 1
        process_file = OUT / ('v83_validation/' + PROCESS_NAMES[version] + '_process.json')
        process = read(process_file)
        expected_exit = 1 if version == 15 else 0
        assert process['status'] == 'terminal' and process['exitCode'] == expected_exit
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        protect(process_file)
        # k15数值已保存后旧UI导出失败；原失败及单独成功重导出均保留，不能冒称原宿主成功。
        if version == 15:
            recovery = OUT / 'v83_validation/recycledexport01_process.json'
            recovered = read(recovery)
            assert recovered['status'] == 'terminal' and recovered['exitCode'] == 0
            assert all(recovered['sourceStableAfter'].values())
            protect(recovery)
        for relative, digest in process['sourcesBefore'].items():
            normalized = relative.replace('\\', '/')
            if normalized.endswith('.m') and (normalized.startswith('internalFlow/') or
                    Path(normalized).name in {'setupInternalFlow.m', 'studyRocketWallVelocityKrylov.m',
                                              'continueRocketWallVelocityKrylov.m'}):
                source(normalized, digest)
        assert any(Path(r).name == 'buildRocketShellGraphPreconditioner.m' for r in process['sourcesBefore'])
        report = core.loadmat(original, variable_names=['report'], simplify_cells=True)['report']
        assert report['engineProcessID'] in process['engineProcessIDs']
        retained = safe(RETAINED / tail)
        original_record = rec(original)
        migration = core.migrate(original, retained)
        assert original.stat().st_mtime_ns == original_stat.st_mtime_ns
        assert rec(original)['sha256'] == original_record['sha256']
        marker = (original.parent / 'OPERATOR_CACHE_PRUNED_20261031.json').relative_to(ROOT).as_posix()
        assert not (ROOT / marker).exists()
        group = dict(original=original_record, retained=rec(retained), process=rec(process_file),
                     originalProcessExitCode=expected_exit, markerPath=marker, **migration)
        groups.append(group)
        print('VERIFIED krylov%02d science=%d removedRawMB=%.3f netMB=%.3f' % (
            version, group['retainedScientificFields'], group['removedRawMatrixBytes']/1e6,
            (group['original']['bytes']-group['retained']['bytes'])/1e6), flush=True)
    # 最新v102研究尚未交付，本轮只核对并保留，绝不混入维护提交。
    for file in (OUT / 'v102_validation').rglob('*.mat'):
        protect(file)
    data = dict(kind='retired-nested-operator-cache-cleanup', scope='retired-v83-graph-operators-20261031',
                status='planned-originals-not-deleted', workspace=str(ROOT),
                createdUTC=datetime.now(timezone.utc).isoformat(),
                gitHead=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                groups=groups, protectedFiles=list(protected.values()), sourceFiles=list(sources.values()),
                originalBytes=sum(g['original']['bytes'] for g in groups),
                retainedBytes=sum(g['retained']['bytes'] for g in groups),
                scientificFields=sum(g['retainedScientificFields'] for g in groups),
                cacheFieldsAllowed=sorted(CACHE_FIELDS), numericalSolvePerformed=False)
    data['netContentBytesReclaimed'] = data['originalBytes']-data['retainedBytes']
    assert data['netContentBytesReclaimed'] > 0
    data['note'] = ('仅移除八份旧结果firstOperator中的图/法向线矩阵、线分解和WLS分解；'
                    '全部其余科学字段的原MAT元素保留。旧审计须在新目录重建缓存与来源；'
                    '当前Q41、v101最终结果和v102全部MAT保留。')
    save(DEST / 'plan.json', data)
    print('PLAN_READY files=%d netBytes=%d' % (len(data['groups']), data['netContentBytesReclaimed']), flush=True)


def archive_sources():
    """模块三：只存缺少Git原字节的小型历史源码，避免依赖唯一outputs副本。"""
    data = read(DEST / 'plan.json')
    commit = data['gitHead']
    archived = []
    bindings = []
    folder = DEST / 'historical_sources'
    for item in data['sourceFiles']:
        path = safe(ROOT / item['path'])
        rec(path, item['sha256'])
        if item['path'].startswith('docs/'):
            bindings.append(dict(item, preservedBy='existing-source-archive'))
            continue
        result = subprocess.run(['git', 'show', commit + ':' + item['path']], cwd=ROOT,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == item['sha256']:
            bindings.append(dict(item, preservedBy='git', commit=commit))
            continue
        folder.mkdir(exist_ok=True)
        target = folder / (path.name + '.' + item['sha256'][:16] + '.raw')
        if not target.exists():
            target.write_bytes(path.read_bytes())
        saved = dict(rec(target, item['sha256']), restorePath=item['restorePath'], originalPath=item['path'])
        archived.append(saved)
        bindings.append(dict(item, preservedBy='new-source-archive', archivePath=saved['path']))
    # 只更新来源保留位置，原始计算身份及SHA不变。
    by_hash = {item['sha256']: item for item in archived}
    for item in data['sourceFiles']:
        if item['sha256'] in by_hash:
            item['originalPath'] = item['path']
            item['path'] = by_hash[item['sha256']]['path']
    save(DEST / 'plan.json', data)
    save(DEST / 'source_archive.json', dict(status='passed', files=archived,
         sourceBindings=bindings, exactSourceVersions=len(data['sourceFiles']),
         createdUTC=datetime.now(timezone.utc).isoformat()))
    print('SOURCES_VERIFIED versions=%d newArchives=%d' % (len(data['sourceFiles']), len(archived)), flush=True)


def audit():
    """模块四：独立重读迁存MAT、来源和全目录身份，不复跑物理求解。"""
    configure()
    core.audit()
    data = read(DEST / 'plan.json')
    archive = read(DEST / 'source_archive.json')
    for item in archive['files']:
        rec(ROOT / item['path'], item['sha256'])
    check = read(DEST / 'postcheck.json')
    check['sourceArchiveFiles'] = len(archive['files'])
    check['scope'] = data['scope']
    save(DEST / 'postcheck.json', check)


def final_audit():
    """模块五：用途文档更新后核对唯一允许的README变更及保留结果。"""
    configure()
    core.final_audit()
    result = read(DEST / 'final_postcheck.json')
    result['currentMaintenanceSources'] = [rec(ROOT / name) for name in (
        'scripts/compactRetiredGraphOperators.py', 'scripts/Invoke-RetiredGraphOperatorCleanup.ps1',
        'tests/testRetiredGraphCache.py')]
    save(DEST / 'final_postcheck.json', result)


if __name__ == '__main__':
    {'plan': plan, 'archive-sources': archive_sources, 'audit': audit,
     'final-audit': final_audit}[sys.argv[1]]()
