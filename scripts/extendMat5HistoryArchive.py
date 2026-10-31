"""追加归档十八份退役MAT；共享旧对象库，完整保留数值并核实实际空间收益。"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys

from archiveMat5History import digest_file, pack, read_manifest
from auditMat5HistoryArchive import audit_manifest
from compactRetiredOutputCaches import no_matlab, read, save, snapshot
from planLegacyOutputCleanup import ROOT, OUT, rec, safe

DEST = ROOT / 'docs/validation/storage_20261031/result_deduplication'
STORE = OUT / 'mat_history_store_20261031'
OLD_ARCHIVE = ROOT / 'docs/validation/storage_20261031/output_review/plan.json'
INPUT_PLANS = (
    ('docs/validation/storage_20261030/nested_operators/plan.json',
     'retained_operator_results_20261030', 10),
    ('docs/validation/storage_20261031/graph_operators/plan.json',
     'retained_graph_results_20261031', 8),
)


def selected_result(relative):
    """模块一：只接受两组已清理矩阵的历史结果，拒绝当前状态和其他目录。"""
    path = safe(ROOT / relative)
    parts = path.relative_to(OUT).parts
    allowed = {entry[1] for entry in INPUT_PLANS}
    assert len(parts) == 4 and parts[0] in allowed and parts[1] == 'v83_validation'
    assert parts[2] in {'krylov%02d' % i for i in list(range(1, 16)) + [17, 18, 19]}
    assert parts[3] == 'wall_velocity_krylov_evidence.mat'
    return path


def object_delta(before, after):
    """模块二：已有共享对象不重复计为本轮新增；拒绝旧对象消失或被覆盖。"""
    for name, old in before.items():
        assert name in after and after[name] == old, '已有归档对象变化：' + name
    return sorted(set(after) - set(before))


def file_identity(path):
    value = safe(path).stat()
    return dict(bytes=value.st_size, mtimeNS=value.st_mtime_ns,
                device=value.st_dev, fileID=value.st_ino)


def store_identities():
    return {p.relative_to(ROOT).as_posix(): file_identity(p) for p in STORE.rglob('*.z')}


def prepare():
    """模块三：绑定旧清单、当前依赖及原文件；先归档和恢复审核，暂不删除。"""
    no_matlab()
    assert not DEST.exists(), '已有本轮记录，不能覆盖或重复准备。'
    assert STORE.is_dir() and (STORE / 'objects').is_dir()
    DEST.mkdir(parents=True)
    before = snapshot()
    with gzip.GzipFile(filename=str(DEST / 'inventory_before.json.gz'), mode='wb', mtime=0) as stream:
        stream.write(json.dumps(before).encode('utf-8'))
    print('INVENTORY %d paths; %.3f GB allocated data' % (
        before['fileCount'], before['compressedDataBytes'] / 1e9), flush=True)
    initial_objects = store_identities()
    previous = read(OLD_ARCHIVE)
    assert set(initial_objects) == {item['path'] for item in previous['archiveFiles']}
    groups, history_plans = [], [rec(OLD_ARCHIVE)]
    for plan_path, allowed, count in INPUT_PLANS:
        history_plans.append(rec(ROOT / plan_path))
        old = read(ROOT / plan_path)
        assert len(old['groups']) == count
        for item in old['groups']:
            original = item['retained']
            file = selected_result(original['path'])
            assert file.relative_to(OUT).parts[0] == allowed
            assert file.stat().st_nlink == 1, '拒绝范围外共享硬链接。'
            original = rec(file, original['sha256'])
            manifest = DEST / 'manifests' / ('v83_' + file.parent.name + '.matmap.z')
            groups.append(dict(original=original, manifestPath=manifest.relative_to(ROOT).as_posix()))
    assert len(groups) == 18 and len({g['manifestPath'] for g in groups}) == 18
    selected = {g['original']['path'] for g in groups}
    protected = {}
    for entry in previous['protectedFiles']:
        path = safe(ROOT / entry['path'])
        assert path.is_file(), '当前保留依赖缺失。'
        protected[entry['path']] = rec(path, entry['sha256'])
    for version in range(94, 104):
        for path in (OUT / ('v%d_validation' % version)).rglob('*.mat'):
            protected[path.relative_to(ROOT).as_posix()] = rec(path)
    assert not selected.intersection(protected), '归档范围涉及当前依赖。'
    source_files = [rec(ROOT / rel) for rel in (
        'scripts/extendMat5HistoryArchive.py', 'scripts/Invoke-Mat5HistoryArchiveExtension.ps1',
        'scripts/archiveMat5History.py', 'scripts/auditMat5HistoryArchive.py',
        'scripts/mat5NestedCache.py', 'scripts/compactRetiredOutputCaches.py',
        'scripts/planLegacyOutputCleanup.py', 'tests/testMat5HistoryArchive.py',
        'tests/testMat5HistoryArchiveExtension.py')]
    scope = dict(groups=groups, protectedFiles=list(protected.values()), sourceFiles=source_files,
                 historicalPlans=history_plans, freeBytesBeforeArchive=before['freeBytes'],
                 initialStoreFiles=initial_objects, currentNumericalRunChanged=False)
    save(DEST / 'scope.json', scope)

    def one(group):
        manifest = safe(ROOT / group['manifestPath'])
        value = pack(selected_result(group['original']['path']), STORE, manifest, group['original']['path'])
        return dict(group, manifest=rec(manifest), originalIdentity=value['original'],
                    verification=value['verification'])

    completed = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(one, group) for group in groups]):
            completed.append(future.result())
            save(DEST / 'progress.json', dict(completed=len(completed), total=18,
                                            groups=completed, originalsDeleted=0))
            print('ARCHIVED %d/18 %s' % (len(completed), completed[-1]['original']['path']), flush=True)
    current_objects = store_identities()
    added = object_delta(initial_objects, current_objects)
    archive_files = [rec(ROOT / relative) for relative in sorted(current_objects)]
    added_objects = [item for item in archive_files if item['path'] in set(added)]
    original_bytes = sum(g['original']['bytes'] for g in completed)
    added_bytes = sum(item['bytes'] for item in added_objects)
    manifest_bytes = sum(g['manifest']['bytes'] for g in completed)
    assert added_bytes + manifest_bytes < original_bytes, '无净存储收益，不允许删除。'
    plan = dict(scope, groups=sorted(completed, key=lambda g: g['original']['path']),
                kind='lossless-retired-MAT-shared-store-extension-20261031',
                store=STORE.relative_to(ROOT).as_posix(), archiveFiles=archive_files,
                addedArchiveFiles=added_objects, originalBytes=original_bytes,
                addedObjectBytes=added_bytes, manifestBytes=manifest_bytes,
                netContentBytes=original_bytes - added_bytes - manifest_bytes,
                createdUTC=datetime.now(timezone.utc).isoformat(),
                gitHead=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())
    save(DEST / 'plan.json', plan)
    print(json.dumps({k: plan[k] for k in (
        'originalBytes', 'addedObjectBytes', 'manifestBytes', 'netContentBytes')}), flush=True)


def audit():
    """模块四：独立实现另进程重构全部原MAT的SHA，逐份公布真实审核进度。"""
    no_matlab()
    plan = read(DEST / 'plan.json')
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(audit_manifest, ROOT / g['manifest']['path'], ROOT / plan['store'])
                   for g in plan['groups']]
        for future in as_completed(futures):
            results.append(future.result())
            save(DEST / 'independent_audit_progress.json', dict(verified=len(results), total=18, results=results))
            print('VERIFIED %d/18 %s' % (len(results), results[-1]['original']), flush=True)
    expected = {g['original']['path']: (g['original']['bytes'], g['original']['sha256']) for g in plan['groups']}
    actual = {r['original']: (r['bytes'], r['sha256']) for r in results}
    assert actual == expected and len(results) == 18
    report = dict(status='passed-independent-exact-original-MAT-SHA', verifiedOriginals=len(results),
                  verifiedOriginalBytes=sum(r['bytes'] for r in results), results=results,
                  planSHA256=digest_file(DEST / 'plan.json'), archiveReplayFunctionUsed=False,
                  auditor=rec(ROOT / 'scripts/auditMat5HistoryArchive.py'),
                  runner=rec(Path(__file__)), freeBytes=shutil.disk_usage(str(OUT)).free)
    save(DEST / 'predelete_archive_audit.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}), flush=True)


def postcheck():
    """模块五：逐路径核对全部其他原文件及共享对象，报告清理后的真实占用。"""
    no_matlab()
    plan = read(DEST / 'plan.json')
    before = json.load(gzip.open(DEST / 'inventory_before.json.gz', 'rt', encoding='utf-8'))
    deleted = {g['original']['path'] for g in plan['groups']}
    readme_update = DEST / 'readme_update.json'
    readme = read(readme_update) if readme_update.exists() else None
    unchanged = 0
    for relative, item in before['files'].items():
        file = safe(ROOT / relative)
        if relative in deleted:
            assert not file.exists(), relative
        elif readme is not None and relative == 'outputs/README.md':
            rec(file, readme['after']['sha256'])
            rec(ROOT / readme['beforeCopy']['path'], readme['beforeCopy']['sha256'])
        else:
            assert file_identity(file) == {k: item[k] for k in ('bytes', 'mtimeNS', 'device', 'fileID')}, relative
            unchanged += 1
    for item in plan['protectedFiles'] + plan['sourceFiles'] + plan['historicalPlans'] + plan['archiveFiles']:
        rec(ROOT / item['path'], item['sha256'])
    for group in plan['groups']:
        rec(ROOT / group['manifest']['path'], group['manifest']['sha256'])
        marker = read((ROOT / group['original']['path']).parent / 'EVIDENCE_ARCHIVED_20261031.json')
        assert marker['manifest'] == group['manifest']
        assert marker['originalSHA256'] == group['original']['sha256']
        assert read_manifest(ROOT / group['manifest']['path'])['original']['sha256'] == group['original']['sha256']
    after = snapshot()
    with gzip.GzipFile(filename=str(DEST / 'inventory_after.json.gz'), mode='wb', mtime=0) as stream:
        stream.write(json.dumps(after).encode('utf-8'))
    storage = {k: v for k, v in after.items() if k != 'files'}
    save(DEST / 'storage_after.json', storage)
    report = dict(status='passed', originalsRemoved=len(deleted), otherOriginalFilesUnchanged=unchanged,
                  protectedFilesVerified=len(plan['protectedFiles']), archiveObjects=len(plan['archiveFiles']),
                  addedObjects=len(plan['addedArchiveFiles']), netContentBytesReclaimed=plan['netContentBytes'],
                  freeBytesBeforeArchive=plan['freeBytesBeforeArchive'], freeBytesAfterPostcheck=after['freeBytes'],
                  freeIncreaseSinceArchiveStart=after['freeBytes'] - plan['freeBytesBeforeArchive'],
                  storage=storage, scientificDataDiscarded=False, currentNumericalRunChanged=False,
                  readmeUpdated=readme is not None)
    save(DEST / ('final_postcheck.json' if readme else 'postcheck.json'), report)
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    {'prepare': prepare, 'audit': audit, 'postcheck': postcheck}[sys.argv[1]]()
