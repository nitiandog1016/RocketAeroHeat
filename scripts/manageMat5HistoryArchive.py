"""固定二十份已退役结果的归档清单；保留当前依赖，删除交给逐路径PowerShell。"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys

from archiveMat5History import digest_file, pack, read_manifest
from compactRetiredOutputCaches import no_matlab, read, save
from planLegacyOutputCleanup import ROOT, OUT, rec, safe

DEST = ROOT / 'docs/validation/storage_20261031/output_review'
STORE = OUT / 'mat_history_store_20261031'
PREVIOUS = ROOT / 'docs/validation/storage_20261031/retired_lines/plan.json'


def prepare():
    """模块一：范围仅来自上轮已核验清单，逐文件确认原结果SHA和文件身份。"""
    no_matlab()
    assert not (DEST / 'plan.json').exists()
    previous = read(PREVIOUS)
    groups = []
    for group in previous['groups']:
        item = group['retained']
        source = safe(ROOT / item['path'])
        assert source.relative_to(OUT).parts[0] == 'retained_line_results_20261031'
        assert source.stat().st_nlink == 1, '不处理当前范围以外的共享对象。'
        assert source.stat().st_size == item['bytes'] and digest_file(source) == item['sha256']
        name = source.parent.name
        manifest = DEST / 'manifests' / ('v83_' + name + '.matmap.z')
        groups.append({'original': item, 'manifestPath': manifest.relative_to(ROOT).as_posix()})
    assert len(groups) == 20
    protected = {}
    # 当前状态/几何/方向等沿用上轮清单；另保留v94—v103全部数值文件。
    for entry in previous['protectedFiles']:
        file = safe(ROOT / entry['path'])
        if file.is_file():
            protected[entry['path']] = rec(file, entry['sha256'])
    for version in range(94, 104):
        folder = OUT / ('v%d_validation' % version)
        if folder.is_dir():
            for file in folder.rglob('*.mat'):
                protected[file.relative_to(ROOT).as_posix()] = rec(file)
    sources = [rec(ROOT / p) for p in (
        'scripts/archiveMat5History.py', 'scripts/manageMat5HistoryArchive.py',
        'scripts/auditMat5HistoryArchive.py', 'scripts/Invoke-Mat5HistoryArchiveCleanup.ps1',
        'tests/testMat5HistoryArchive.py', 'scripts/mat5NestedCache.py')]
    save(DEST / 'scope.json', {'groups': groups, 'protectedFiles': list(protected.values()),
                              'sourceFiles': sources, 'previousPlan': rec(PREVIOUS),
                              'freeBytesBeforeArchive': read(DEST / 'pilot.json')['freeBeforeArchive']})

    def one(group):
        source = safe(ROOT / group['original']['path'])
        manifest = safe(ROOT / group['manifestPath'])
        if manifest.exists():
            value = read_manifest(manifest)
            assert value['original']['sha256'] == group['original']['sha256']
            assert value['verification']['status'] == 'passed-exact-original-file-SHA'
        else:
            value = pack(source, STORE, manifest, group['original']['path'])
        return dict(group, manifest=rec(manifest), originalIdentity=value['original'],
                    verification=value['verification'])

    completed = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(one, g) for g in groups]):
            value = future.result()
            completed.append(value)
            save(DEST / 'progress.json', {'completed': len(completed), 'total': 20,
                                         'groups': completed, 'originalsDeleted': 0})
            print('ARCHIVED %d/20 %s' % (len(completed), value['original']['path']), flush=True)
    archive_files = [rec(p) for p in STORE.rglob('*.z')]
    manifests = [item['manifest'] for item in completed]
    original_bytes = sum(g['original']['bytes'] for g in completed)
    archive_bytes = sum(p['bytes'] for p in archive_files) + sum(p['bytes'] for p in manifests)
    assert archive_bytes < original_bytes, '不能减少占用，拒绝删除原文件。'
    plan = dict(read(DEST / 'scope.json'), groups=sorted(completed, key=lambda g: g['original']['path']),
                kind='lossless-retired-MAT-history-archive-20261031', archiveFiles=archive_files,
                createdUTC=datetime.now(timezone.utc).isoformat(), originalBytes=original_bytes,
                archiveBytes=archive_bytes, netContentBytes=original_bytes - archive_bytes,
                store=STORE.relative_to(ROOT).as_posix(),
                gitHead=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())
    save(DEST / 'plan.json', plan)
    print(json.dumps({k: plan[k] for k in ('originalBytes', 'archiveBytes', 'netContentBytes')}, ensure_ascii=False))


def postcheck():
    """模块二：确认其他历史文件身份及当前依赖不变；不把路径累计长度当磁盘实占。"""
    no_matlab()
    plan = read(DEST / 'plan.json')
    before = json.load(gzip.open(DEST / 'inventory_before.json.gz', 'rt', encoding='utf-8'))
    deleted = {g['original']['path'] for g in plan['groups']}
    unchanged = 0
    for relative, item in before['files'].items():
        file = safe(ROOT / relative)
        if relative in deleted:
            assert not file.exists()
            continue
        st = file.stat()
        assert (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns) == (
            item['device'], item['fileID'], item['bytes'], item['mtimeNS']), relative
        unchanged += 1
    for item in plan['protectedFiles'] + plan['sourceFiles'] + plan['archiveFiles']:
        assert digest_file(ROOT / item['path']) == item['sha256'], item['path']
    for group in plan['groups']:
        marker = read((ROOT / group['original']['path']).parent / 'EVIDENCE_ARCHIVED_20261031.json')
        assert marker['manifest'] == group['manifest']
        assert marker['originalSHA256'] == group['original']['sha256']
    report = {'status': 'passed', 'originalsRemoved': len(deleted), 'otherOriginalFilesUnchanged': unchanged,
              'protectedFilesVerified': len(plan['protectedFiles']), 'archiveObjects': len(plan['archiveFiles']),
              'netContentBytesReclaimed': plan['netContentBytes'],
              'freeBytesBeforeArchive': plan['freeBytesBeforeArchive'],
              'freeBytesAfterPostcheck': shutil.disk_usage(str(OUT)).free,
              'scientificDataDiscarded': False, 'currentMATLABRunChanged': False}
    report['freeIncreaseSinceArchiveStart'] = report['freeBytesAfterPostcheck'] - report['freeBytesBeforeArchive']
    save(DEST / 'postcheck.json', report)
    print(json.dumps(report, ensure_ascii=False))


if __name__ == '__main__':
    {'prepare': prepare, 'postcheck': postcheck}[sys.argv[1]]()
