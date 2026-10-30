"""回收旧MAT内嵌构造缓存；其余变量的原压缩记录逐字节迁存，不加载大型矩阵。"""

import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys

from scipy.io import loadmat, whosmat
from planLegacyOutputCleanup import ROOT, OUT, rec, safe, sha

DEST = ROOT / 'docs/validation/storage_20261030/embedded_caches'
RETAINED = OUT / 'retained_evidence_20261030'
# 固定白名单：只移除可由原Q、几何和精确历史源码重建的构造变量。
TARGETS = (
    ('coarserocket01/volume_coarse_real_evidence.mat', {'coarse'}, 'coarserocket01', 'studyRocketVolumeCoarse'),
    ('coarselsrocket01/least_squares_coarse_evidence.mat', {'coarse'}, 'coarsels01', 'studyRocketLeastSquaresCoarse'),
    ('layercoarserocket01/coarse_layer_evidence.mat', {'coarse'}, 'layercoarse01', 'studyRocketCoarseLayerEnrichment'),
    ('basiscandidatesrocket01/primitive/coarse_candidate_evidence.mat', {'coarse'}, 'basiscandidates01', 'studyRocketCoarseBasisCandidates'),
    ('basiscandidatesrocket01/spatial/coarse_candidate_evidence.mat', {'coarse'}, 'basiscandidates01', 'studyRocketCoarseBasisCandidates'),
    ('mixedrocket01/mixed_coarse_evidence.mat', {'coarse'}, 'mixedstudy01', 'studyRocketMixedCoarse'),
    ('mixedrocket02/mixed_coarse_evidence.mat', {'coarse'}, 'mixedstudy02', 'studyRocketMixedCoarse'),
    ('localrocket02/local_coarse_evidence.mat', {'coarse'}, 'localstudy03', 'studyRocketLocalCoarse'),
    ('reservoirjacreal01/riemann_selected_direction_real_evidence.mat', {'cache'}, 'reservoirjacreal01', 'studyRiemannReservoirSelectedDirection'),
    ('schurcheckpointreal01/checkpoint/schur_assembly_evidence.mat',
     {'momentumMatrix', 'schurMatrix', 'couplingB', 'couplingE'}, 'schurcheckpointreal01', 'studyRocketSchurFailureCheckpoint'),
    ('schurcheckpointreal01/checkpoint/momentum_factor_evidence.mat', {'M'}, 'schurcheckpointreal01', 'studyRocketSchurFailureCheckpoint'),
)
DRIVERS = {entry for _, _, _, entry in TARGETS} | {
    'studyRocketSteadyGraph', 'studyRocketBlockILU', 'studyRocketBlockILU0',
    'studyRocketWallVelocityKrylov', 'continueRocketWallVelocityKrylov',
}


def save(path, value):
    safe(path).write_bytes((json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))


def read(path):
    return json.loads(safe(path).read_text(encoding='utf-8-sig'))


def no_matlab():
    state = subprocess.check_output(['tasklist', '/FI', 'IMAGENAME eq MATLAB.exe', '/FO', 'CSV'], text=True)
    assert '"MATLAB.exe"' not in state, 'MATLAB运行期间不能维护计算数据。'


def snapshot():
    """模块一：逐路径记录身份，按NTFS文件ID去重并计入透明压缩。"""
    function = ctypes.WinDLL('kernel32', use_last_error=True).GetCompressedFileSizeW
    function.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
    function.restype = wintypes.DWORD
    files, unique = {}, {}
    for base, dirs, names in os.walk(str(OUT), followlinks=False):
        for name in dirs + names:
            p = Path(base) / name
            assert not (getattr(p.lstat(), 'st_file_attributes', 0) & 0x400), '拒绝重解析路径。'
        for name in names:
            p = Path(base) / name
            s = p.stat()
            value = dict(bytes=s.st_size, mtimeNS=s.st_mtime_ns, fileID=s.st_ino, device=s.st_dev, links=s.st_nlink)
            files[p.relative_to(ROOT).as_posix()] = value
            key = (s.st_dev, s.st_ino)
            if key not in unique:
                high = wintypes.DWORD(0)
                ctypes.set_last_error(0)
                low = function(str(p), ctypes.byref(high))
                assert low != 0xFFFFFFFF or ctypes.get_last_error() == 0
                unique[key] = (s.st_size, (high.value << 32) | low)
    return dict(files=files, fileCount=len(files), pathBytes=sum(v['bytes'] for v in files.values()),
                uniqueFileCount=len(unique), uniqueContentBytes=sum(v[0] for v in unique.values()),
                compressedDataBytes=sum(v[1] for v in unique.values()), freeBytes=shutil.disk_usage(str(OUT)).free)


def copy_records(source, target, removed):
    """模块二：保留MAT头与每个非缓存变量的原压缩字节，不重新编码科学数据。"""
    names = whosmat(source)
    assert len({v[0] for v in names}) == len(names) and removed < {v[0] for v in names}
    target.parent.mkdir(parents=True, exist_ok=True)
    records = []
    with safe(source).open('rb') as old, safe(target).open('xb') as new:
        header = old.read(128)
        assert header.startswith(b'MATLAB 5.0 MAT-file,') and header[-2:] == b'IM'
        new.write(header)
        for name, shape, kind in names:
            offset = old.tell()
            tag = old.read(8)
            record_kind, size = struct.unpack('<II', tag)
            assert record_kind == 15 and size > 0, '仅接受已核实的独立压缩MAT5记录。'
            keep = name not in removed
            retained_offset = new.tell() if keep else None
            if keep:
                new.write(tag)
            digest = hashlib.sha256(tag)
            remaining = size
            while remaining:
                block = old.read(min(remaining, 4 * 1024**2))
                assert block
                digest.update(block)
                if keep:
                    new.write(block)
                remaining -= len(block)
            records.append(dict(variable=name, shape=shape, matlabType=kind, originalOffset=offset,
                                retainedOffset=retained_offset, bytes=size + 8,
                                rawRecordSHA256=digest.hexdigest(), retained=keep))
        assert old.tell() == source.stat().st_size
    assert whosmat(target) == [v for v in names if v[0] not in removed]
    # 第二次独立打开逐块比较原字节，避免只依靠复制过程中的SHA。
    with source.open('rb') as old, target.open('rb') as new:
        assert old.read(128) == new.read(128)
        for record in records:
            if not record['retained']:
                continue
            old.seek(record['originalOffset'])
            assert new.tell() == record['retainedOffset']
            remaining = record['bytes']
            while remaining:
                count = min(remaining, 4 * 1024**2)
                assert old.read(count) == new.read(count)
                remaining -= count
        assert not new.read(1)
    return records


def plan():
    """模块三：固定范围、保留状态及准确历史源码；此阶段不删除原文件。"""
    no_matlab()
    assert not DEST.exists() and not RETAINED.exists(), '禁止覆盖维护记录或既有迁存数据。'
    before = snapshot()
    DEST.mkdir(parents=True)
    with gzip.GzipFile(filename='', mode='wb', fileobj=(DEST/'inventory_before.json.gz').open('wb'), mtime=0) as stream:
        stream.write(json.dumps(before, ensure_ascii=False).encode('utf-8'))
    all_targets = {(OUT/'v83_validation'/name).relative_to(ROOT).as_posix() for name, _, _, _ in TARGETS}
    historical, output_sources = {}, {}
    for file in (ROOT/'docs/validation').rglob('*.raw'):
        historical.setdefault(sha(file), file)
    for name in before['files']:
        if name.endswith('.m'):
            output_sources.setdefault(Path(name).name, []).append(ROOT/name)
    protected, versions, groups, processes = {}, {}, [], {}

    def protect(file, expected=None):
        item = rec(file, expected)
        assert item['path'] not in all_targets
        protected[item['path']] = item
        return item

    def source(relative, expected):
        relative = relative.replace('\\', '/')
        key = (relative, expected)
        if key in versions:
            return
        current = safe(ROOT/relative)
        candidate = current if current.is_file() and sha(current) == expected else historical.get(expected)
        if candidate is None:
            candidate = next((p for p in output_sources.get(current.name, []) if sha(p) == expected), None)
        if candidate is None:
            revisions = subprocess.check_output(['git', 'log', '--all', '--format=%H', '--', relative], cwd=ROOT, text=True).splitlines()
            for revision in revisions:
                result = subprocess.run(['git', 'show', revision+':'+relative], cwd=ROOT, capture_output=True)
                if result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == expected:
                    folder = DEST/'historical_sources'
                    folder.mkdir(exist_ok=True)
                    candidate = folder/(current.name+'.'+expected[:16]+'.raw')
                    if not candidate.exists():
                        candidate.write_bytes(result.stdout)
                    break
        assert candidate is not None, '缺少精确构造源码：'+relative
        versions[key] = dict(rec(candidate, expected), restorePath=relative)

    for version in (8, 9, 10, 11, 12, 17, 20, 41):
        protect(OUT/('v83_validation/krylov%02d/wall_velocity_krylov_evidence.mat'%version))
    for relative in (
        'v80_validation/domain02/shell_domain_evidence.mat',
        'v93_validation/continuous01/step01/active_direction.mat',
        'v93_validation/continuous01/step01/direction_context.mat',
        'v100_validation/full04/axisymmetric_patch_evidence.mat',
        'v83_validation/schurcheckpointreal01/checkpoint/density_energy_failure.mat',
        'v83_validation/schurcheckpointreal01/checkpoint/density_energy_failure.mat.json',
        'v83_validation/blockilureal01/CACHE_PRUNED.json',
        'v83_validation/blockilureal01_process.json',
    ):
        protect(OUT/relative)
    protect(ROOT/'docs/validation/storage_20261029/factors/plan.json')
    protect(ROOT/'docs/validation/storage_20261029/factors/postcheck.json')

    for name, removed, job, entry in TARGETS:
        original = safe(OUT/'v83_validation'/name)
        assert original.stat().st_nlink == 1, '本轮只处理无硬链接别名的非冻结文件。'
        assert '/delivery/' not in str(original).replace('\\', '/') and '/frozen/' not in str(original).replace('\\', '/')
        process_file = OUT/'v83_validation'/(job+'_process.json')
        process = read(process_file)
        assert process['status'] == 'terminal' and process['exitCode'] == 0
        assert process['sourceStableAfter'] and all(process['sourceStableAfter'].values())
        protect(process_file)
        if job not in processes:
            bindings = process.get('sourcesBefore', process.get('sourceSHA256'))
            assert any(Path(p).stem == entry for p in bindings), entry
            for relative, digest in bindings.items():
                if relative.endswith('.m') and (not relative.startswith('tests/') or Path(relative).stem in DRIVERS):
                    source(relative, digest)
            processes[job] = rec(process_file)
        retained = RETAINED/'v83_validation'/name
        original_record = rec(original)
        records = copy_records(original, retained, removed)
        metadata = loadmat(retained, variable_names=['report', 'sourceFile', 'parentFile', 'spatialFile', 'primitiveFile', 'sourceBinding', 'sourceIdentity', 'assemblyFile'], simplify_cells=True)
        report = metadata.get('report', {})
        if 'engineProcessID' in report:
            assert report['engineProcessID'] in process['engineProcessIDs']
        if 'sourceSHA256' in report:
            protect(report['sourceFile'], report['sourceSHA256'])
        for key in ('sourceFile', 'parentFile', 'spatialFile', 'primitiveFile', 'assemblyFile'):
            value = metadata.get(key)
            if not isinstance(value, str) or not value:
                continue
            candidate = safe(Path(value))
            if candidate.relative_to(ROOT).as_posix() not in all_targets:
                protect(candidate)
        binding = metadata.get('sourceBinding', {})
        if 'sourceFile' in binding:
            protect(binding['sourceFile'], binding['sourceFileSHA256'])
        for stem, digest in metadata.get('sourceIdentity', {}).items():
            matching = [(p, h) for p, h in process['sourcesBefore'].items() if Path(p).stem == stem and h == digest]
            assert len(matching) == 1, '原构造身份与进程快照不同。'
            source(*matching[0])
        groups.append(dict(original=original_record, retained=rec(retained), removedVariables=sorted(removed),
                           variableRecords=records, rebuildEntry='tests/'+entry+'.m', process=rec(process_file),
                           markerPath=(original.parent/'EMBEDDED_CACHE_PRUNED_20261030.json').relative_to(ROOT).as_posix()))
        print('VERIFIED '+name+' removedRecordBytes='+str(sum(v['bytes'] for v in records if not v['retained'])), flush=True)
    result = dict(status='planned-originals-not-deleted', kind='retired-embedded-MAT-cache-records', workspace=str(ROOT),
                  createdUTC=datetime.now(timezone.utc).isoformat(), gitHead=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  groups=groups, protectedFiles=list(protected.values()), sourceFiles=list(versions.values()),
                  plannedOriginalBytes=sum(g['original']['bytes'] for g in groups),
                  retainedResultBytes=sum(g['retained']['bytes'] for g in groups),
                  plannedRemovedCacheBytes=sum(v['bytes'] for g in groups for v in g['variableRecords'] if not v['retained']),
                  note='所有非缓存变量保留原压缩字节。独有状态、方向、残差、差分证据、失败记录和当前系统保留；历史缓存需按精确源码在新目录重建，旧容器SHA不承诺复现。')
    assert result['plannedOriginalBytes']-result['retainedResultBytes']==result['plannedRemovedCacheBytes']
    save(DEST/'plan.json', result)
    print('PLAN_READY originals='+str(len(groups))+' removedCacheBytes='+str(result['plannedRemovedCacheBytes']), flush=True)


def audit():
    """模块四：删除后独立检查每条保留记录及全部其余输出的文件身份。"""
    no_matlab()
    p = read(DEST/'plan.json')
    e = read(DEST/'execution.json')
    assert e['status'] == 'completed-embedded-cache-cleanup' and e['planSHA256'] == sha(DEST/'plan.json')
    with gzip.open(str(DEST/'inventory_before.json.gz'), 'rb') as stream:
        before = json.loads(stream.read().decode('utf-8'))
    deleted = {g['original']['path'] for g in p['groups']}
    after = snapshot()
    assert set(before['files']) - set(after['files']) == deleted
    for name, identity in before['files'].items():
        if name not in deleted:
            assert after['files'][name] == identity, '非计划原文件变化：'+name
    count = 0
    for group in p['groups']:
        assert not (ROOT/group['original']['path']).exists()
        rec(ROOT/group['retained']['path'], group['retained']['sha256'])
        with (ROOT/group['retained']['path']).open('rb') as stream:
            for item in group['variableRecords']:
                if not item['retained']:
                    continue
                stream.seek(item['retainedOffset'])
                remaining, digest = item['bytes'], hashlib.sha256()
                while remaining:
                    block = stream.read(min(remaining, 4*1024**2))
                    assert block
                    digest.update(block)
                    remaining -= len(block)
                assert digest.hexdigest() == item['rawRecordSHA256']
                count += 1
        expected = [(r['variable'], tuple(r['shape']), r['matlabType']) for r in group['variableRecords'] if r['retained']]
        assert whosmat(ROOT/group['retained']['path']) == expected
    for item in p['protectedFiles']+p['sourceFiles']:
        rec(ROOT/item['path'], item['sha256'])
    save(DEST/'postcheck.json', dict(status='passed', retainedCompressedRecordsVerified=count,
          otherOriginalFilesUnchanged=len(before['files'])-len(deleted), protectedFiles=len(p['protectedFiles']),
          exactSourceVersions=len(p['sourceFiles']), originalMatContainerReproducibilityClaimed=False,
          summaryBefore={k:v for k,v in before.items() if k!='files'},
          summaryAfter={k:v for k,v in after.items() if k!='files'},
          deleteWindowFreeIncreaseBytes=e['actualFreeIncreaseBytes'],
          measuredNetFreeIncreaseBytes=after['freeBytes']-before['freeBytes'],
          netCompressedDataReductionBytes=before['compressedDataBytes']-after['compressedDataBytes']))
    print('AUDIT_PASSED retainedRecords='+str(count)+' actualFreeIncreaseBytes='+str(e['actualFreeIncreaseBytes']), flush=True)


if __name__ == '__main__':
    assert len(sys.argv)==2 and sys.argv[1] in ('plan','audit')
    (plan if sys.argv[1]=='plan' else audit)()
