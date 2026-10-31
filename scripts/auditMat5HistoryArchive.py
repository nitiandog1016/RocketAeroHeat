"""独立读取归档对象并重构原MAT的SHA；不调用归档器的重放函数，不删除文件。"""

import base64
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zlib


def audit_manifest(path, store):
    """独立按记录顺序重压缩；完整原文件SHA相同才接受迁存。"""
    value = json.loads(zlib.decompress(Path(path).read_bytes()).decode('utf-8'))
    assert value['format'] == 'rocket-mat5-lossless-dedup-v1'
    assert value['zlibRuntime'] == zlib.ZLIB_RUNTIME_VERSION
    whole, count = hashlib.sha256(), 0

    def accept(raw):
        nonlocal count
        whole.update(raw)
        count += len(raw)

    def raw_piece(piece):
        if 'inline' in piece:
            raw = base64.b64decode(piece['inline'], validate=True)
        else:
            digest = piece['sha256']
            assert len(digest) == 64 and all(c in '0123456789abcdef' for c in digest)
            encoded = (Path(store) / 'objects' / digest[:2] / (digest + '.z')).read_bytes()
            raw = zlib.decompress(encoded)
            assert hashlib.sha256(raw).hexdigest() == digest
        assert len(raw) == piece['bytes']
        return raw

    accept(base64.b64decode(value['header'], validate=True))
    for record in value['records']:
        if record['mode'] == 'raw':
            for piece in record['pieces']:
                accept(raw_piece(piece))
        else:
            assert record['mode'] == 'zlib6'
            accept(base64.b64decode(record['tag'], validate=True))
            encoder = zlib.compressobj(6)
            inflated = hashlib.sha256()
            for piece in record['pieces']:
                raw = raw_piece(piece)
                inflated.update(raw)
                accept(encoder.compress(raw))
            accept(encoder.flush())
            assert inflated.hexdigest() == record['inflatedSHA256']
    assert count == value['original']['bytes']
    assert whole.hexdigest() == value['original']['sha256']
    return {'original': value['original']['path'], 'sha256': whole.hexdigest(), 'bytes': count}


def main():
    root = Path(__file__).resolve().parents[1]
    dest = root / 'docs/validation/storage_20261031/output_review'
    plan = json.loads((dest / 'plan.json').read_text(encoding='utf-8'))
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(lambda g: audit_manifest(root / g['manifest']['path'], root / plan['store']), plan['groups']))
    report = {'status': 'passed-independent-exact-original-MAT-SHA', 'verifiedOriginals': len(results),
              'verifiedOriginalBytes': sum(x['bytes'] for x in results), 'results': results,
              'freeBytes': shutil.disk_usage(str(root)).free, 'archiveReplayFunctionUsed': False}
    (dest / ('postdelete_archive_audit.json' if '--postdelete' in sys.argv else 'predelete_archive_audit.json')).write_bytes(
        (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
