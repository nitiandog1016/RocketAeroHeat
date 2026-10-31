"""旧MAT结果无损去重：保留所有数值，并在删除前验证可恢复原文件SHA。"""

import argparse
import base64
import hashlib
import json
from pathlib import Path
import struct
import threading
import zlib

from mat5NestedCache import InflatedReader, descriptor, element

BLOCK = 1024 * 1024
INLINE = 8192
STORE_LOCK = threading.Lock()


def digest_file(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(BLOCK), b''):
            value.update(block)
    return value.hexdigest()


def plain_path(path):
    """模块一：拒绝链接跳转；写入仅使用新文件，不覆盖已有结果。"""
    path = Path(path).absolute()
    for parent in (path, *path.parents):
        if parent.exists():
            assert not (getattr(parent.lstat(), 'st_file_attributes', 0) & 0x400), str(parent)
            assert not parent.is_symlink(), str(parent)
    return path


class Store:
    """模块二：相同原字节只存一次，压缩对象逐项校验SHA和长度。"""

    def __init__(self, folder):
        self.folder = plain_path(folder)
        self.known = set()
        self.added_bytes = 0

    def object_path(self, digest):
        assert len(digest) == 64 and all(c in '0123456789abcdef' for c in digest)
        return plain_path(self.folder / 'objects' / digest[:2] / (digest + '.z'))

    def put(self, raw):
        if len(raw) < INLINE:
            return {'inline': base64.b64encode(raw).decode('ascii'), 'bytes': len(raw)}
        # 并行迁存时串行发布完整对象，另一线程不得读到尚未写完的数据。
        with STORE_LOCK:
            return self.put_large(raw)

    def put_large(self, raw):
        digest = hashlib.sha256(raw).hexdigest()
        target = self.object_path(digest)
        if digest not in self.known:
            if target.exists():
                assert self.get({'sha256': digest, 'bytes': len(raw)}) == raw
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                encoded = zlib.compress(raw, 1)
                with target.open('xb') as stream:
                    stream.write(encoded)
                self.added_bytes += len(encoded)
            self.known.add(digest)
        return {'sha256': digest, 'bytes': len(raw)}

    def get(self, ref):
        assert isinstance(ref['bytes'], int) and 0 <= ref['bytes'] <= BLOCK
        if 'inline' in ref:
            raw = base64.b64decode(ref['inline'], validate=True)
            assert 'sha256' not in ref and len(raw) < INLINE
        else:
            with self.object_path(ref['sha256']).open('rb') as stream:
                encoded = stream.read(2 * BLOCK)
                assert not stream.read(1)
            decoder = zlib.decompressobj()
            raw = decoder.decompress(encoded, BLOCK + 1)
            assert decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail
            assert hashlib.sha256(raw).hexdigest() == ref['sha256'], '归档对象SHA不符。'
        assert len(raw) == ref['bytes'], '归档对象长度不符。'
        return raw


def split_matrix(reader, store, pieces, depth=0):
    """模块三：按MAT结构/元胞和数组边界分块；原类型、形状、NaN及符号位不变。"""
    assert depth < 64, '嵌套层数超出保守读取范围。'
    item = descriptor(reader)
    pieces.append(store.put(item['prefix']))
    count = 1
    for size in item['shape']:
        assert size >= 0
        count *= size
    if item['matlabClass'] == 2:
        kind, raw_stride, stride_element = element(reader)
        assert kind == 5 and len(raw_stride) == 4
        stride = struct.unpack('<i', raw_stride)[0]
        kind, names, names_element = element(reader)
        assert kind == 1 and 0 < stride <= 256 and len(names) % stride == 0
        pieces.append(store.put(stride_element + names_element))
        children = count * (len(names) // stride)
    elif item['matlabClass'] == 1:
        children = count
    else:
        children = None
    if children is not None:
        assert children < 1000000
        for _ in range(children):
            split_matrix(reader, store, pieces, depth + 1)
    else:
        while reader.position < item['end']:
            pieces.append(store.put(reader.read(min(BLOCK, item['end'] - reader.position))))
    assert reader.position == item['end'], 'MAT元素范围不符。'
    return item


def replay(manifest, folder, sink):
    """模块四：重放原记录；原result使用同版zlib/级别恢复，并核对整个原MAT的SHA。"""
    assert manifest['format'] == 'rocket-mat5-lossless-dedup-v1'
    assert manifest['zlibRuntime'] == zlib.ZLIB_RUNTIME_VERSION
    store = Store(folder)
    whole, total = hashlib.sha256(), 0

    def write(raw):
        nonlocal total
        whole.update(raw)
        total += len(raw)
        sink(raw)

    write(base64.b64decode(manifest['header'], validate=True))
    for record in manifest['records']:
        if record['mode'] == 'raw':
            for ref in record['pieces']:
                write(store.get(ref))
        else:
            assert record['mode'] == 'zlib6'
            write(base64.b64decode(record['tag'], validate=True))
            encoder = zlib.compressobj(6)
            raw_sha = hashlib.sha256()
            encoded_sha, encoded_size = hashlib.sha256(), 0
            for ref in record['pieces']:
                raw = store.get(ref)
                raw_sha.update(raw)
                encoded = encoder.compress(raw)
                encoded_sha.update(encoded)
                encoded_size += len(encoded)
                write(encoded)
            encoded = encoder.flush()
            encoded_sha.update(encoded)
            encoded_size += len(encoded)
            write(encoded)
            assert raw_sha.hexdigest() == record['inflatedSHA256']
            assert encoded_size == record['compressedBytes']
            assert encoded_sha.hexdigest() == record['compressedSHA256'], '不能逐字节恢复原压缩记录。'
    assert total == manifest['original']['bytes']
    assert whole.hexdigest() == manifest['original']['sha256'], '不能恢复整个原文件SHA。'
    return {'status': 'passed-exact-original-file-SHA', 'bytes': total, 'sha256': whole.hexdigest()}


def pack(source, folder, manifest_path, original_relative=None):
    """模块五：只处理明确范围的MAT5；其他记录保留原压缩字节，原文件此步不删除。"""
    source, folder, manifest_path = map(plain_path, (source, folder, manifest_path))
    assert not manifest_path.exists()
    identity = source.stat()
    original = {'path': original_relative or str(source), 'bytes': identity.st_size,
                'sha256': digest_file(source), 'fileID': identity.st_ino,
                'device': identity.st_dev, 'mtimeNS': identity.st_mtime_ns, 'links': identity.st_nlink}
    store = Store(folder)
    records = []
    with source.open('rb') as stream:
        header = stream.read(128)
        assert header.startswith(b'MATLAB 5.0 MAT-file') and header[-2:] == b'IM'
        while stream.tell() < original['bytes']:
            start = stream.tell()
            tag = stream.read(8)
            kind, size = struct.unpack('<II', tag)
            assert kind == 15 and size > 0, '仅接受独立压缩的MAT5顶层记录。'
            reader = InflatedReader(stream, size)
            item = descriptor(reader)
            stream.seek(start + 8)
            if item['name'] == 'result':
                assert item['matlabClass'] == 2
                pieces = []
                reader = InflatedReader(stream, size)
                split_matrix(reader, store, pieces)
                reader.finish()
                raw_sha = hashlib.sha256()
                for ref in pieces:
                    raw_sha.update(store.get(ref))
                stream.seek(start + 8)
                compressed_sha = hashlib.sha256()
                remaining = size
                while remaining:
                    block = stream.read(min(remaining, BLOCK))
                    assert block
                    compressed_sha.update(block)
                    remaining -= len(block)
                record = {'name': item['name'], 'mode': 'zlib6', 'pieces': pieces,
                          'tag': base64.b64encode(tag).decode('ascii'), 'compressedBytes': size,
                          'compressedSHA256': compressed_sha.hexdigest(), 'inflatedSHA256': raw_sha.hexdigest()}
            else:
                stream.seek(start)
                pieces, remaining = [], size + 8
                while remaining:
                    block = stream.read(min(remaining, BLOCK))
                    assert block
                    pieces.append(store.put(block))
                    remaining -= len(block)
                record = {'name': item['name'], 'mode': 'raw', 'pieces': pieces}
            records.append(record)
            assert stream.tell() == start + size + 8
    assert any(r['mode'] == 'zlib6' for r in records)
    now = source.stat()
    assert (now.st_dev, now.st_ino, now.st_size, now.st_mtime_ns) == (
        identity.st_dev, identity.st_ino, identity.st_size, identity.st_mtime_ns)
    value = {'format': 'rocket-mat5-lossless-dedup-v1', 'original': original,
             'header': base64.b64encode(header).decode('ascii'), 'records': records,
             'zlibRuntime': zlib.ZLIB_RUNTIME_VERSION, 'newObjectBytes': store.added_bytes}
    # 原文件仍在：先验证能够恢复其全部字节，验证失败不得产生可删除凭证。
    value['verification'] = replay(value, folder, lambda _: None)
    assert digest_file(source) == original['sha256']
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    encoded = zlib.compress((json.dumps(value, ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf-8'), 6)
    with manifest_path.open('xb') as target:
        target.write(encoded)
    return value


def read_manifest(path):
    return json.loads(zlib.decompress(plain_path(path).read_bytes()).decode('utf-8'))


def restore(manifest, folder, target):
    """模块六：恢复到新路径；已有文件绝不覆盖，失败输出不用于续算。"""
    target = plain_path(target)
    assert not target.exists(), '恢复目标已存在，请改用新的路径。'
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('xb') as stream:
        report = replay(manifest, folder, stream.write)
    assert digest_file(target) == manifest['original']['sha256']
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['pack', 'verify', 'restore'])
    parser.add_argument('--store', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--source')
    parser.add_argument('--target')
    args = parser.parse_args()
    if args.action == 'pack':
        report = pack(args.source, args.store, args.manifest)['verification']
    elif args.action == 'verify':
        report = replay(read_manifest(args.manifest), args.store, lambda _: None)
    else:
        report = restore(read_manifest(args.manifest), args.store, args.target)
    print(json.dumps(report, ensure_ascii=False))


if __name__ == '__main__':
    main()
