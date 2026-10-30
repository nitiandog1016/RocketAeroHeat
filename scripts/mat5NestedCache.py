"""流式处理MAT5标量结构体，只删白名单构造字段，其余矩阵原字节保留。"""

import hashlib
import struct
import zlib


class InflatedReader:
    """模块一：限制解压缓冲，避免为旧GB级MAT加载整套数值数组。"""

    def __init__(self, file, compressed_bytes):
        self.file = file
        self.remaining = compressed_bytes
        self.decoder = zlib.decompressobj()
        self.pending = b''
        self.buffer = b''
        self.position = 0

    def fill(self):
        while not self.buffer:
            if not self.pending:
                if not self.remaining:
                    return False
                self.pending = self.file.read(min(self.remaining, 256 * 1024))
                assert self.pending, '压缩记录提前结束。'
                self.remaining -= len(self.pending)
            self.buffer = self.decoder.decompress(self.pending, 1024 * 1024)
            self.pending = self.decoder.unconsumed_tail
        return True

    def read(self, count):
        assert 0 <= count <= 1024 * 1024, '单次读取超过流式预算。'
        blocks = []
        while count:
            assert self.fill(), '解压记录提前结束。'
            size = min(count, len(self.buffer))
            blocks.append(self.buffer[:size])
            self.buffer = self.buffer[size:]
            self.position += size
            count -= size
        return b''.join(blocks)

    def consume(self, count, digest=None, sink=None):
        while count:
            block = self.read(min(count, 256 * 1024))
            if digest is not None:
                digest.update(block)
            if sink is not None:
                sink(block)
            count -= len(block)

    def finish(self):
        assert not self.fill(), '矩阵结束后存在未列明原始数据。'
        assert not self.remaining and not self.pending and self.decoder.eof
        assert not self.decoder.unused_data, '压缩记录存在额外尾部。'


def tag(reader):
    raw = reader.read(8)
    value, size = struct.unpack('<II', raw)
    if value >> 16:
        assert 0 < value >> 16 <= 4
        return value & 65535, value >> 16, raw[4:4 + (value >> 16)], raw
    return value, size, None, raw


def element(reader):
    kind, size, small, raw = tag(reader)
    assert size <= 1024 * 1024, '结构元数据超出预算。'
    if small is None:
        data = reader.read(size)
        raw += data + reader.read((-size) % 8)
    else:
        data = small
    return kind, data, raw


def descriptor(reader):
    """模块二：读取MAT类型、形状、名字与绝对记录范围。"""
    start = reader.position
    kind, size, small, raw = tag(reader)
    assert kind == 14 and small is None and size % 8 == 0, '仅接受完整miMATRIX记录。'
    flag_kind, flags, flag_raw = element(reader)
    dim_kind, dims, dim_raw = element(reader)
    name_kind, name, name_raw = element(reader)
    assert flag_kind == 6 and len(flags) == 8 and dim_kind == 5 and len(dims) % 4 == 0
    assert name_kind == 1
    return dict(start=start, end=start + 8 + size, bytes=size + 8,
                matlabClass=struct.unpack('<II', flags)[0] & 255,
                shape=list(struct.unpack('<' + 'i' * (len(dims) // 4), dims)),
                name=name.decode('utf-8'), prefix=raw + flag_raw + dim_raw + name_raw)


def struct_fields(reader, item):
    assert item['matlabClass'] == 2 and item['shape'] == [1, 1], '只处理标量结构体。'
    kind, stride_raw, stride_element = element(reader)
    assert kind == 5 and len(stride_raw) == 4
    stride = struct.unpack('<i', stride_raw)[0]
    field_start = reader.position
    kind, names, field_raw = element(reader)
    assert kind == 1 and 0 < stride <= 256 and len(names) % stride == 0
    fields = [names[i:i + stride].split(b'\0')[0].decode('utf-8')
              for i in range(0, len(names), stride)]
    assert len(fields) == len(set(fields)) and all(fields)
    return fields, dict(start=field_start, bytes=len(field_raw), stride=stride,
                        prefixRaw=stride_element + field_raw, valuesStart=reader.position)


def scan_result(file, offset, compressed_bytes, cache_fields, require_removed=True):
    """模块三：给每个科学字段计算原矩阵SHA；白名单缓存另行列明。"""
    file.seek(offset)
    reader = InflatedReader(file, compressed_bytes)
    root = descriptor(reader)
    assert root['name'] == 'result'
    fields, root_meta = struct_fields(reader, root)
    assert 'firstOperator' in fields
    kept, removed, operator, metadata = {}, {}, None, None

    def leaf(item, key, target):
        digest = hashlib.sha256(item['prefix'])
        reader.consume(item['end'] - reader.position, digest=digest)
        target[key] = {k: item[k] for k in ('start', 'end', 'bytes', 'matlabClass', 'shape')}
        target[key]['sha256'] = digest.hexdigest()

    for name in fields:
        item = descriptor(reader)
        if name != 'firstOperator':
            leaf(item, 'result.' + name, kept)
            continue
        operator = item
        children, metadata = struct_fields(reader, item)
        assert not require_removed or set(cache_fields).intersection(children), '没有找到白名单矩阵。'
        for child in children:
            value = descriptor(reader)
            key = 'result.firstOperator.' + child
            if child in cache_fields:
                assert value['matlabClass'] == (2 if child == 'wlsNormalFactor' else 6)
                if child != 'wlsNormalFactor':
                    assert len(value['shape']) == 3 and value['shape'][:2] == [5, 5]
                if child != 'wlsNormalFactor':
                    leaf(value, key, removed)
                else:
                    factors, factor_meta = struct_fields(reader, value)
                    assert set(factors) == {'inverse', 'upper1', 'upper2', 'effectiveLower1'}
                    digest = hashlib.sha256(value['prefix'] + factor_meta['prefixRaw'])
                    for factor in factors:
                        block = descriptor(reader)
                        assert block['matlabClass'] == 6 and len(block['shape']) == 3 and block['shape'][:2] == [5, 5]
                        digest.update(block['prefix'])
                        reader.consume(block['end'] - reader.position, digest=digest)
                    assert reader.position == value['end']
                    removed[key] = {k: value[k] for k in ('start', 'end', 'bytes', 'matlabClass', 'shape')}
                    removed[key].update(sha256=digest.hexdigest(), factorFields=factors)
            else:
                leaf(value, key, kept)
        assert reader.position == operator['end']
        metadata['fields'] = children
    assert reader.position == root['end']
    reader.finish()
    return dict(root={k: root[k] for k in ('start', 'end', 'bytes')},
                operator={k: operator[k] for k in ('start', 'end', 'bytes')},
                fields=fields, operatorFields=metadata['fields'],
                fieldTable={k: v for k, v in metadata.items() if k != 'prefixRaw'},
                kept=kept, removed=removed)


def patches(scan):
    """模块四：仅重写两级结构长度和字段表；保留字段的原数据元素不变。"""
    metadata = scan['fieldTable']
    removed_names = {key.rsplit('.', 1)[1] for key in scan['removed']}
    names = [name for name in scan['operatorFields'] if name not in removed_names]
    stride = metadata['stride']
    data = b''.join(name.encode('utf-8').ljust(stride, b'\0') for name in names)
    new_table = struct.pack('<II', 1, len(data)) + data + b'\0' * ((-len(data)) % 8)
    decrease = sum(item['bytes'] for item in scan['removed'].values()) + metadata['bytes'] - len(new_table)
    edits = [(scan['root']['start'] + 4, 4, struct.pack('<I', scan['root']['bytes'] - 8 - decrease)),
             (scan['operator']['start'] + 4, 4, struct.pack('<I', scan['operator']['bytes'] - 8 - decrease)),
             (metadata['start'], metadata['bytes'], new_table)]
    edits += [(value['start'], value['bytes'], b'') for value in scan['removed'].values()]
    edits.sort()
    assert all(a[0] + a[1] <= b[0] for a, b in zip(edits, edits[1:]))
    return edits


def prune_result(source, destination, offset, compressed_bytes, scan):
    """模块五：流式重压缩修改过的结构；不重算或重新序列化任何科学数组。"""
    source.seek(offset)
    reader = InflatedReader(source, compressed_bytes)
    begin = destination.tell()
    destination.write(b'\0' * 8)
    encoder = zlib.compressobj(6)

    def write(raw):
        data = encoder.compress(raw)
        if data:
            destination.write(data)

    for position, count, replacement in patches(scan):
        reader.consume(position - reader.position, sink=write)
        reader.consume(count)
        write(replacement)
    reader.consume(scan['root']['end'] - reader.position, sink=write)
    reader.finish()
    destination.write(encoder.flush())
    end = destination.tell()
    destination.seek(begin)
    destination.write(struct.pack('<II', 15, end - begin - 8))
    destination.seek(end)
    return dict(offset=begin, bytes=end - begin)


def same_science(original, retained):
    """模块六：独立扫描输出后，比较每个保留字段的原字节、类型及形状。"""
    assert set(original['kept']) == set(retained['kept'])
    assert not retained['removed']
    for key, value in original['kept'].items():
        actual = retained['kept'][key]
        assert all(value[k] == actual[k] for k in ('bytes', 'matlabClass', 'shape', 'sha256')), key
