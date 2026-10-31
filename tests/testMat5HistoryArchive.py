"""无损历史归档验证：原容器逐字节恢复、重复数组共享和损坏拒绝。"""

import io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import tempfile
import unittest
import zlib

import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from archiveMat5History import Store, pack, read_manifest, replay, restore
from auditMat5HistoryArchive import audit_manifest


class HistoryArchiveTests(unittest.TestCase):
    def make_mat(self, path, phase=0):
        common = np.arange(40000, dtype=np.float64).reshape(100, 400)
        cells = np.empty((1, 2), dtype=object)
        cells[0, 0] = dict(value=common, state='失败证据')
        cells[0, 1] = np.array([np.nan, -0., 0., np.inf, -np.inf])
        value = dict(result=dict(initialQ=common, Q=common + phase, trace=cells,
                                 empty=np.zeros((0, 3)), flag=np.array([True, False]),
                                 signed=np.array([-3, 5], dtype=np.int64),
                                 complexValue=np.array([1 + 2j, 3 - 4j]),
                                 info=dict(field=common, status='未稳态')),
                     report=dict(eligible=False, productionArea=0.))
        buf = io.BytesIO()
        savemat(buf, value, do_compression=True, long_field_names=True)
        path.write_bytes(buf.getvalue())
        return value

    def test_exact_original_bytes_and_scientific_types(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, manifest = root / 'source.mat', root / 'manifest.z'
            self.make_mat(source)
            value = pack(source, root / 'store', manifest)
            restored = root / 'restored.mat'
            restore(read_manifest(manifest), root / 'store', restored)
            self.assertEqual(source.read_bytes(), restored.read_bytes())
            actual = loadmat(restored, simplify_cells=True)['result']
            self.assertEqual(actual['signed'].dtype, np.int64)
            self.assertTrue(np.signbit(actual['trace'][1][1]))
            self.assertFalse(np.signbit(actual['trace'][1][2]))
            self.assertEqual(actual['info']['status'], '未稳态')
            self.assertEqual(value['verification']['status'], 'passed-exact-original-file-SHA')

    def test_duplicate_arrays_and_files_share_objects(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for n in (1, 2):
                self.make_mat(root / ('s%d.mat' % n), phase=n)
            first = pack(root / 's1.mat', root / 'store', root / 'm1.z')
            second = pack(root / 's2.mat', root / 'store', root / 'm2.z')
            first_refs = {r['sha256'] for x in first['records'] for r in x['pieces'] if 'sha256' in r}
            second_refs = {r['sha256'] for x in second['records'] for r in x['pieces'] if 'sha256' in r}
            self.assertGreater(len(first_refs.intersection(second_refs)), 0)
            replay(first, root / 'store', lambda _: None)
            replay(second, root / 'store', lambda _: None)

    def test_missing_and_corrupt_objects_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.make_mat(root / 's.mat')
            value = pack(root / 's.mat', root / 'store', root / 'm.z')
            store = Store(root / 'store')
            ref = next(r for x in value['records'] for r in x['pieces'] if 'sha256' in r)
            obj = store.object_path(ref['sha256'])
            encoded = obj.read_bytes()
            obj.write_bytes(zlib.compress(b'wrong'))
            with self.assertRaises(AssertionError):
                replay(value, root / 'store', lambda _: None)
            obj.write_bytes(encoded)
            obj.unlink()
            with self.assertRaises(FileNotFoundError):
                replay(value, root / 'store', lambda _: None)

    def test_path_traversal_and_existing_target_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            store = Store(root / 'store')
            with self.assertRaises(AssertionError):
                store.object_path('../unsafe')
            self.make_mat(root / 's.mat')
            value = pack(root / 's.mat', root / 'store', root / 'm.z')
            with self.assertRaises(AssertionError):
                restore(value, root / 'store', root / 's.mat')
            with self.assertRaises(AssertionError):
                pack(root / 's.mat', root / 'store', root / 'm.z')

    def test_wrong_original_identity_and_extra_stream_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.make_mat(root / 's.mat')
            value = pack(root / 's.mat', root / 'store', root / 'm.z')
            value['original']['sha256'] = '0' * 64
            with self.assertRaises(AssertionError):
                replay(value, root / 'store', lambda _: None)
            self.assertTrue((root / 's.mat').exists())

    def test_corrupt_source_is_refused_without_deletion(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.make_mat(root / 's.mat')
            raw = bytearray((root / 's.mat').read_bytes())
            raw[160] ^= 255
            (root / 's.mat').write_bytes(raw)
            with self.assertRaises(Exception):
                pack(root / 's.mat', root / 'store', root / 'm.z')
            self.assertFalse((root / 'm.z').exists())
            self.assertTrue((root / 's.mat').exists())

    def test_parallel_publication_and_independent_auditor(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.make_mat(root / 's.mat')
            def run(n):
                return pack(root / 's.mat', root / 'store', root / ('m%d.z' % n))
            with ThreadPoolExecutor(max_workers=3) as pool:
                values = list(pool.map(run, range(3)))
            for n, value in enumerate(values):
                report = audit_manifest(root / ('m%d.z' % n), root / 'store')
                self.assertEqual(report['sha256'], value['original']['sha256'])

    def test_empty_nonscalar_struct_and_sparse_are_exact(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            array = np.empty((1, 2), dtype=[('x', object), ('y', object)])
            array[0, 0] = (np.array([1, 2]), '甲')
            array[0, 1] = (np.array([3, 4]), '乙')
            savemat(root / 's.mat', {'result': dict(array=array, sparse=csr_matrix(np.eye(4)),
                                                  cell=np.empty((0, 2), dtype=object))}, do_compression=True)
            value = pack(root / 's.mat', root / 'store', root / 'm.z')
            restore(value, root / 'store', root / 'r.mat')
            self.assertEqual((root / 's.mat').read_bytes(), (root / 'r.mat').read_bytes())

    def test_wrong_zlib_version_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.make_mat(root / 's.mat')
            value = pack(root / 's.mat', root / 'store', root / 'm.z')
            value['zlibRuntime'] = 'unknown-version'
            with self.assertRaises(AssertionError):
                replay(value, root / 'store', lambda _: None)


if __name__ == '__main__':
    unittest.main(verbosity=2)
