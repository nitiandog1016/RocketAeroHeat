"""核对嵌套MAT缓存清理：特殊数值、字段类型、损坏记录与白名单拒绝。"""

import io
from pathlib import Path
import struct
import sys
import unittest

import numpy as np
from scipy.io import loadmat, savemat, whosmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from mat5NestedCache import prune_result, same_science, scan_result

CACHES = {'diagonal', 'lower', 'upper', 'wlsNormalFactor'}


class NestedCacheTests(unittest.TestCase):
    def fixture(self, extra_factor=False):
        matrix = np.arange(100, dtype=np.float64).reshape(5, 5, 4)
        factors = dict(inverse=matrix, upper1=matrix, upper2=matrix, effectiveLower1=matrix)
        if extra_factor:
            factors['uniqueResult'] = np.array([999.])
        operator = dict(kind='cached-operator', diagonal=matrix, lower=-matrix, upper=matrix + 1,
                        residualScale=np.array([1., 2., 3., 4., 5.]), wlsNormalFactor=factors,
                        report=dict(sweeps=np.int32(2), physicalSolution=False))
        result = dict(initialQ=np.array([[np.nan, np.inf, -np.inf, -0., 0.]]),
                      Q=np.array([[1., 2., 3., 4., 5.]]),
                      info=dict(pressure=np.array([2., 3.]), status='未稳态'),
                      firstOperator=operator, failures=np.array([-1, 0, 3], dtype=np.int64))
        file = io.BytesIO()
        savemat(file, dict(result=result), do_compression=True, long_field_names=True)
        return file.getvalue(), result

    def test_exact_scientific_records_and_native_reader(self):
        data, expected = self.fixture()
        kind, length = struct.unpack('<II', data[128:136])
        self.assertEqual(kind, 15)
        original = scan_result(io.BytesIO(data), 136, length, CACHES)
        destination = io.BytesIO(data[:128])
        destination.seek(128)
        record = prune_result(io.BytesIO(data), destination, 136, length, original)
        retained = scan_result(io.BytesIO(destination.getvalue()), record['offset'] + 8,
                               record['bytes'] - 8, CACHES, require_removed=False)
        same_science(original, retained)
        actual = loadmat(io.BytesIO(destination.getvalue()), simplify_cells=True)['result']
        np.testing.assert_array_equal(actual['Q'], expected['Q'].ravel())
        np.testing.assert_array_equal(actual['initialQ'], expected['initialQ'].ravel())
        self.assertTrue(np.signbit(actual['initialQ'][3]))
        self.assertFalse(np.signbit(actual['initialQ'][4]))
        self.assertEqual(actual['info']['status'], '未稳态')
        self.assertEqual(actual['failures'].dtype, np.int64)
        self.assertEqual(set(actual['firstOperator']), {'kind', 'residualScale', 'report'})
        self.assertEqual(whosmat(io.BytesIO(data)), whosmat(io.BytesIO(destination.getvalue())))

    def test_extra_wls_result_is_refused(self):
        data, _ = self.fixture(extra_factor=True)
        length = struct.unpack('<II', data[128:136])[1]
        with self.assertRaises(AssertionError):
            scan_result(io.BytesIO(data), 136, length, CACHES)

    def test_corrupt_deflate_is_refused(self):
        data, _ = self.fixture()
        changed = bytearray(data)
        changed[-1] ^= 1
        length = struct.unpack('<II', data[128:136])[1]
        with self.assertRaises(Exception):
            scan_result(io.BytesIO(changed), 136, length, CACHES)

    def test_absent_cache_cannot_authorize_deletion(self):
        data, _ = self.fixture()
        length = struct.unpack('<II', data[128:136])[1]
        with self.assertRaises(AssertionError):
            scan_result(io.BytesIO(data), 136, length, {'absent'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
