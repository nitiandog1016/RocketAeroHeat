"""核对图预条件缓存迁存：线分解可删，科学值与分组身份必须原字节保留。"""

import io
from pathlib import Path
import struct
import sys
import unittest

import numpy as np
from scipy.io import loadmat, savemat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compactRetiredGraphOperators import CACHE_FIELDS
from mat5NestedCache import scan_result, prune_result, same_science


class RetiredGraphCacheTests(unittest.TestCase):
    def fixture(self, invalid_shape=False):
        blocks = np.arange(75, dtype=np.float64).reshape(5, 5, 3)
        operator = dict(diagonal=blocks, lineInverse=blocks+2, lineUpperFactor=blocks-1,
                        sideFaceIDs=np.array([11, 12, 13], dtype=np.int32),
                        colorGroups=np.array([1, 2, 1], dtype=np.int32),
                        maximumFactorDefect=3e-14,
                        report=dict(physicalSolution=False, status='未稳态'))
        if invalid_shape:
            operator['lineInverse'] = np.ones((4, 4, 3))
        result = dict(Q=np.array([[np.nan, -0., 0., np.inf, -np.inf]]),
                      info=dict(wallPressure=np.array([11., 13.]),
                                wallHeatFlux=np.array([np.nan, 7.])),
                      firstOperator=operator, failures=np.array([2, 9], dtype=np.int64))
        source = io.BytesIO()
        savemat(source, dict(result=result), do_compression=True, long_field_names=True)
        return source.getvalue()

    def test_graph_blocks_removed_scientific_records_exact(self):
        data = self.fixture()
        count = struct.unpack('<II', data[128:136])[1]
        before = scan_result(io.BytesIO(data), 136, count, CACHE_FIELDS)
        self.assertEqual(set(before['removed']), {
            'result.firstOperator.diagonal', 'result.firstOperator.lineInverse',
            'result.firstOperator.lineUpperFactor'})
        output = io.BytesIO(data[:128]); output.seek(128)
        record = prune_result(io.BytesIO(data), output, 136, count, before)
        after = scan_result(io.BytesIO(output.getvalue()), record['offset']+8,
                            record['bytes']-8, CACHE_FIELDS, require_removed=False)
        same_science(before, after)
        actual = loadmat(io.BytesIO(output.getvalue()), simplify_cells=True)['result']
        self.assertTrue(np.signbit(actual['Q'][1]))
        self.assertFalse(np.signbit(actual['Q'][2]))
        self.assertTrue(np.isnan(actual['info']['wallHeatFlux'][0]))
        self.assertEqual(set(actual['firstOperator']), {
            'sideFaceIDs', 'colorGroups', 'maximumFactorDefect', 'report'})
        self.assertEqual(actual['firstOperator']['sideFaceIDs'].dtype, np.int32)
        self.assertEqual(actual['firstOperator']['report']['status'], '未稳态')
        self.assertEqual(actual['failures'].dtype, np.int64)

    def test_non_5_by_5_line_cache_is_refused(self):
        data = self.fixture(invalid_shape=True)
        count = struct.unpack('<II', data[128:136])[1]
        with self.assertRaises(AssertionError):
            scan_result(io.BytesIO(data), 136, count, CACHE_FIELDS)


if __name__ == '__main__':
    unittest.main(verbosity=2)
