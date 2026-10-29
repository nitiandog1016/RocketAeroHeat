# -*- coding: utf-8 -*-
"""对真实修改前后MAT执行完整公共字段比较，排除且只排除两项新增分辨率元数据。"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.io import loadmat
from auditThermalPressureSurfaceFlow import sha


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('before'); parser.add_argument('after'); parser.add_argument('--report', required=True)
    args = parser.parse_args(); target = Path(args.report); assert not target.exists()
    a = loadmat(args.before, simplify_cells=True); b = loadmat(args.after, simplify_cells=True)
    excluded = {'entryDerivativeResolution', 'maximumRelativeEntryDerivativeDifference'}; count = 0

    def compare(x, y, path):
        nonlocal count
        if isinstance(x, dict):
            assert isinstance(y, dict) and set(x)-excluded == set(y)-excluded, path
            for key in x:
                if key not in excluded: compare(x[key], y[key], path+'.'+key)
        elif isinstance(x, (list, tuple)) or isinstance(x, np.ndarray) and x.dtype == object:
            aa = list(x.flat) if isinstance(x, np.ndarray) else x; bb = list(y.flat) if isinstance(y, np.ndarray) else y
            assert len(aa) == len(bb), path
            for i, (u, v) in enumerate(zip(aa, bb)): compare(u, v, path+'['+str(i)+']')
        else:
            try: same = np.array_equal(x, y, equal_nan=True)
            except TypeError: same = np.array_equal(x, y)
            assert same, path; count += 1

    compare(a['cases'], b['cases'], 'cases'); compare(a['prefix'], b['prefix'], 'prefix')
    report = dict(status='passed-exact-physical-and-case-payload-comparison-before-after-inlet-resolution-gate',
                  savedObjectsExactlyEqual=count, excludedOnlyNewMetadata=sorted(excluded), noNonzeroNumericToleranceUsed=True,
                  beforeSHA256=sha(args.before), afterSHA256=sha(args.after), scriptSHA256=sha(__file__))
    target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf8'))
    print(json.dumps(report))


if __name__ == '__main__':
    main()
