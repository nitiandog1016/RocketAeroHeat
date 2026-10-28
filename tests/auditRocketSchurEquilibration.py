# -*- coding: utf-8 -*-
"""独立遍历真实Schur的全部系数，核对可逆R/S/C及原单位还原；不证明数值满秩。"""
from pathlib import Path
import os, sys, json, hashlib
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
case = sys.argv[1]
assert case.replace('_', '').isalnum()
target = VAL / (case + '_independent_equilibration.json')
assert not target.exists()

def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            result.update(block)
    return result.hexdigest()

evidence = VAL / case / 'schur_equilibration_evidence.mat'
saved = loadmat(evidence, simplify_cells=True)
assembly = Path(saved['assemblyFile'])
assert sha(assembly) == saved['assemblySHA256']
source = ROOT / 'internalFlow/cutCell/equilibrateRocketSchurScales.m'
assert sha(source) == saved['scaleSourceSHA256']
data = loadmat(assembly, simplify_cells=True, variable_names=['schurMatrix', 'sourceBinding'])
assert data['sourceBinding'] == saved['sourceBinding']
m = data['schurMatrix']
n = int(m['cellCount'])
column = m['column'].astype(np.int64) - 1
pointer = m['rowPointer'].astype(np.int64) - 1
owners = np.repeat(np.arange(n, dtype=np.int32), np.diff(pointer))
values = m['values'].transpose(2, 0, 1)
row_max = np.zeros((n, 2))
column_max = np.zeros((n, 2))
for a in range(2):
    for b in range(2):
        np.maximum.at(row_max[:, a], owners, abs(values[:, a, b]))
assert np.all(np.isfinite(row_max) & (row_max > 0))
R = 1 / row_max
for a in range(2):
    for b in range(2):
        np.maximum.at(column_max[:, b], column, abs(values[:, a, b]) * R[owners, a])
assert np.all(np.isfinite(column_max) & (column_max > 0))
C = 1 / column_max
np.testing.assert_allclose(saved['rowScale'], R, rtol=1e-14, atol=0)
np.testing.assert_allclose(saved['columnScale'], C, rtol=1e-14, atol=0)
scaled_rows = np.zeros((n, 2))
scaled_columns = np.zeros((n, 2))
max_roundtrip = 0.
original_nonzeros = scaled_nonzeros = 0
for first in range(0, len(column), 4096):
    ids = slice(first, first + 4096)
    raw = values[ids]
    balanced = raw * R[owners[ids], :, None] * C[column[ids], None, :]
    assert np.all(np.isfinite(balanced))
    np.testing.assert_array_equal(raw != 0, balanced != 0)
    restored = balanced / R[owners[ids], :, None] / C[column[ids], None, :]
    mask = raw != 0
    max_roundtrip = max(max_roundtrip, float(np.max(abs((restored[mask] - raw[mask]) / raw[mask]))))
    original_nonzeros += int(np.count_nonzero(raw))
    scaled_nonzeros += int(np.count_nonzero(balanced))
    for field in range(2):
        np.maximum.at(scaled_rows[:, field], owners[ids], np.max(abs(balanced[:, field, :]), axis=1))
        np.maximum.at(scaled_columns[:, field], column[ids], np.max(abs(balanced[:, :, field]), axis=1))
assert max_roundtrip < 1e-14
assert np.max(abs(scaled_rows - 1)) < 1e-14 and np.max(abs(scaled_columns - 1)) < 1e-14
assert sha(assembly) == saved['assemblySHA256'] and sha(source) == saved['scaleSourceSHA256']
result = dict(
    status='passed-independent-complete-real-Schur-row-column-equilibration-and-original-unit-roundtrip',
    processID=os.getpid(), cellCount=n, scalarUnknowns=2*n, allInputBlocks=len(column),
    originalScalarNonzeros=original_nonzeros, scaledScalarNonzeros=scaled_nonzeros,
    rowScaleRange=[float(R.min()), float(R.max())], columnScaleRange=[float(C.min()), float(C.max())],
    maximumRelativeCoefficientRoundtripDifference=max_roundtrip,
    equilibratedRowInfinityNormRange=[float(scaled_rows.min()), float(scaled_rows.max())],
    equilibratedColumnInfinityNormRange=[float(scaled_columns.min()), float(scaled_columns.max())],
    coefficientsAddedOrDeleted=False, numericalFullRankOrStableFactorizationProven=False,
    newSelectedJApplications=0, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    physicalSolution=False, inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (assembly, evidence, source)},
    scriptSHA256=sha(__file__))
target.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
print(json.dumps(result))
