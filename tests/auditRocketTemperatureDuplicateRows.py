"""独立核对原行精确重构、QR右映射、全约束投影及稠密受限目标。"""
from pathlib import Path
import json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from screenRocketFullFLeastSquares import ROOT, sha

folder = Path(sys.argv[1]).resolve()
source = folder / 'temperature_duplicate_unit.mat'
before = sha(source)
s = loadmat(source, simplify_cells=True)
cases = s['cases'] if isinstance(s['cases'], list) else list(np.atleast_1d(s['cases']))
checks = []
for c in cases:
    projector = c['projector']; H = c['rows']['H'].toarray(); n = int(c['factor']['cellCount']); N = 5 * n
    reps = np.atleast_1d(projector['exactRepresentativeRowIndices']).astype(int) - 1
    mapping = np.atleast_1d(projector['originalToUniqueRow']).astype(int) - 1
    np.testing.assert_array_equal(H, H[reps[mapping]])
    assert len(H) == 7 and len(reps) == 3 and projector['constraintCount'] == 7
    assert projector['allOriginalRowsExactlyRepresented'] and not projector['rankTruncationUsed']
    assert c['nonduplicateDependenceRejected']
    P = np.zeros((N, N))
    buckets = c['factor']['buckets']
    if isinstance(buckets, dict): buckets = [buckets]
    for k, bucket in enumerate(buckets, 1):
        count = np.size(bucket['patchIDs'])
        if not count: continue
        members = np.asarray(bucket['members'], int).reshape(k, count, order='F') - 1
        maps = np.asarray(bucket['map']).reshape(5*k, 5*k, count, order='F')
        for page in range(count):
            ids = (members[:, page, None] + n*np.arange(5)).ravel()
            P[np.ix_(ids, ids)] = maps[:, :, page]
    np.testing.assert_allclose(P, c['P'], rtol=1e-14, atol=1e-14)
    B = P.T @ H[reps].T
    U, _ = np.linalg.qr(B / np.linalg.norm(B, axis=0), mode='reduced')
    Pi = np.eye(N) - U @ U.T
    error = float(np.linalg.norm(Pi-c['Pi']))
    assert error < 1e-12 and np.linalg.norm(H@P@Pi) < 1e-10*np.linalg.norm(H@P)
    base = c['source']; rhs = base['base']['rhs'].ravel(order='F')
    A = np.row_stack([base['base']['C'], c['factor']['penaltyInOriginalResidualUnits']*base['W']])
    Z = null_space(H); dense = Z @ np.linalg.lstsq(A@Z, np.r_[rhs, np.zeros(N)], rcond=None)[0]
    value = c['value'].ravel(order='F')
    dense_error = float(np.linalg.norm(value-dense)/max(1., np.linalg.norm(dense)))
    assert dense_error < 1e-8
    np.testing.assert_allclose(P@Pi@c['coordinates'].ravel(order='F'), value, rtol=1e-10, atol=1e-10)
    assert c['detail']['adjointIdentityRelativeDifference'] < 1e-12
    checks.append(dict(allSevenOriginalRowsExactlyRepresented=True, independentConditions=3,
                       symmetricProjectorDifference=error, denseConstrainedDifference=dense_error,
                       nonduplicateLinearDependenceStillRejected=True))
assert sha(source) == before
out = dict(status='passed-independent-exact-duplicate-constraints-full-P-and-dense-objective',
           checks=checks, inputSHA256=before, scriptSHA256=sha(__file__),
           physicalSolution=False, newProductionSolvedAreaM2=0)
target = folder / 'independent_temperature_duplicates.json'; assert not target.exists()
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
