"""独立核对超旧上限的全部行、完整局部映射和紧凑投影，不赋予物理资格。"""
from pathlib import Path
import os, sys, json
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from screenRocketFullFLeastSquares import ROOT, sha
from auditRocketTemperatureProjection import row_reference


def full_map(factor):
    n = int(factor['cellCount']); P = np.zeros((5*n, 5*n))
    buckets = factor['buckets']
    if isinstance(buckets, dict): buckets = [buckets]
    for k, b in enumerate(buckets, 1):
        count = np.size(b['patchIDs'])
        if not count: continue
        cells = np.asarray(b['members'], int).reshape(k, count, order='F') - 1
        maps = np.asarray(b['map']).reshape(5*k, 5*k, count, order='F')
        for page in range(count):
            ids = (cells[:, page, None] + n*np.arange(5)).ravel()
            P[np.ix_(ids, ids)] = maps[:, :, page]
    return P


def main(folder):
    folder = Path(folder).resolve(); source = folder/'temperature_capacity_unit.mat'
    before = sha(source); s = loadmat(source, simplify_cells=True); checks = []
    for c in s['cases']:
        p = c['projector']; n = int(c['factor']['cellCount']); N = 5*n
        P = full_map(c['factor']); H = c['rows']['H'].toarray()
        ids = np.atleast_1d(p['basisCoordinateIndices']).astype(int)-1
        U = np.zeros((N, int(p['independentConstraintCount'])))
        U[ids] = np.asarray(p['basis']).reshape(len(ids), -1, order='F')
        reps = np.atleast_1d(p['exactRepresentativeRowIndices']).astype(int)-1
        mapping = np.atleast_1d(p['originalToUniqueRow']).astype(int)-1
        np.testing.assert_array_equal(H, H[reps[mapping]])
        B = P.T@H[reps].T
        outside = np.setdiff1d(np.arange(N), ids)
        assert np.all(B[outside] == 0), '紧凑支撑外不得存在被忽略的P转置H转置项'
        span = float(np.linalg.norm(B-U@(U.T@B))/max(1., np.linalg.norm(B)))
        assert span < 1e-12
        Pi = np.eye(N)-U@U.T
        assert np.linalg.norm(Pi-c['Pi']) < 1e-12
        assert np.linalg.norm(H@P@Pi) < 1e-10*max(1., np.linalg.norm(H@P))
        base = c['source']; rhs = base['base']['rhs'].ravel(order='F')
        A = np.row_stack([base['base']['C'], c['factor']['penaltyInOriginalResidualUnits']*base['W']])
        Z = null_space(H); dense = Z@np.linalg.lstsq(A@Z, np.r_[rhs, np.zeros(N)], rcond=None)[0]
        value = c['value'].ravel(order='F')
        error = float(np.linalg.norm(value-dense)/max(1., np.linalg.norm(dense)))
        assert error < 1e-8 and c['nonduplicateDependenceRejected']
        np.testing.assert_allclose(P@Pi@c['coordinates'].ravel(order='F'), value, rtol=1e-10, atol=1e-10)
        assert p['additionalConstructionArrayBoundBytes'] <= 448*1024**2
        checks.append(dict(originalRows=len(H), distinctConditions=len(reps), supportCoordinates=len(ids),
                           completeRightMapSpanDifference=span, denseConstrainedDifference=error,
                           allOmittedCoordinatesExactlyZero=True, nonduplicateDependenceStillRejected=True))
    reference = row_reference(s['realRows']) if 'realRows' in s else None
    if reference is not None: assert len(np.atleast_1d(s['realRows']['faceIDs'])) == 17
    assert sha(source) == before
    out = dict(status='passed-independent-compact-support-full-P-dense-objective-and-public-current-rows',
               cases=checks, publicCurrentRows=reference, inputSHA256=before, scriptSHA256=sha(__file__),
               physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_temperature_capacity.json'; assert not target.exists()
    target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8'); print(json.dumps(out))


if __name__ == '__main__': main(sys.argv[1])
