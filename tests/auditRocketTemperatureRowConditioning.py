"""以精确有理行关系、独立完整P转置和稠密方向核对二进制行预处理。"""
from pathlib import Path
import os
import sys
import json
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from riemannReservoirReference import ROOT, sha
from rocketTemperatureProjectorReference import conditioned_rows_reference, certify_projector_span


def unit(folder):
    file = folder/'temperature_row_conditioning_unit.mat'
    data = loadmat(file, simplify_cells=True); checks = []
    for c in data['cases']:
        working, proof = conditioned_rows_reference(c['H'], c['conditioning'])
        assert (working-c['workingRows']).nnz == 0
        assert c['proportionalConditionsStillRejected'] and np.any(c['oldNormalizedQRDiagonal'] <= 1e-10)
        p = c['projector']
        assert p['constraintCount'] == p['independentConstraintCount'] == 2 and not p['rankTruncationUsed']
        assert np.all(p['normalizedQRDiagonal'] > 1e-10)
        A = np.row_stack([c['C'], c['penalty']*c['W']]); n = len(c['value'])
        Z = null_space(working.toarray())
        rhs = np.r_[c['rhs'].ravel(order='F'), np.zeros(5*n)]
        truth = Z @ np.linalg.lstsq(A @ Z, rhs, rcond=None)[0]
        error = float(np.linalg.norm(c['value'].ravel(order='F')-truth)/max(1., np.linalg.norm(truth)))
        assert error < 1e-8
        np.testing.assert_allclose(c['H'] @ c['value'].ravel(order='F'), c['detail']['constraintDerivativeK'], atol=1e-8, rtol=1e-8)
        proof['denseDirectionRelativeDifference'] = error; checks.append(proof)
    return dict(cases=checks, inputSHA256={file.relative_to(ROOT).as_posix(): sha(file)})


def real(folder):
    file = folder/'temperature_row_conditioning_real.mat'
    saved = loadmat(file, simplify_cells=True)
    context = ROOT/saved['report']['contextFile']; preview = ROOT/saved['report']['previewFile']
    assert sha(context) == saved['report']['contextSHA256'] and sha(preview) == saved['report']['previewSHA256']
    factor = loadmat(context, simplify_cells=True, variable_names=['factor'])['factor']; n = int(factor['cellCount'])
    plans = []; covered = np.zeros(n, bool)
    for b in factor['buckets']:
        maps = np.asarray(b['map']); size = maps.shape[0]//5
        members = np.asarray(b['members']).astype(int).reshape(size, -1, order='F')-1
        assert not np.any(covered[members.ravel()]); covered[members.ravel()] = True
        count = members.shape[1]; maps = maps.reshape(5*size, 5*size, count, order='F')
        ids = (members.T[:, :, None]+n*np.arange(5)).reshape(count, 5*size).T
        plans.append((ids, maps))
    assert np.all(covered)
    def transpose(raw):
        result = np.zeros(5*n)
        for ids, maps in plans:
            values = raw[ids]; pages = np.flatnonzero(np.any(values != 0, axis=0))
            if len(pages):
                result[ids[:, pages]] = np.einsum('abp,ap->bp', maps[:, :, pages], values[:, pages])
        return result
    projector = saved['projector']; H = saved['rows']['H']
    proof = certify_projector_span(projector, H, transpose)
    reps = np.atleast_1d(projector['exactRepresentativeRowIndices']).astype(int)-1
    working, exact = conditioned_rows_reference(H[reps], saved['conditioning'])
    assert (working-saved['workingRows']).nnz == 0
    return dict(projector=proof, exactPhysicalRows=exact,
                inputSHA256={p.relative_to(ROOT).as_posix(): sha(p) for p in (file, context, preview)})


def main(mode, folder):
    folder = Path(folder).resolve(); result = {'unit': unit, 'real': real}[mode](folder)
    result.update(status='passed-independent-exact-rational-all-condition-row-preconditioning-'+mode,
                  scriptSHA256=sha(__file__), originalRankThreshold=1e-10, newActualNSEvaluations=0,
                  newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_row_conditioning.json'; assert not target.exists()
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print(json.dumps(result))


if __name__ == '__main__':
    main(*sys.argv[1:])
