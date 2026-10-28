# -*- coding: utf-8 -*-
"""以独立完整 F 正向/转置和固定局部逆筛查有界 LSQR，不改变 MATLAB 生产求解。"""
from pathlib import Path
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import bsr_matrix
from scipy.sparse.linalg import LinearOperator, lsqr

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(part)
    return h.hexdigest()


def matrix(graph):
    n = int(graph['cellCount'])
    np.testing.assert_array_equal(np.atleast_1d(graph['order']), np.arange(1, n+1))
    ptr = np.atleast_1d(graph['rowPointer']).astype(np.int64)-1
    col = np.atleast_1d(graph['column']).astype(np.int64)-1
    values = np.asarray(graph['values']).reshape(5, 5, len(col), order='F').transpose(2, 0, 1)
    return bsr_matrix((values, col, ptr), shape=(5*n, 5*n))


def arrays(factor):
    entries = [factor['buckets']] if isinstance(factor['buckets'], dict) else list(np.atleast_1d(factor['buckets']).ravel())
    plan = []
    for bucket in entries:
        pages = np.size(bucket['groupIDs'])
        if not pages:
            continue
        k = np.asarray(bucket['members']).size // pages
        members = np.asarray(bucket['members']).reshape(k, pages, order='F').T.astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape(5*k, 5*k, pages, order='F').transpose(2, 0, 1)
        plan.append((members, inverse))
    return plan


def local(plan, vector, transpose=False):
    rhs = np.asarray(vector).reshape(-1, 5)
    out = np.zeros_like(rhs)
    for members, inverse in plan:
        ng, k = members.shape
        b = rhs[members].reshape(ng, 5*k)
        x = np.einsum('pji,pj->pi' if transpose else 'pij,pj->pi', inverse, b)
        out[members] = x.reshape(ng, k, 5)
    return out.ravel()


def operator(A, plan):
    counts = dict(F=0, FT=0, M=0, MT=0)

    def forward(v):
        counts['F'] += 1
        counts['M'] += 1
        return A @ local(plan, v)

    def transpose(v):
        counts['FT'] += 1
        counts['MT'] += 1
        return local(plan, A.T @ v, transpose=True)

    return LinearOperator(A.shape, matvec=forward, rmatvec=transpose, dtype=np.float64), counts


def adjoint_error(B):
    rng = np.random.default_rng(830726)
    v, w = rng.normal(size=B.shape[1]), rng.normal(size=B.shape[0])
    forward, transpose = B @ v, B.rmatvec(w)
    left, right = float(np.dot(forward, w)), float(np.dot(v, transpose))
    return abs(left-right) / max(1., np.linalg.norm(forward)*np.linalg.norm(w), np.linalg.norm(transpose)*np.linalg.norm(v))


def main():
    folder = VAL / 'fullFlsqrscreen01'
    assert not folder.exists()
    folder.mkdir()
    bindings = {}

    def read(path, variables):
        bindings[path.relative_to(ROOT).as_posix()] = sha(path)
        return loadmat(path, simplify_cells=True, variable_names=variables)

    small = read(VAL/'overlapFunit01/restricted_overlap_unit_evidence.mat', ['cases'])
    small_reports = []
    for case in small['cases']:
        A = matrix(case['graph'])
        plan = arrays(case['plainFactor'])
        B, counts = operator(A, plan)
        adjoint = adjoint_error(B)
        assert adjoint < 1e-12
        eye = np.eye(A.shape[1])
        M = np.column_stack([local(plan, eye[:, j]) for j in range(len(eye))])
        dense = A.toarray() @ M
        probe = np.sin(np.arange(len(eye)))
        assert np.linalg.norm(B @ probe - dense @ probe) / max(1., np.linalg.norm(dense @ probe)) < 1e-12
        assert np.linalg.norm(B.rmatvec(probe) - dense.T @ probe) / max(1., np.linalg.norm(dense.T @ probe)) < 1e-12
        b = case['rhs'].ravel()
        result = lsqr(B, b, atol=1e-13, btol=1e-13, conlim=1e12, iter_lim=300)
        direct = np.linalg.lstsq(dense, b, rcond=None)[0]
        residual = float(np.linalg.norm(b-B @ result[0])/np.linalg.norm(b))
        error = float(np.linalg.norm(result[0]-direct)/max(1., np.linalg.norm(direct)))
        assert residual < 1e-9 and error < 1e-8
        zero = lsqr(B, np.zeros_like(b), iter_lim=32)
        assert np.all(zero[0] == 0) and zero[2] == 0
        small_reports.append(dict(adjointIdentityDifference=adjoint, denseSolutionRelativeDifference=error,
            fullFRelativeResidual=residual, iterations=int(result[2]), zeroRHSExactlyZero=True))
    del small, case, A, B, plan, M, dense, eye, result
    gc.collect()
    print('SMALL_FORWARD_ADJOINT_DENSE_LSQR_REFERENCES_PASSED', flush=True)

    base_file = VAL/'fullFreal01/full_F_real_evidence.mat'
    s = read(base_file, ['assemblyFile', 'rhsF', 'Fimage', 'valueF', 'sourceBinding'])
    assembly = Path(s['assemblyFile'])
    audit = json.loads((VAL/'spatialreductionreal01_independent_spatial_reduction.json').read_text())
    assert sha(assembly) == audit['inputSHA256'][str(assembly.relative_to(ROOT))]
    raw = read(assembly, ['reduction'])
    plan = arrays(raw['reduction']['Ffactor'])
    del raw
    gc.collect()
    graph_file = VAL/'fullFreal01/full_F_graph_evidence.mat'
    ga = json.loads((VAL/'fullFreal01_independent_full_F.json').read_text())
    assert sha(graph_file) == ga['inputSHA256'][str(graph_file.relative_to(ROOT))]
    graph = read(graph_file, ['graph'])['graph']
    A = matrix(graph)
    B, counts = operator(A, plan)
    adjoint = adjoint_error(B)
    assert adjoint < 1e-12
    image = (A @ s['valueF'].ravel()).reshape(s['rhsF'].shape)
    assert np.linalg.norm(image-s['Fimage'])/max(1., np.linalg.norm(s['Fimage'])) < 1e-12
    b = s['rhsF'].ravel()
    r0 = b-image.ravel()
    gradient = B.rmatvec(r0)
    gradient_image = B @ gradient
    alpha = float(np.dot(r0, gradient_image)/np.dot(gradient_image, gradient_image))
    gradient_residual = float(np.linalg.norm(r0-alpha*gradient_image)/np.linalg.norm(r0))
    # 固定32次双对角化预算；按原F右端记录最终实际剩余，不改原GMRES或NS门槛。
    initial_counts = counts.copy()
    result = lsqr(B, r0, atol=0., btol=1e-3*np.linalg.norm(b)/np.linalg.norm(r0),
        conlim=1e12, iter_lim=32)
    solver_counts = {key: counts[key]-initial_counts[key] for key in counts}
    correction = local(plan, result[0])
    final = s['valueF'].ravel()+correction
    final_image = A @ final
    relative = float(np.linalg.norm(b-final_image)/np.linalg.norm(b))
    remaining_relative = float(np.linalg.norm(r0-B @ result[0])/np.linalg.norm(r0))
    assert np.all(np.isfinite(final)) and abs(relative-remaining_relative*np.linalg.norm(r0)/np.linalg.norm(b)) < 1e-12
    savemat(folder/'full_F_lsqr_screen_evidence.mat', dict(valueF=final.reshape(s['valueF'].shape),
        Fimage=final_image.reshape(s['rhsF'].shape), rhsF=s['rhsF'], scaledCoordinates=result[0],
        sourceBinding=s['sourceBinding'], graphEvidenceSHA256=sha(graph_file), originalFactorAssemblySHA256=sha(assembly)), do_compression=True)
    assert all(sha(ROOT/name) == digest for name, digest in bindings.items())
    out = dict(status='tested-independent-full-F-forward-adjoint-and-bounded-right-preconditioned-LSQR-screen',
        processID=os.getpid(), recordedUTC=datetime.now(timezone.utc).isoformat(), smallReferences=small_reports,
        method='LSQR on exact A_FF times certified original within-patch inverse; adjoint is M_F^T A_FF^T',
        realAdjointIdentityDifference=adjoint, maximumBidiagonalizationSteps=32,
        originalSavedFullFResidual=float(np.linalg.norm(r0)/np.linalg.norm(b)),
        optimalAdjointGradientStepOnSavedRemaining=gradient_residual, finalRelativeFullFResidual=relative,
        relativeCorrectionEquationResidual=remaining_relative, fullFInternalTargetMet=relative <= 1e-3,
        nativeStopCode=int(result[1]), actualIterations=int(result[2]),
        nativeResidualNorm=float(result[3]), nativeNormalResidualNorm=float(result[7]),
        nativeEstimatedCondition=float(result[6]), solverApplications=solver_counts,
        totalCompositeApplications=counts, separateOriginalFImageEvaluations=2,
        separateFinalLocalMapApplications=1, newLocalInverseConstructions=0,
        nativeGMRESArnoldiColumnChecksClaimed=False, independentAllRealBidiagonalizationColumnsClaimed=False,
        originalCOrSchurIdentityClaimed=False, matlabProductionSolverModified=False,
        physicalStateChanged=False, newSelectedJApplications=0, newActualNSEvaluations=0,
        newFlowAcceptances=0, newProductionSolvedAreaM2=0, inputSHA256=bindings,
        processPeakMemoryClaimed=False, scriptSHA256=sha(__file__))
    (folder/'full_F_lsqr_screen_summary.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k: v for k, v in out.items() if k != 'inputSHA256'}), flush=True)


if __name__ == '__main__':
    main()
