# -*- coding: utf-8 -*-
"""在真实当前系统上比较完整增广像的分块QR右预处理，不更新流场。"""
from pathlib import Path
from datetime import datetime, timezone
import os, sys, json, gc, time
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import bsr_matrix
from scipy.sparse.linalg import LinearOperator, lsqr
from screenRocketFullFLeastSquares import ROOT, sha, local
from screenRocketParentImageScale import setup


def metric_blocks(metric, order):
    """独立按已保存链式系数构造每个CV的5×5度量，含冻结行权重。"""
    n = len(order)
    W = np.zeros((n, 5, 5))
    for j in range(5):
        W[:, j, j] = metric['diagonal'][order, j]
    W[:, 1:4, 0] = metric['densityToVelocity'][order]
    W[:, 4, 0] = metric['densityToTemperature'][order]
    W[:, 4, 1:4] = metric['momentumToTemperature'][order]
    # reweight已将权重乘进上述系数，peakRowWeights只用于报告原始幅度。
    return W


def block_weight(W, x, transpose=False):
    return np.einsum('pji,pj->pi' if transpose else 'pij,pj->pi', W, x.reshape(-1, 5)).ravel()


def build_qr(A, plan, W, penalty):
    """只在每组至多40个原未知量上正交化；所有原C行均进入该组像。"""
    n = len(W)
    counts = np.bincount(A.indices, minlength=n)
    by_column = np.argsort(A.indices, kind='stable')
    pointer = np.r_[0, np.cumsum(counts)]
    owners = np.repeat(np.arange(n), np.diff(A.indptr))
    new_plan = []
    used = 0
    max_error = 0.
    min_diagonal = 1.
    old_local_min_ratio = 1.
    started = time.perf_counter()
    for members, old_inverse in plan:
        ng, k = members.shape
        inverse = np.empty((ng, 5*k, 5*k))
        for page, group in enumerate(members):
            ids = np.concatenate([by_column[pointer[j]:pointer[j+1]] for j in group])
            cols = np.concatenate([np.full(counts[j], a) for a, j in enumerate(group)])
            rows, local_rows = np.unique(owners[ids], return_inverse=True)
            panel = np.zeros((len(rows), 5, k, 5))
            panel[local_rows, :, cols, :] = A.data[ids]
            panel = panel.reshape(5*len(rows), 5*k)
            bottom = np.zeros((k, 5, k, 5))
            for j, cell in enumerate(group):
                bottom[j, :, j, :] = penalty*W[cell]
            augmented = np.row_stack([panel, bottom.reshape(5*k, 5*k)])
            norms = np.linalg.norm(augmented, axis=0)
            assert np.all(np.isfinite(norms)) and np.all(norms > 0)
            # 先归一再QR；不构造正规方程，也不丢弃小奇异值或补对角。
            _, R = np.linalg.qr(augmented/norms, mode='reduced')
            P = np.linalg.solve(R, np.eye(5*k))/norms[:, None]
            image = augmented@P
            error = np.linalg.norm(image.T@image-np.eye(5*k), ord=np.inf)
            assert np.isfinite(error) and error < 1e-8
            inverse[page] = P
            max_error = max(max_error, float(error))
            min_diagonal = min(min_diagonal, float(np.min(abs(np.diag(R)))))
            # 旧组像列的最小/最大奇异值只作局部诊断，不代表全局条件数。
            old = augmented@old_inverse[page]
            old /= np.linalg.norm(old, axis=0)
            sv = np.linalg.svd(old, compute_uv=False)
            old_local_min_ratio = min(old_local_min_ratio, float(sv[-1]/sv[0]))
            used += len(ids)
        new_plan.append((members, inverse))
    assert used == len(A.indices)
    return new_plan, dict(originalBlocksUsedExactlyOnce=used,
                         maximumAugmentedLocalGramDefect=max_error,
                         minimumNormalizedQRDiagonal=min_diagonal,
                         oldMinimumWithinPatchImageSingularRatio=old_local_min_ratio,
                         constructionSeconds=time.perf_counter()-started,
                         globalConditionNumberClaimed=False)


def augmented_operator(A, AT, plan, W, penalty, scale=None):
    d = np.ones(A.shape[1]) if scale is None else scale.ravel()
    counts = dict(forward=0, transpose=0)
    def forward(x):
        counts['forward'] += 1
        v = local(plan, d*x)
        return np.r_[A@v, penalty*block_weight(W, v)]
    def transpose(x):
        counts['transpose'] += 1
        back = AT@x[:A.shape[0]]+penalty*block_weight(W, x[A.shape[0]:], True)
        return d*local(plan, back, True)
    return LinearOperator((2*A.shape[0], A.shape[1]), matvec=forward, rmatvec=transpose, dtype=float), counts


def adjoint(B):
    x = np.sin(np.arange(B.shape[1])*.17)
    y = np.cos(np.arange(B.shape[0])*.23)
    bx, bty = B@x, B.rmatvec(y)
    return abs(float(bx@y-x@bty))/max(1., np.linalg.norm(bx)*np.linalg.norm(y), np.linalg.norm(bty)*np.linalg.norm(x))


def small_reference():
    packed_file = ROOT/'outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat'
    metric_file = ROOT/'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat'
    packed_cases = loadmat(packed_file, simplify_cells=True, variable_names=['cases'])['cases']
    metric_cases = loadmat(metric_file, simplify_cells=True, variable_names=['cases'])['cases']
    checks = []
    for a, m in zip(packed_cases, metric_cases):
        A, plan, order = setup(a['packed'], a['factor'])
        W = metric_blocks(m['metric'], order)
        n = len(order);N = 5*n
        # 稠密参考使用独立MATLAB已保存的完整W，跨两种物理/代数排列核对。
        ids = (order[:, None]+n*np.arange(5)[None, :]).ravel()
        expected = m['W'][np.ix_(ids, ids)]
        actual = np.column_stack([block_weight(W, np.eye(N)[:, j]) for j in range(N)])
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-14)
        rhs = m['base']['rhs'][order].ravel()
        penalty = m['detail']['penaltyInOriginalResidualUnits']
        qplan, construction = build_qr(A, plan, W, penalty)
        B, _ = augmented_operator(A, A.T.tobsr(), qplan, W, penalty)
        adj = adjoint(B);assert adj < 1e-12
        augmented = np.row_stack([A.toarray(), penalty*expected])
        result = lsqr(B, np.r_[rhs, np.zeros(N)], atol=1e-14, btol=1e-14, conlim=1e14, iter_lim=900)
        value = local(qplan, result[0])
        dense = np.linalg.lstsq(augmented, np.r_[rhs, np.zeros(N)], rcond=None)[0]
        error = float(np.linalg.norm(value-dense)/max(1., np.linalg.norm(dense)))
        assert error < 1e-8
        checks.append(dict(denseSolutionRelativeDifference=error, adjointRelativeError=adj, construction=construction))
    return checks, {str(packed_file): sha(packed_file), str(metric_file): sha(metric_file)}


def main(folder, preparation):
    folder, preparation = Path(folder).resolve(), Path(preparation).resolve()
    assert not folder.exists()
    folder.mkdir(parents=True)
    p = json.loads(preparation.read_text(encoding='utf-8-sig'))
    small_checks, bindings = small_reference()
    print('SMALL_REFERENCE', json.dumps(small_checks), flush=True)
    bindings[str(preparation)] = sha(preparation)
    def read(path, variables):
        path = Path(path)
        bindings[str(path)] = sha(path)
        return loadmat(path, simplify_cells=True, variable_names=variables)
    for a in p['artifacts'].values():
        assert sha(a['file']) == a['sha256']
    assert sha(p['sourceBinding']['sourceFile']) == p['sourceBinding']['sourceFileSHA256']
    for identity in ('methodIdentity', 'sourceIdentity'):
        for name, digest in p[identity].items():
            path = ROOT/'internalFlow/cutCell'/(name+'.m')
            if name == 'resolveNSCharacteristicBoundary3D':
                path = ROOT/'internalFlow'/(name+'.m')
            assert sha(path) == digest
            bindings[str(path)] = digest
    raw = read(p['artifacts']['packed']['file'], ['packed'])
    factor = read(p['artifacts']['factor']['file'], ['factor'])
    A, plan, order = setup(raw['packed'], factor['factor'])
    del raw, factor
    gc.collect()
    m = read(p['artifacts']['metric']['file'], ['metric', 'rhs', 'metricNormalization'])
    W = metric_blocks(m['metric'], order)
    rhs = m['rhs'][order].ravel()
    penalty = .1*np.linalg.norm(rhs)/p['metricNormalization']
    column_scale = read(p['artifacts']['scale']['file'], ['columnScale'])['columnScale'][order]
    qplan, construction = build_qr(A, plan, W, penalty)
    print('QR_CONSTRUCTION', json.dumps(construction), flush=True)
    AT = A.T.tobsr()
    reports = []
    for label, mapping, d in [('old-patch-inverse-image-column-scale', plan, column_scale),
                              ('full-augmented-image-patch-QR', qplan, None)]:
        B, counts = augmented_operator(A, AT, mapping, W, penalty, d)
        adj = adjoint(B)
        assert adj < 1e-12
        before = counts.copy()
        start = time.perf_counter()
        result = lsqr(B, np.r_[rhs, np.zeros_like(rhs)], atol=1e-12, btol=1e-12, conlim=1e14, iter_lim=128)
        elapsed = time.perf_counter()-start
        coordinates = result[0]
        v = local(mapping, coordinates if d is None else d.ravel()*coordinates)
        image = A@v
        weighted = block_weight(W, v)
        relative = float(np.linalg.norm(rhs-image)/np.linalg.norm(rhs))
        objective = relative**2+(penalty*np.linalg.norm(weighted)/np.linalg.norm(rhs))**2
        report = dict(method=label, originalCRelativeResidual=relative, normalizedAugmentedObjective=float(objective),
                      nativeIterations=int(result[2]), nativeStopCode=int(result[1]),
                      actualOperatorCalls={k: counts[k]-before[k] for k in counts},
                      elapsedSeconds=elapsed, adjointRelativeError=adj, originalCTargetMet=relative <= 1e-3)
        reports.append(report)
        print(json.dumps(report), flush=True)
        physical = np.empty((len(order), 5));physical[order] = v.reshape(-1, 5)
        physical_image = np.empty_like(physical);physical_image[order] = image.reshape(-1, 5)
        savemat(folder/(label+'_evidence.mat'), dict(value=physical, image=physical_image,
                rhs=m['rhs'], coordinatesAlgebraic=coordinates, report=report,
                preparationFile=str(preparation), preparationSHA256=sha(preparation)), do_compression=True)
    assert all(sha(path) == digest for path, digest in bindings.items())
    out = dict(status='tested-independent-current-system-same-objective-128-step-right-preconditioning-screen',
               sourceBinding=p['sourceBinding'], preparationFile=str(preparation), preparationSHA256=sha(preparation),
               originalCVs=len(order), unknowns=5*len(order), regularizationStrength=.1, penalty=penalty,
               allOriginalEquationsKept=True, smallReferences=small_checks, construction=construction, comparisons=reports,
               newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
               physicalSolution=False, processID=os.getpid(), recordedUTC=datetime.now(timezone.utc).isoformat(),
               inputSHA256=bindings, scriptSHA256=sha(__file__))
    (folder/'qr_screen.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main(*sys.argv[1:3])
