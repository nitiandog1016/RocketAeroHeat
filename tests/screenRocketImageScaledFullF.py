# -*- coding: utf-8 -*-
"""逐组完整计算A_FF*M_F像列范数，右侧归一后作同端点/同预算LSQR筛查。"""
from pathlib import Path
from datetime import datetime, timezone
import os
import json
import gc
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse.linalg import LinearOperator, lsqr
from screenRocketFullFLeastSquares import ROOT, VAL, sha, arrays, local, matrix, operator, adjoint_error


def column_norms(A, plan):
    """模块一：每次仅形成一组完整像列；所有输入块按列所属组恰读取一次。"""
    n = A.shape[0]//5
    ptr, col = A.indptr, A.indices
    owner = np.repeat(np.arange(n), np.diff(ptr))
    entries = np.argsort(col, kind='stable')
    counts = np.bincount(col, minlength=n)
    colptr = np.r_[0, np.cumsum(counts)]
    result = np.empty((n, 5))
    coverage = np.zeros(n, np.int64)
    used = groups = maximum_panel_bytes = maximum_support = 0
    for nodes, inverses in plan:
        for ids, inverse in zip(nodes, inverses):
            blocks = np.concatenate([entries[colptr[node]:colptr[node+1]] for node in ids])
            rows, local_rows = np.unique(owner[blocks], return_inverse=True)
            local_cols = np.repeat(np.arange(len(ids)), counts[ids])
            panel = np.zeros((5*len(rows), 5*len(ids)))
            rr = 5*local_rows[:, None]+np.arange(5)[None, :]
            cc = 5*local_cols[:, None]+np.arange(5)[None, :]
            panel[rr[:, :, None], cc[:, None, :]] = A.data[blocks]
            image = panel @ inverse
            norms = np.linalg.norm(image, axis=0)
            assert np.all(np.isfinite(norms) & (norms > 0))
            result[ids] = norms.reshape(len(ids), 5)
            coverage[ids] += 1
            groups += 1
            used += len(blocks)
            maximum_panel_bytes = max(maximum_panel_bytes, panel.nbytes+image.nbytes)
            maximum_support = max(maximum_support, len(rows))
            if groups % 5000 == 0:
                print('COMPLETE_IMAGE_NORM_GROUPS='+str(groups), flush=True)
    assert np.all(coverage == 1) and used == len(col)
    return result, dict(groups=groups, originalBlocksUsedExactlyOnce=used,
        allOriginalUnknownColumnsCoveredExactlyOnce=True, maximumSupportedOriginalRows=maximum_support,
        maximumDensePanelAndImageArrayBytes=maximum_panel_bytes,
        fullGlobalPreconditionedMatrixAllocated=False, processPeakMemoryClaimed=False)


def scaled_operator(B, scale):
    d = scale.ravel()
    return LinearOperator(B.shape, matvec=lambda v: B @ (d*v),
        rmatvec=lambda w: d*B.rmatvec(w), dtype=np.float64)


def main():
    folder = VAL/'fullFimagescalescreen01'
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
        norms, construction = column_norms(A, plan)
        eye = np.eye(A.shape[1])
        dense = np.column_stack([B @ eye[:, j] for j in range(len(eye))])
        expected = np.linalg.norm(dense, axis=0).reshape(norms.shape)
        error = np.max(abs(norms-expected)/expected)
        assert error < 1e-12
        scale = 1/norms
        S = scaled_operator(B, scale)
        adjoint = adjoint_error(S)
        assert adjoint < 1e-12
        b = case['rhs'].ravel()
        result = lsqr(S, b, atol=0., btol=1e-11, conlim=1e12, iter_lim=300)
        residual = float(np.linalg.norm(b-S @ result[0])/np.linalg.norm(b))
        assert residual < 1e-10
        small_reports.append(dict(columnNormMaximumRelativeDifference=float(error), adjointIdentityDifference=adjoint,
            fullFRelativeResidual=residual, construction=construction))
    del small, case, A, plan, B, dense, eye, S, result
    gc.collect()
    print('SMALL_COMPLETE_IMAGE_NORMS_AND_SCALED_ADJOINT_PASSED', flush=True)

    source = read(VAL/'fullFreal01/full_F_real_evidence.mat', ['assemblyFile','rhsF','valueF','Fimage','sourceBinding'])
    assembly_file = Path(source['assemblyFile'])
    old_audit = json.loads((VAL/'spatialreductionreal01_independent_spatial_reduction.json').read_text())
    assert sha(assembly_file) == old_audit['inputSHA256'][str(assembly_file.relative_to(ROOT))]
    raw = read(assembly_file, ['reduction'])
    plan = arrays(raw['reduction']['Ffactor'])
    del raw
    gc.collect()
    graph_file = VAL/'fullFreal01/full_F_graph_evidence.mat'
    audit = json.loads((VAL/'fullFreal01_independent_full_F.json').read_text())
    assert sha(graph_file) == audit['inputSHA256'][str(graph_file.relative_to(ROOT))]
    graph = read(graph_file, ['graph'])['graph']
    A = matrix(graph)
    B, counts = operator(A, plan)
    norms, construction = column_norms(A, plan)
    assert construction['groups'] == 28902 and construction['originalBlocksUsedExactlyOnce'] == 1285801
    scale = 1/norms
    assert np.all(np.isfinite(scale) & (scale > 0))
    S = scaled_operator(B, scale)
    adjoint = adjoint_error(S)
    assert adjoint < 1e-12
    # 独立完整作用检查极端列及各分位数，构造没有使用这些全向量探针。
    ranking = np.argsort(norms.ravel())
    probes = np.unique(ranking[np.linspace(0, len(ranking)-1, 17).astype(np.int64)])
    maximum_probe_error = 0.
    for index in probes:
        v = np.zeros(A.shape[1]);v[index] = 1
        actual = np.linalg.norm(B @ v)
        error = abs(actual-norms.ravel()[index])/norms.ravel()[index]
        assert error < 1e-12
        maximum_probe_error = max(maximum_probe_error, float(error))
    scale_file = folder/'full_F_image_scale_evidence.mat'
    savemat(scale_file, dict(columnNorm=norms, columnScale=scale, sourceBinding=source['sourceBinding'],
        graphEvidenceSHA256=sha(graph_file), originalFactorAssemblySHA256=sha(assembly_file)), do_compression=True)
    b = source['rhsF'].ravel()
    image0 = A @ source['valueF'].ravel()
    assert np.linalg.norm(image0-source['Fimage'].ravel())/max(1., np.linalg.norm(image0)) < 1e-12
    residual0 = b-image0
    before = counts.copy()
    result = lsqr(S, residual0, atol=0., btol=1e-3*np.linalg.norm(b)/np.linalg.norm(residual0), conlim=1e12, iter_lim=32)
    native_counts = {key: counts[key]-before[key] for key in counts}
    correction = local(plan, scale.ravel()*result[0])
    final = source['valueF'].ravel()+correction
    image = A @ final
    final_relative = float(np.linalg.norm(b-image)/np.linalg.norm(b))
    predicted = float(np.linalg.norm(residual0-S @ result[0])/np.linalg.norm(b))
    assert abs(final_relative-predicted) < 1e-12
    unscaled_file = VAL/'fullFlsqrscreen01/full_F_lsqr_screen_summary.json'
    unscaled = json.loads(unscaled_file.read_text())
    bindings[unscaled_file.relative_to(ROOT).as_posix()] = sha(unscaled_file)
    savemat(folder/'full_F_image_scaled_screen_evidence.mat', dict(rhsF=source['rhsF'], valueF=final.reshape(source['rhsF'].shape),
        Fimage=image.reshape(source['rhsF'].shape), initial=source['valueF'], sourceBinding=source['sourceBinding']), do_compression=True)
    assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
    out = dict(status='tested-complete-preconditioned-F-image-column-normalization-and-bounded-LSQR-screen',
        processID=os.getpid(), recordedUTC=datetime.now(timezone.utc).isoformat(), smallReferences=small_reports,
        construction=construction, columnNormRange=[float(norms.min()),float(norms.max())],
        columnNormSpan=float(norms.max()/norms.min()), columnNormQuantiles=np.quantile(norms,[0,.25,.5,.75,1]).tolist(),
        columnScaleRange=[float(scale.min()),float(scale.max())], columnScaleFile=scale_file.relative_to(ROOT).as_posix(), columnScaleSHA256=sha(scale_file),
        realAdjointIdentityDifference=adjoint, completeColumnProbes=len(probes), maximumColumnProbeRelativeDifference=maximum_probe_error,
        originalSavedFullFResidual=float(np.linalg.norm(residual0)/np.linalg.norm(b)), maximumBidiagonalizationSteps=32,
        actualIterations=int(result[2]), nativeStopCode=int(result[1]), finalRelativeFullFResidual=final_relative,
        unscaledSameEndpointSameIterationBudgetResidual=unscaled['finalRelativeFullFResidual'], fullFInternalTargetMet=final_relative<=1e-3,
        nativeEstimatedCondition=float(result[6]), nativeSolverApplications=native_counts, totalCompositeApplications=counts,
        separateOriginalFImageEvaluations=2, separateFinalLocalMapApplications=1,
        fixedPositiveRightScalingOnly=True, originalFResidualNormAndRightHandSideUnchanged=True,
        oldApproximateSchurIdentityClaimed=False, physicalStateChanged=False, matlabProductionSolverModified=False,
        independentAllRealBidiagonalizationColumnsClaimed=False, newSelectedJApplications=0, newActualNSEvaluations=0,
        newFlowAcceptances=0, newProductionSolvedAreaM2=0, inputSHA256=bindings,
        helperSourceSHA256=sha(ROOT/'tests/screenRocketFullFLeastSquares.py'), scriptSHA256=sha(__file__))
    (folder/'full_F_image_scaled_screen_summary.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','smallReferences')}), flush=True)


if __name__ == '__main__':
    main()
