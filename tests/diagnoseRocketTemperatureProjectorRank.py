"""独立重建失败检查点的完整支撑右坐标列，诊断拒绝原因；不删除或截断约束。"""
from pathlib import Path
import os
import sys
import json
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from riemannReservoirReference import ROOT, sha


def main(step_folder, output_folder):
    folder = Path(step_folder).resolve(); out = Path(output_folder).resolve()
    preview_file = folder/'linear_checkpoints/repreview00.mat'
    context_file = folder/'direction_context.mat'
    bindings = {p.relative_to(ROOT).as_posix(): sha(p) for p in (preview_file, context_file)}
    saved = loadmat(preview_file, simplify_cells=True)
    rows = saved['rows']; selected = np.atleast_1d(saved['summary']['newActiveFaceIDs']).astype(int)
    ids = np.atleast_1d(rows['faceIDs']).astype(int); loc = np.searchsorted(ids, selected)
    np.testing.assert_array_equal(ids[loc], selected)
    H = rows['H'][loc].tocsr(); n = int(rows['cellCount'])
    seen = {}; representatives = []; mapping = []
    for j in range(H.shape[0]):
        h = H.getrow(j); key = (tuple(h.indices), h.data.tobytes())
        if key not in seen:
            seen[key] = len(representatives); representatives.append(j)
        mapping.append(seen[key])
    assert (H-H[np.asarray(representatives)[mapping]]).nnz == 0
    factor = loadmat(context_file, simplify_cells=True, variable_names=['factor'])['factor']
    buckets = factor['buckets'] if isinstance(factor['buckets'], (list, np.ndarray)) else [factor['buckets']]
    touched = np.zeros(n, bool); touched[H.indices % n] = True
    covered = np.zeros(n, bool); pages = []
    for b in buckets:
        size = np.asarray(b['map']).shape[0]//5
        members = np.asarray(b['members']).astype(int).reshape(size, -1, order='F')-1
        b['members0'] = members
        selected_pages = np.flatnonzero(np.any(touched[members], axis=0))
        pages.append(selected_pages)
        assert not np.any(covered[members[:, selected_pages].ravel()])
        covered[members[:, selected_pages].ravel()] = True
    assert np.all(covered[touched])
    support = np.flatnonzero(covered); lookup = np.full(n, -1, int)
    lookup[support] = np.arange(len(support))
    columns = np.zeros((5*len(support), len(representatives)))
    for b, selected_pages in zip(buckets, pages):
        for page in selected_pages:
            members = b['members0'][:, page]
            original = (members[:, None]+np.arange(5)*n).ravel()
            compact = (lookup[members, None]+np.arange(5)*len(support)).ravel()
            matrix = np.asarray(b['map']); matrix = matrix if matrix.ndim == 2 else matrix[:, :, page]
            columns[compact] = matrix.T @ H[representatives][:, original].toarray().T
    norms = np.linalg.norm(columns, axis=0); assert np.all(norms > 0)
    unit = columns/norms
    _, R = np.linalg.qr(unit, mode='reduced'); diagonal = abs(np.diag(R))
    failed = np.flatnonzero(diagonal <= 1e-10); singular = np.linalg.svd(unit, compute_uv=False)
    failures = []
    for j in failed:
        i = representatives[j]; h = H.getrow(i); matches = []
        for earlier in representatives[:j]:
            a = H.getrow(earlier)
            if np.array_equal(h.indices, a.indices):
                ratio = h.data[0]/a.data[0]
                error = float(np.max(abs(h.data-ratio*a.data))/max(np.max(abs(h.data)), 1e-300))
                matches.append(dict(earlierFace=int(selected[earlier]), ratio=float(ratio),
                                    relativeDifference=error, bitwiseScalarReconstruction=bool(np.array_equal(h.data, ratio*a.data))))
        failures.append(dict(representative=int(j+1), faceID=int(selected[i]), pivot=float(diagonal[j]),
                             ownerCV=int(np.atleast_1d(rows['ownerCV'])[loc[i]]),
                             neighbourCV=int(np.atleast_1d(rows['neighbourCV'])[loc[i]]),
                             region=int(np.atleast_1d(rows['sampleRegion'])[loc[i]]),
                             sameSupportComparisons=sorted(matches, key=lambda x: x['relativeDifference'])[:3]))
    assert len(failed), '独立QR应重现本次原门槛拒绝，不能默认失败归因。'
    assert all(sha(ROOT/p) == h for p, h in bindings.items())
    result = dict(status='independent-Q40-compact-projector-rank-failure-diagnostic',
                  activeOriginalRows=H.shape[0], exactDistinctRows=len(representatives), supportCV=len(support),
                  minimumNormalizedQRDiagonal=float(diagonal.min()), failedPivots=failures,
                  singularValuesLast10=singular[-10:].tolist(), originalRankGateUnchanged=1e-10,
                  noRankDeletionOrNewFlowUpdate=True, inputSHA256=bindings, scriptSHA256=sha(__file__),
                  newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0)
    out.mkdir(parents=True, exist_ok=True)
    target = out/'projector_rank_diagnostic.json'; assert not target.exists()
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    savemat(out/'projector_rank_reference.mat', dict(unitColumns=unit, normalizedQRDiagonal=diagonal,
            selectedFaces=selected, representatives=np.asarray(representatives)+1, supportCells=support+1))
    print(json.dumps(result))


if __name__ == '__main__':
    main(*sys.argv[1:])
