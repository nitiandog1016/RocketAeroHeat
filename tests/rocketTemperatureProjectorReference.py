"""读取完整或紧凑QR坐标基，以原物理排列独立计算投影。"""
import numpy as np


def project_coordinates(projector, coordinates):
    value = np.asarray(coordinates).ravel(order='F').copy()
    if 'basisCoordinateIndices' in projector:
        ids = np.atleast_1d(projector['basisCoordinateIndices']).astype(int)-1
        assert len(np.unique(ids)) == len(ids) and np.all((ids >= 0) & (ids < len(value)))
    else:
        ids = np.arange(len(value))
    U = np.asarray(projector['basis']).reshape(len(ids), -1, order='F')
    np.testing.assert_allclose(U.T@U, np.eye(U.shape[1]), atol=1e-12, rtol=0)
    value[ids] -= U@(U.T@value[ids])
    return value.reshape(np.shape(coordinates), order='F')


def certify_projector_span(projector, H, transpose_right_map):
    """逐个独立原行计算完整P转置，不把紧凑支撑正确性作为输入假设。"""
    reps = np.atleast_1d(projector['exactRepresentativeRowIndices']).astype(int)-1
    mapping = np.atleast_1d(projector['originalToUniqueRow']).astype(int)-1
    assert (H-H[reps[mapping]]).nnz == 0
    ids = (np.atleast_1d(projector['basisCoordinateIndices']).astype(int)-1
           if 'basisCoordinateIndices' in projector else np.arange(H.shape[1]))
    U = np.asarray(projector['basis']).reshape(len(ids), len(reps), order='F')
    np.testing.assert_allclose(U.T@U, np.eye(len(reps)), atol=1e-12, rtol=0)
    outside = np.ones(H.shape[1], dtype=bool); outside[ids] = False
    errors = []; norms = []
    for row in reps:
        column = transpose_right_map(H.getrow(row).toarray().ravel())
        assert np.all(column[outside] == 0), '完整P转置H转置在声称支撑外不为零'
        norm = np.linalg.norm(column); assert np.isfinite(norm) and norm > 0
        error = float(np.linalg.norm(column[ids]-U@(U.T@column[ids]))/norm)
        assert error < 1e-10
        errors.append(error); norms.append(norm)
    np.testing.assert_allclose(np.asarray(norms)[mapping], np.atleast_1d(projector['constraintColumnNorms']), rtol=1e-12)
    if 'basisCoordinateIndices' in projector:
        assert U.nbytes == projector['basisArrayBytes']
        assert projector['additionalConstructionArrayBoundBytes'] <= projector['originalConstructionBudgetBytes'] == 448*1024**2
    return dict(originalRows=H.shape[0], distinctConditions=len(reps), storedCoordinateCount=len(ids),
                compactStorage='basisCoordinateIndices' in projector,
                maximumCompleteRightMapSpanDifference=max(errors), allOmittedCoordinatesExactlyZero=True)
