# -*- coding: utf-8 -*-
"""独立复核旧外缘拒绝及首个局部域原面集合；不生成计算结果。"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(out):
    """模块一：以Heron面积及原MAT标志，独立核对既有CSV传播诊断。"""
    started = time.monotonic()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    mat = ROOT/'outputs/v104_validation/baseline02/current_case.mat'
    csv_path = ROOT/'outputs/v189_validation/analysis02/original_face_diagnosis.csv'
    protection = ROOT/'outputs/v176_validation/protected_initial.json'
    inputs = [mat, csv_path, protection, Path(__file__), ROOT/'solveExternalFlow.m',
              ROOT/'solveThermallyPerfectPrandtlMeyerExpansion.m']
    before = {str(p): digest(p) for p in inputs}
    result = dict(passed=False, newSolvedAreaM2=0, newPhysicalQualification=False)
    try:
        m = loadmat(mat, simplify_cells=True)['current']
        v = m['V']; f = m['F'].astype(int)-1; n = len(f)
        xyz = v[f]; center = xyz.mean(axis=1)
        lengths = np.linalg.norm(xyz-np.roll(xyz, 1, axis=1), axis=2)
        half = lengths.sum(axis=1)/2
        area = np.sqrt(half*np.prod(half[:, None]-lengths, axis=1))
        with csv_path.open(encoding='utf-8', newline='') as stream:
            rows = list(csv.DictReader(stream))
        original_area = np.array([float(r['areaM2']) for r in rows])
        assert len(rows) == n == 19752
        assert np.array_equal([int(r['faceID']) for r in rows], np.arange(n)+1)
        error = float(np.max(np.abs(area-original_area)))
        assert error < 1e-14
        assert abs(area.sum()-1.9646292565715551) < 2e-13
        roots = np.array([int(r['directStateRootFace']) for r in rows])
        edge = np.array([int(r['edgePrimaryCandidate']) for r in rows], dtype=bool)
        thermal = np.asarray(m['externalFlow']['edge']['thermalExpansionOutOfRangeFace'], dtype=bool)
        direct = np.asarray(m['continuumRegime']['directStateUnresolvedFace'], dtype=bool)
        root_mask = np.zeros(n, dtype=bool)
        root_mask[np.unique(roots[edge])-1] = True
        assert np.all(roots[edge] > 0)
        assert np.array_equal(root_mask, thermal)
        assert np.all(direct[thermal])
        assert np.array_equal(edge, np.isin(m['heatingCoverage']['reasonCode'], [4, 14]))
        predecessor = m['surfaceStreamlines']['predecessorFace'].astype(int)
        # 独立验证CSV根的每条前驱关系，不重复原深度优先搜索。
        for i in np.flatnonzero(edge):
            if thermal[i]:
                assert roots[i] == i+1
            else:
                parent = predecessor[i]-1
                assert 0 <= parent < n and edge[parent] and roots[parent] == roots[i]
        angle = m['externalFlow']['edge']['incidenceAngleDeg']
        valid_expansion = (~thermal) & (angle < 0) & m['externalFlow']['edge']['thermalStateValidFace'].astype(bool)
        def describe(mask):
            ids = np.flatnonzero(mask)
            return dict(faceCount=int(len(ids)), areaM2=float(area[ids].sum()),
                        faceIDs=(ids+1).tolist(),
                        centerXRangeM=[float(center[ids, 0].min()), float(center[ids, 0].max())])
        result['thermalRoots'] = describe(thermal)
        result['thermalRoots']['incidenceAngleRangeDeg'] = [float(angle[thermal].min()), float(angle[thermal].max())]
        result['mostNegativeAcceptedExpansionAngleDeg'] = float(angle[valid_expansion].min())
        result['rootIdentity'] = 'All 336 distinct traced roots equal the saved thermalExpansionOutOfRangeFace mask.'
        result['interpretation'] = '80 K is the old local simple-wave model limit, not a measured global-flow temperature.'
        result['independentHeronMaxAreaDifferenceM2'] = error
        # 模块二：在既有完整环处选域。包含鼻部，避免凭theta/H伪造入口剖面。
        rings = np.unique(v[:, 0])
        cut = float(rings[np.argmin(abs(rings+0.9))])
        extension = float(rings[np.argmin(abs(rings+1.2))])
        trusted = float(rings[np.argmin(abs(rings+0.75))])
        local = xyz[:, :, 0].min(axis=1) >= cut-1e-10
        extended = xyz[:, :, 0].min(axis=1) >= extension-1e-10
        trusted_mask = xyz[:, :, 0].min(axis=1) >= trusted-1e-10
        crossing = (xyz[:, :, 0].min(axis=1) < cut-1e-10) & (xyz[:, :, 0].max(axis=1) > cut+1e-10)
        assert not crossing.any() and np.all(local[thermal]) and np.all(local[trusted_mask])
        # 切口是流动出口，绝不生成虚假物理端盖；全周向原壁面保持。
        geometry = dict(schema='RAH.HybridRegionGeometry.v1', status='design-only-not-admitted',
                        cutXM=cut, extensionTestCutXM=extension, trustedEndXM=trusted,
                        cutRingVertices=int(np.sum(abs(v[:, 0]-cut)<1e-10)),
                        originalAreaM2=float(area.sum()), fullCircumference=True,
                        localWall=describe(local), extendedWall=describe(extended),
                        provisionalTrustedWall=describe(trusted_mask),
                        rootWall=describe(thermal), crossingOriginalFaces=int(crossing.sum()),
                        otherOriginalFacesRetained=True, artificialPhysicalCap=False,
                        newSolvedAreaM2=0, newPhysicalQualification=False)
        (out/'region_geometry.json').write_text(json.dumps(geometry, indent=2)+'\n', encoding='utf-8')
        result['candidateGeometry'] = {k:geometry[k] for k in ('cutXM','extensionTestCutXM','trustedEndXM','cutRingVertices','crossingOriginalFaces')}
        result['candidateGeometry'].update(localFaces=int(local.sum()), localAreaM2=float(area[local].sum()),
                                          extendedFaces=int(extended.sum()), extendedAreaM2=float(area[extended].sum()))
        # 模块三：验证并发保护文件及本轮输入没有被修改。
        protected = json.loads(protection.read_text(encoding='utf-8'))
        changed = [name for name, expected in protected.items() if digest(ROOT/name) != expected]
        result['protectedFiles'] = len(protected); result['changedProtectedFiles'] = changed
        assert not changed
        result['sourceSHA256'] = before
        result['inputStable'] = all(digest(p) == h for p, h in before.items())
        assert result['inputStable']
        result['passed'] = True
    except Exception as exc:
        result['error'] = repr(exc)
        raise
    finally:
        result['elapsedSeconds'] = time.monotonic()-started
        (out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('passed','candidateGeometry','protectedFiles','elapsedSeconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    verify(parser.parse_args().out)
