# -*- coding: utf-8 -*-
"""用独立平面交点和三角重心方程，核对真实箭体最近内法线命中。"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.io import loadmat

from auditThermalPressureSurfaceFlow import sha


def audit(inlet_file, mapping_file, process_file, report_file):
    assert not report_file.exists()
    source = loadmat(inlet_file, simplify_cells=True); saved = loadmat(mapping_file, simplify_cells=True)
    process = json.loads(process_file.read_text(encoding='utf-8-sig'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
    V = source['geometry']['V']; F = source['geometry']['F'].astype(int)-1; m = saved['mapping']
    triangles = V[F]; a = triangles[:, 1]-triangles[:, 0]; b = triangles[:, 2]-triangles[:, 0]
    N = np.cross(a, b); plane = np.einsum('ij,ij->i', N, triangles[:, 0])
    aa = np.einsum('ij,ij->i', a, a); bb = np.einsum('ij,ij->i', b, b); ab = np.einsum('ij,ij->i', a, b)
    determinant = aa*bb-ab*ab
    errors = []
    for j, (origin, direction) in enumerate(zip(m['rayPosition'], m['rayDirection'])):
        dot = N@direction
        t = np.divide(plane-N@origin, dot, out=np.full(len(F), np.nan), where=abs(dot)>64*np.finfo(float).eps*np.linalg.norm(N, axis=1))
        point = origin+t[:, None]*direction; v = point-triangles[:, 0]
        va = np.einsum('ij,ij->i', v, a); vb = np.einsum('ij,ij->i', v, b)
        u = (va*bb-vb*ab)/determinant; w = (vb*aa-va*ab)/determinant
        tol = m['barycentricTolerance']
        keep = (t>128*np.finfo(float).eps*np.linalg.norm(np.ptp(V, axis=0))) & (u>=-tol) & (w>=-tol) & (u+w<=1+tol)
        assert keep.any()
        nearest = t[keep].min(); ties = np.flatnonzero(keep & (abs(t-nearest)<=m['distanceToleranceM']))
        index = ties.min(); assert index+1 == m['originalFaceID'][j]
        assert N[index]@direction < 0
        expected = point[index]
        errors.append(float(np.linalg.norm(expected-m['bodyPosition'][j])))
        assert abs(t[index]-m['distanceM'][j]) < 1e-11
        np.testing.assert_allclose([1-u[index]-w[index], u[index], w[index]], m['barycentric'][j], rtol=0, atol=2e-10)
    assert max(errors)<1e-11
    assert not saved['miss']['hit'].any() and np.isnan(saved['miss']['distanceM']).all()
    assert not m['flowValuesAssignedToBody'] and not m['bodyMatchedByFlowSolution'] and not m['productionQualified']
    report = dict(status='passed-independent-original-face-plane-intersection-audit', rays=len(errors),
                  originalFaces=len(F), maximumPositionDifferenceM=max(errors),
                  sourceSHA256=sha(inlet_file), mappingSHA256=sha(mapping_file), processSHA256=sha(process_file),
                  flowValuesAssignedToBody=False, newProductionSolvedAreaM2=0)
    report_file.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf-8'));print(json.dumps(report))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('inlet',type=Path);p.add_argument('mapping',type=Path)
    p.add_argument('process',type=Path);p.add_argument('--report',required=True,type=Path)
    a=p.parse_args();audit(a.inlet,a.mapping,a.process,a.report)
