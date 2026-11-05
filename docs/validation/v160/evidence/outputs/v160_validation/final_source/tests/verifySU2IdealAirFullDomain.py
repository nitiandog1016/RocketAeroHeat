"""独立核对真实变比热整箭的原面CSV、保守热负荷和周向交线。"""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from auditSU2SSTState import vtk_tokens
from su2ThermalRuntime import read, sha


def main():
    """按原节点和三角形逐项复算，输出身份/未知温度/物理资格一起核对。"""
    base = ROOT/'outputs/v160_validation'
    view, case = base/'main02', base/'real01'
    s = loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'), squeeze_me=True, struct_as_record=False)['s']
    V, F = s.V, s.F.astype(int)-1
    av = np.cross(V[F[:, 1]]-V[F[:, 0]], V[F[:, 2]]-V[F[:, 0]])/2
    area = np.linalg.norm(av, axis=1)
    xyz, native = vtk_tokens(case/'surface_flow.vtk')
    distance, ids = cKDTree(V).query(xyz)
    assert max(distance) < 1e-6 and len(set(ids)) == len(V)
    with (view/'fields.csv').open(encoding='utf8') as stream:
        rows = list(csv.DictReader(stream))
    def column(name):
        return np.array([float(row[name]) for row in rows])
    assert len(rows) == 19752 and np.array_equal(column('originalFaceID'), np.arange(1, len(F)+1))
    errors = {}
    node_q = np.empty(len(V))
    node_q[ids] = native['Heat_Flux'][:, 0]
    for output, key in [('pressureRaw_Pa', 'Pressure'), ('heatFluxRaw_W_m2', 'Heat_Flux'),
                        ('fluidNodeTemperatureMean_K', 'Temperature')]:
        node = np.empty(len(V))
        node[ids] = native[key][:, 0]
        expected = node[F].mean(axis=1)
        error = float(max(abs(column(output)-expected)/np.maximum(1, abs(expected))))
        assert error < 1e-12
        errors[output] = error
    node_area = np.zeros(len(V))
    node_vector = np.zeros_like(V)
    for j in range(3):
        np.add.at(node_area, F[:, j], area/3)
        np.add.at(node_vector, F[:, j], av/3)
    native_area = np.linalg.norm(node_vector, axis=1)
    load = (node_q*native_area/node_area)[F].mean(axis=1)
    errors['nativeLoadDensity_W_m2'] = float(max(abs(column('nativeLoadDensity_W_m2')-load)/np.maximum(1, abs(load))))
    assert errors['nativeLoadDensity_W_m2'] < 1e-12
    assert max(abs(column('area_m2')/area-1)) < 1e-12
    assert np.all(column('givenWallTemperature_K') == 300) and np.all(np.isnan(column('solvedWallTemperature_K')))
    assert not np.any(column('productionQualified')) and not np.any(column('heatFluxIsBoundaryInput'))
    with (view/'section.csv').open(encoding='utf-8-sig') as stream:
        section = list(csv.reader(stream))
    assert len(section) == 93
    endpoint_error, heat_error = 0., 0.
    for row in section[1:]:
        x, face = float(row[0]), int(row[2])-1
        assert x == -1.5
        a = np.array([x, float(row[5]), float(row[6])])
        b = np.array([x, float(row[7]), float(row[8])])
        tri = V[F[face]]
        assert abs(np.linalg.norm(a-b)-float(row[9])) < 1e-12
        for point in (a, b):
            matrix = (tri[1:]-tri[0]).T
            uv = np.linalg.lstsq(matrix, point-tri[0], rcond=None)[0]
            error = np.linalg.norm(matrix@uv-(point-tri[0]))
            endpoint_error = max(endpoint_error, float(error))
            assert min(uv) > -1e-10 and sum(uv) < 1+1e-10 and error < 1e-10
        heat_error = max(heat_error, abs(float(row[10])-column('heatFluxRaw_W_m2')[face])/max(1, abs(float(row[10]))))
        assert 'SST-2003m' in row[12] and 'OTH' in row[12] and '未授予' in row[11]
    assert heat_error < 1e-12
    for name in ('fields.csv.json', 'section.csv.json'):
        metadata = read(view/name)
        audit = metadata['convergenceAudit']
        assert audit['checkpointIdentityConfirmed'] and audit['checkpointIteration'] == 149
        assert audit['completedUpdatesBeforeRun'] == 0 and audit['cumulativeCheckpointIteration'] == 149
        assert 'OTH' in metadata['algorithmRuntime']['method'] and not metadata['productionQualified']
    report = dict(passed=True, originalFaces=len(F), originalAreaM2=float(sum(area)), sectionSegments=92,
        rawFieldRelativeErrors=errors, sectionHeatRelativeError=heat_error, maximumEndpointTriangleDistanceM=endpoint_error,
        physicalSurfaceHeatPowerW=float(area@column('heatFluxRaw_W_m2')),
        nativeConservativeHeatPowerW=float(area@column('nativeLoadDensity_W_m2')),
        inputSHA256={str(p): sha(p) for p in [Path(__file__).resolve(), case/'surface_flow.vtk',
            view/'fields.csv', view/'section.csv', view/'fields.csv.json', view/'section.csv.json']},
        physicsQualified=False, newQualifiedArea=0)
    out = base/'export_check.json'
    assert not out.exists()
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    print(json.dumps({key: value for key, value in report.items() if key != 'inputSHA256'}, indent=2))


if __name__ == '__main__':
    main()
