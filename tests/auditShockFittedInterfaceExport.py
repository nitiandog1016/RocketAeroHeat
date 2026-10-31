"""从原MAT独立核对当前帧壁压、诊断时刻、真实截线和逐面导出。"""

from pathlib import Path
import argparse
import csv
import io
import json

import numpy as np
from scipy.io import loadmat

from auditShockFittedALEKernel import geometry
from auditShockFittedCap import euler_ale
from auditShockFittedDynamics import primitive
from auditShockFittedDynamicsExport import FIELDS, numeric, intersect, contours
from auditShockFittedInterface import current_wall_pressure
from auditThermalPressureSurfaceFlow import sha

ROOT = Path(__file__).resolve().parents[1]


def audit(folder, sections, target):
    """模块一：所有主体场量取本帧，旧步诊断另核输入态/平均几何及五个通量分量。"""
    assert not target.exists()
    data = json.loads((folder / 'ale-data.json').read_text(encoding='utf-8'))
    assert data['kind'] == 'shock-fitted-interface-study-v103' and not data['productionQualified']
    assert data['newProductionSolvedAreaM2'] == 0 and not data['experimentalValidation']
    assert (folder / 'ale-data.js').read_bytes() == b'window.ALE_STUDY=' + (folder / 'ale-data.json').read_bytes() + b';\n'
    original = data['originalGeometrySource']
    assert sha(Path(original['evidenceFile'])) == original['evidenceSHA256']
    legacy = loadmat(original['evidenceFile'], simplify_cells=True)
    np.testing.assert_array_equal(data['originalGeometry']['vertices'], legacy['geometry']['V'])
    np.testing.assert_array_equal(data['originalGeometry']['faces'], legacy['geometry']['F'] - 1)
    checked = 0
    maximum_wall_error, maximum_defect_error = 0., 0.
    for case in data['cases']:
        source = case['source']
        assert sha(Path(source['evidenceFile'])) == source['evidenceSHA256']
        assert sha(Path(source['auditFile'])) == source['auditSHA256']
        saved = loadmat(source['evidenceFile'], simplify_cells=True)['study']
        shell, gas = saved['cap']['shell'], saved['gas']
        mesh = shell['mesh']
        boundary = np.asarray(mesh['boundaryFaces'], int) - 1
        owner = np.asarray(mesh['owner'], int)[boundary] - 1
        np.testing.assert_array_equal(case['faces'], mesh['faces'][boundary] - 1)
        np.testing.assert_array_equal(case['boundaryKind'], saved['boundary']['kind'])
        np.testing.assert_array_equal(case['sourceFaceIDs'], shell['cellOriginalFaceIDs'][owner])
        assert case['motionMethod'] == saved.get('motionMethod', 'face-riemann-L2')
        assert case['shockFluxMethod'] == saved.get('shockFluxMethod', 'RH')
        assert not saved['productionQualified'] and not case['productionQualified']
        for j, frame in enumerate(case['frames']):
            Q = saved['Q'] if j else saved['initialQ']
            X = saved['nodes'] if j else saved['initialNodes']
            expected_time = saved['timeS'] if j else 0
            assert frame['timeS'] == frame['wallPressureTimeS'] == expected_time
            np.testing.assert_array_equal(frame['position'], X)
            state = primitive(Q, gas)
            expected = {key: np.full(len(owner), np.nan) for key in FIELDS}
            expected.update(pressurePa=state[owner, 0] * gas['R'] * state[owner, 4],
                            temperatureK=state[owner, 4], speedMS=np.linalg.norm(state[owner, 1:4], axis=1),
                            wallPressurePa=current_wall_pressure(saved, Q=Q, nodes=X))
            for key in FIELDS:
                actual = numeric(frame[key])
                np.testing.assert_allclose(actual, expected[key], rtol=3e-12, atol=1e-12, equal_nan=True)
                checked += len(actual)
            wall_error = np.nanmax(abs(numeric(frame['wallPressurePa']) / expected['wallPressurePa'] - 1))
            maximum_wall_error = max(maximum_wall_error, float(wall_error))
            d = frame['diagnostics']
            assert bool(d['available']) == bool(j)
            if not j:
                continue
            dt = saved['history'][-1, 1]
            assert d['inputTimeS'] == expected_time - dt
            np.testing.assert_array_equal(d['geometryTimeSpanS'], [expected_time - dt, expected_time])
            for key, reference in (
                ('cellPressureMismatch', saved['lastStep']['shockInteriorPressureRelativeMismatch']),
                ('stepInputWallPressurePa', saved['lastStep']['wallPressurePa'])):
                np.testing.assert_allclose(numeric(d[key]), reference, rtol=4e-15, atol=0, equal_nan=True)
            outer = np.asarray(shell['outerBoundaryMask'], bool)
            for key, reference in (('requestedFaceSpeedMS', 'requestedFaceNormalSpeedMS'),
                                   ('projectedFaceSpeedMS', 'faceNormalSpeedMS'),
                                   ('faceSpeedResidualMS', 'projectionResidualMS')):
                expected_array = np.full(len(boundary), np.nan)
                expected_array[outer] = saved['lastMotion'][reference]
                np.testing.assert_allclose(numeric(d[key]), expected_array, rtol=4e-15, atol=0, equal_nan=True)
            geo = geometry(mesh, saved['lastNodes0'], saved['nodes'], dt)
            area = np.linalg.norm(geo['Abar'], axis=1)
            normals = geo['Abar'] / area[:, None]
            D = geo['swept'] / (dt * area)
            average = np.full(len(boundary), np.nan)
            average[outer] = D[boundary[outer]]
            np.testing.assert_allclose(numeric(d['stepAverageGridNormalSpeedMS']), average,
                                       rtol=2e-10, atol=2e-8, equal_nan=True)
            reconstructed = np.full(len(boundary), np.nan)
            if case['motionMethod'] == 'nodal-characteristic':
                count = len(shell['originalVertexIDs'])
                outer_faces = np.asarray(mesh['faces'], int)[boundary[outer]] - 1 - int(shell['layers']) * count
                nodal = saved['lastMotion']['nodalTrace']['primitiveRP']
                face_pressure = nodal[outer_faces, 4].mean(axis=1)
                reconstructed[outer] = face_pressure / saved['lastStep']['shockPressurePa'][outer] - 1
            np.testing.assert_allclose(numeric(d['reconstructedPressureMismatch']), reconstructed,
                                       rtol=2e-11, atol=2e-12, equal_nan=True)
            up = saved['upstream']
            U = np.asarray(up['velocityVector'])
            Qup = np.r_[up['density'], up['density'] * U,
                        up['pressure'] / (gas['gamma'] - 1) + up['density'] * (U @ U) / 2]
            speed = max(np.linalg.norm(U), np.sqrt(gas['gamma'] * gas['R'] * up['temperature']))
            scale = up['density'] * np.r_[speed, np.full(3, speed ** 2), speed ** 3]
            defect = np.full((len(boundary), 5), np.nan)
            for row in np.flatnonzero(outer):
                face = boundary[row]
                defect[row] = (saved['lastStep']['faceALEFlux'][face] -
                               euler_ale(Qup, normals[face], D[face], gas['gamma'])) / scale
            actual = np.asarray([[np.nan if v is None else v for v in r] for r in d['RHFluxDefect5']])
            np.testing.assert_allclose(actual, defect, rtol=2e-10, atol=2e-12, equal_nan=True)
            maximum_defect_error = max(maximum_defect_error, float(np.nanmax(abs(actual - defect))))
            checked += int(np.count_nonzero(outer)) * 5

    # 模块二：以独立重心约束/SVD和KD树核对全部周向导出，不复用JS求交。
    manifest = json.loads((sections / 'node_audit.json').read_text(encoding='utf-8'))
    rows_checked, max_endpoint = 0, 0.
    for item in manifest['sectionFiles']:
        file = sections / item['name']
        assert sha(file) == item['sha256']
        case = data['cases'][item['caseIndex']]
        frame = case['frames'][item['frameIndex']]
        ids = [i for i, kind in enumerate(case['boundaryKind'])
               if item['view'] == 'all' or (kind == 2 if item['view'] == 'wall' else kind in (3, 5))]
        faces = np.asarray(case['faces'], int)[ids]
        vertices = np.asarray(frame['position'])
        used = vertices[np.unique(faces)]
        size, mag = np.ptp(used, axis=0).max(), np.abs(used).max()
        eps = max(size * 1e-10, mag * np.finfo(float).eps * 128)
        tolerance = max(size * 1e-11, mag * np.finfo(float).eps * 512)
        x = item['x']
        side = -1 if x >= used[:, 0].max() - eps else 1
        segments = {}
        for k, face in enumerate(faces):
            segment = intersect(vertices[face], x, eps, side)
            if segment is not None:
                segments[k] = segment
        topology = contours(segments, x, eps)
        rows = list(csv.reader(io.StringIO(file.read_text(encoding='utf-8-sig'))))
        assert len(rows[0]) == 18 and len(rows) - 1 == item['rows'] == len(segments)
        observed = set()
        groups = {}
        for row in rows[1:]:
            assert len(row) == 18 and float(row[0]) == x
            local = int(row[2]) - 1
            assert local in segments and local not in observed
            observed.add(local)
            actual = np.array([[x, float(row[5]), float(row[6])], [x, float(row[7]), float(row[8])]])
            segment = segments[local]
            error = min(np.max(np.linalg.norm(actual - segment, axis=1)),
                        np.max(np.linalg.norm(actual - segment[::-1], axis=1)))
            assert error <= tolerance
            max_endpoint = max(max_endpoint, float(error))
            value = frame[item['field']][ids[local]]
            assert row[10] == '' if value is None else float(row[10]) == value
            assert row[12] == case['methods'][item['field']] and row[13] == case['scope']
            assert case['source']['evidenceSHA256'] in row[14]
            assert 't=' + str(frame['timeS']) in row[14] or float(row[14].split('t=')[1].split(' s')[0]) == frame['timeS']
            contour = topology[local]
            groups.setdefault(contour['group'], set()).add(int(row[1]))
            center = contour['center']
            a, b = actual[:, 1:] - center[1:]
            span = b - a
            fraction = np.clip(-a @ span / (span @ span), 0, 1)
            valid = (contour['closed'] and not contour['branched'] and
                     np.linalg.norm(a + fraction * span) > eps and
                     abs(np.cross(a, b)) > eps * (np.linalg.norm(a) + np.linalg.norm(b)))
            if valid:
                assert row[17] == '有效周向角'
                angles = np.degrees(np.arctan2(actual[:, 2] - center[2], actual[:, 1] - center[1])) % 360
                assert np.max(abs((np.array([float(row[3]), float(row[4])]) - angles + 180) % 360 - 180)) < 2e-9
            else:
                assert row[3] == row[4] == ''
            rows_checked += 1
        assert observed == set(segments) and all(len(v) == 1 for v in groups.values())

    # 模块三：绑定真正部署到查看器的资源，不把未提交的并发界面改动混入本次验收。
    assets = {}
    for asset in ('ale-study.js', 'ale-study.css', 'ale-section-core.js', 'geometry/geometry.js',
                  'geometry/geometry.css', 'geometry/ale-study-geometry.js'):
        assert (folder / asset).read_bytes() == (ROOT / 'web' / Path(asset).name).read_bytes()
        assets[asset] = sha(folder / asset)
    assert (folder / 'index.html').read_bytes() == (ROOT / 'web/ale-study.html').read_bytes()
    result = dict(status='passed-independent-current-frame-fields-wall-times-five-fluxes-and-SVD-sections',
                  dataSHA256=sha(folder / 'ale-data.json'), nodeTestSHA256=sha(ROOT / 'tests/testShockFittedInterfaceStudy.js'),
                  nodeAuditSHA256=sha(sections / 'node_audit.json'), cases=len(data['cases']),
                  frames=2 * len(data['cases']), checkedValues=checked,
                  independentWallRelativeDifference=maximum_wall_error,
                  independentRHDefectDifference=maximum_defect_error,
                  sectionFiles=len(manifest['sectionFiles']), sectionRows=rows_checked,
                  maximumEndpointErrorM=max_endpoint, originalFaces=len(data['originalGeometry']['faces']),
                  deployedAssets=assets, productionQualified=False, newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(result, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for key in ('folder', 'sections', 'target'):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    audit(args.folder, args.sections, args.target)
