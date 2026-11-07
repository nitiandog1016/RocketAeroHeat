"""独立验证Eilmer原SU2混合域准备；只准备小输入，不启动任何原生程序。

模块一：建立四类单元全face连通小域，解析体积17/8，独立检查参考Jacobian。
模块二：显式重现原RAH SU2棱柱次序，不从生产准备器借用拓扑或体积公式。
模块三：逐项核对原边界顺序/身份、源路径引用、无复制及错误输入拒绝。
"""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

# VTK/标准参考单元的外向面；由测试固定，不导入生产面表。
FACES = {
    10: ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
    12: ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
         (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
    13: ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
    14: ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
}
WEDGE_INVERSE = (1, 0, 2, 4, 3, 5)


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def det(a, b, c):
    return dot(a, cross(b, c))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixture():
    """四单元分别共享一个完整面；没有重叠，也不是四个不连通体。"""
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
              (2, .5, .5), (2, .5, -.5),
              (1.5, -1, 0), (2.5, -.5, -.5), (2.5, -.5, .5)]
    canonical = [(12, (0, 1, 2, 3, 4, 5, 6, 7)),
                 (14, (1, 2, 6, 5, 8)),
                 (10, (1, 2, 8, 9)),
                 (13, (1, 9, 8, 10, 11, 12))]
    source = [(kind, tuple(nodes[i] for i in WEDGE_INVERSE) if kind == 13 else nodes)
              for kind, nodes in canonical]
    owners = {}
    for cell_id, (kind, nodes) in enumerate(canonical):
        for face in FACES[kind]:
            outward = tuple(nodes[i] for i in face)
            owners.setdefault(tuple(sorted(outward)), []).append((cell_id, outward))
    assert sum(len(v) == 2 for v in owners.values()) == 3
    assert all(len(v) in (1, 2) for v in owners.values())
    outside = [v[0][1] for v in owners.values() if len(v) == 1]
    assert len(outside) == 14
    # 使用非自然顺序并交错两种标记，使sort/reverse/丢面都可被检出。
    wall_ids = [13, 0, 8, 4, 10, 2, 6]
    far_ids = [i for i in reversed(range(len(outside))) if i not in wall_ids]
    markers = [('wall', [outside[i] for i in wall_ids]),
               ('farfield', [outside[i] for i in far_ids])]
    return points, canonical, source, markers


def signed_surface_volume(points, kind, nodes):
    total = 0.0
    for face in FACES[kind]:
        ids = [nodes[i] for i in face]
        for j in range(1, len(ids)-1):
            total += det(points[ids[0]], points[ids[j]], points[ids[j+1]]) / 6.0
    return total


def reference_jacobians(points, kind, nodes):
    """各类型解析形函数导数；棱柱/四面体无生产median体积依赖。"""
    xyz = [points[i] for i in nodes]
    if kind == 10:
        return [det(sub(xyz[1], xyz[0]), sub(xyz[2], xyz[0]), sub(xyz[3], xyz[0]))]
    samples = []
    if kind == 12:
        signs = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                 (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
        for r in (-1, 0, 1):
            for s in (-1, 0, 1):
                for t in (-1, 0, 1):
                    derivatives = [(a*(1+b*s)*(1+c*t)/8, b*(1+a*r)*(1+c*t)/8,
                                    c*(1+a*r)*(1+b*s)/8) for a, b, c in signs]
                    samples.append(derivatives)
    elif kind == 13:
        for r, s in ((0, 0), (1, 0), (0, 1), (1/3, 1/3)):
            for t in (-1, 0, 1):
                bary = (1-r-s, r, s)
                dr, ds = (-1, 1, 0), (-1, 0, 1)
                samples.append([(dr[j]*(1-t)/2, ds[j]*(1-t)/2, -bary[j]/2) for j in range(3)] +
                               [(dr[j]*(1+t)/2, ds[j]*(1+t)/2, bary[j]/2) for j in range(3)])
    elif kind == 14:
        for r in (-1, 0, 1):
            for s in (-1, 0, 1):
                # 金字塔塌缩参考坐标在顶点t=1的零Jacobian是参数化退化，不采顶点。
                for t in (0, .25, .5, .75):
                    d = [(a*(1+b*s)*(1-t)/4, b*(1+a*r)*(1-t)/4,
                          -(1+a*r)*(1+b*s)/4) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
                    samples.append(d + [(0, 0, 1)])
    result = []
    for derivatives in samples:
        columns = [tuple(sum(xyz[j][axis]*derivatives[j][direction] for j in range(len(xyz)))
                         for axis in range(3)) for direction in range(3)]
        result.append(det(*columns))
    return result


def write_mesh(path, points, cells, markers, dimension=3):
    rows = ['NDIME= %d' % dimension, 'NELEM= %d' % len(cells)]
    for i, (kind, ids) in enumerate(cells):
        rows.append(' '.join(str(x) for x in (kind,) + tuple(ids) + (i,)))
    rows.append('NPOIN= %d' % len(points))
    for i, p in enumerate(points):
        rows.append(' '.join(format(float(x), '.17g') for x in p) + ' %d' % i)
    rows.append('NMARK= %d' % len(markers))
    for tag, faces in markers:
        rows += ['MARKER_TAG= '+tag, 'MARKER_ELEMS= %d' % len(faces)]
        for ids in faces:
            rows.append(' '.join(str(x) for x in ((5 if len(ids) == 3 else 9),) + tuple(ids)))
    Path(path).write_text('\n'.join(rows)+'\n', encoding='ascii')


def verify(output_dir):
    from prepareEilmerRocketCase import prepare_case
    started = time.monotonic()
    out = Path(output_dir).resolve()
    if out.exists():
        raise ValueError('测试输出目录已存在，不能覆盖准确失败或旧证据')
    out.mkdir(parents=True)
    points, canonical, source_cells, markers = fixture()
    expected = {12: 1.0, 14: 1/3, 10: 1/6, 13: 5/8}
    independent = {}
    for kind, ids in canonical:
        volume = signed_surface_volume(points, kind, ids)
        jac = reference_jacobians(points, kind, ids)
        assert abs(volume-expected[kind]) < 3e-15 and min(jac) > 0
        independent[str(kind)] = dict(volumeM3=volume, minimumSampledReferenceJacobian=min(jac), samples=len(jac))
    assert abs(sum(v['volumeM3'] for v in independent.values())-17/8) < 3e-15
    raw_wedge = source_cells[-1][1]
    raw_volume = signed_surface_volume(points, 13, raw_wedge)
    raw_jac = reference_jacobians(points, 13, raw_wedge)
    assert abs(raw_volume+5/8) < 3e-15 and max(raw_jac) < 0
    assert tuple(raw_wedge[i] for i in WEDGE_INVERSE) == canonical[-1][1]

    condition = dict(pressurePa=1197.0279007383137, temperatureK=226.50908361133003,
                     velocityMps=[-1400.0, 500.0, 0.0], wallTemperatureK=300.0,
                     tkeM2s2=2.0, omegaPerS=120.0, maxNewtonSteps=4,
                     stopRelativeResidual=1e-8)
    mesh_path = out/'mixed_original.su2'
    write_mesh(mesh_path, points, source_cells, markers)
    source_before = sha(mesh_path)
    positive_dir = out/'positive'
    result = prepare_case(mesh_path, positive_dir, condition, source_order='rah-su2-v1')
    assert result['status'] == 'prepared-only'
    assert result['nativeVerified'] is False and result['productionQualified'] is False
    mesh = result['mesh']
    assert mesh['nodes'] == 13 and mesh['cells'] == 4
    assert {str(k): v for k, v in mesh['cellTypeCounts'].items()} == {'10': 1, '12': 1, '13': 1, '14': 1}
    assert mesh['wallFaces'] == mesh['farfieldFaces'] == 7
    assert abs(mesh['discreteVolumeM3']-17/8) < 2e-14
    assert math.isfinite(mesh['minimumMedianSubVolumeM3']) and mesh['minimumMedianSubVolumeM3'] > 0
    assert mesh['sourceSHA256'] == source_before and Path(mesh['sourcePath']).resolve() == mesh_path
    assert mesh['sourceOrder'] == 'rah-su2-v1' and mesh['expectGmshOrderForWedges'] is True
    assert sha(mesh_path) == source_before
    assert not list(positive_dir.rglob('*.su2')), '准备不能复制原网格'
    for filename in ('gas.inp', 'job.lua', 'grid.lua', 'case_manifest.json', 'PREPARATION.md'):
        assert (positive_dir/filename).is_file()
    assert json.loads((positive_dir/'case_manifest.json').read_text('utf8')) == result
    lua = (positive_dir/'grid.lua').read_text('utf8')
    assert re.search(r'expect_gmsh_order_for_wedges\s*=\s*false', lua)
    assert str(mesh_path).replace('\\', '/') in lua.replace('\\\\', '/').replace('\\', '/')
    native_condition = dict(flowMode='nativeMach', pressurePa=condition['pressurePa'],
        temperatureK=condition['temperatureK'], wallTemperatureK=300., Mach=5., alphaDeg=20., betaDeg=0.,
        turbulenceIntensity=.01, eddyViscosityRatio=10., maxNewtonSteps=4, stopRelativeResidual=1e-8)
    native_result = prepare_case(mesh_path, out/'native_mach', native_condition)
    assert native_result['condition']['Mach']==5 and native_result['condition']['alphaDeg']==20

    rejected = []

    def reject(name, test_mesh=mesh_path, test_condition=None, source_order='rah-su2-v1'):
        destination = out/('reject_'+name)
        inputs = copy.deepcopy(condition if test_condition is None else test_condition)
        try:
            prepare_case(test_mesh, destination, inputs, source_order=source_order)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(name+' unexpectedly accepted')
        assert not destination.exists(), name+' should fail before creating target'

    bad = out/'unknown_marker.su2'
    write_mesh(bad, points, source_cells, [('unknown-wall', markers[0][1]), markers[1]])
    reject('unknown_marker', bad)
    bad = out/'wrong_dimension.su2'
    write_mesh(bad, points, source_cells, markers, dimension=2)
    reject('wrong_dimension', bad)
    bad_points = list(points); bad_points[0] = (float('nan'), 0, 0)
    bad = out/'nonfinite_node.su2'; write_mesh(bad, bad_points, source_cells, markers)
    reject('nonfinite_node', bad)
    bad = out/'already_converted_wedge.su2'
    write_mesh(bad, points, canonical, markers)
    reject('already_converted_wedge', bad)
    bad_cells = list(source_cells); bad_cells[2] = (10, (1, 2, 8, 8))
    bad = out/'degenerate_tetra.su2'; write_mesh(bad, points, bad_cells, markers)
    reject('degenerate_tetra', bad)
    reject('unknown_source_order', source_order='unknown-order')
    for name, key, value in [('nonfinite_condition', 'wallTemperatureK', float('nan')),
                             ('velocity_dimension', 'velocityMps', [-1400, 500]),
                             ('zero_velocity', 'velocityMps', [0, 0, 0]),
                             ('noninteger_steps', 'maxNewtonSteps', 1.5)]:
        changed = copy.deepcopy(condition); changed[key] = value
        reject(name, test_condition=changed)
    changed = copy.deepcopy(condition); del changed['pressurePa']
    reject('missing_condition', test_condition=changed)
    snapshot = {p.name: sha(p) for p in positive_dir.iterdir() if p.is_file()}
    try:
        prepare_case(mesh_path, positive_dir, condition)
    except ValueError:
        rejected.append('existing_destination')
    else:
        raise AssertionError('existing destination overwritten')
    assert snapshot == {p.name: sha(p) for p in positive_dir.iterdir() if p.is_file()}
    assert sha(mesh_path) == source_before

    report = dict(passed=True, originalNodes=13, originalCells=4, connectedCellInterfaces=3,
                  originalBoundaryFaces=14, expectedVolumeM3=17/8, independentlyAuditedTypes=independent,
                  rawRAHSU2WedgeSignedVolumeM3=raw_volume, rawRAHSU2WedgeMaximumJacobian=max(raw_jac),
                  auditReferenceWedgeInverse=list(WEDGE_INVERSE), gmshOrderSwitch=False,
                  noMeshCopy=True, sourceUnchanged=True, existingDestinationUnchanged=True,
                  rejected=rejected, nativeRuns=0, physicsQualified=False,
                  elapsedSeconds=time.monotonic()-started,
                  inputSHA256={str(mesh_path): source_before, str(Path(__file__).resolve()): sha(__file__),
                               str(ROOT/'scripts/prepareEilmerRocketCase.py'): sha(ROOT/'scripts/prepareEilmerRocketCase.py')})
    (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', '--out', dest='output', required=True)
    print(json.dumps(verify(parser.parse_args().output), allow_nan=False))
