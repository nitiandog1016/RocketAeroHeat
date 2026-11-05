"""原生轴对称壁面采样与完整体场配对，保留精确网格和输出截位界。"""
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from auditSU2SSTState import vtk_tokens
from su2AxisymmetricAnalysis import fields, wall_values, sample_stations, integrate_scalar


def surface_frame(path, points, markers, contract):
    """核验全部原壁节点及有向连接，再用精确坐标计算测站和旋转热功率。"""
    path = Path(path)
    xyz, values = vtk_tokens(path)
    ids = np.unique(markers['wall'])
    assert len(xyz) == len(ids) and np.isfinite(xyz).all()
    distance, match = cKDTree(points[ids]).query(xyz[:, :2])
    assert max(distance) <= 6e-8 and np.all(xyz[:, 2] == 0)
    assert len(np.unique(match)) == len(ids), '壁面节点配对不唯一'
    global_ids = ids[match]
    # 既有VTK解析器读取场量，这里额外核验原生线段的有向连接。
    words = path.read_text().split()
    offset = words.index('CELLS')
    count, size = map(int, words[offset+1:offset+3])
    assert count == len(markers['wall']) and size == 3*count
    cells = np.array(words[offset+3:offset+3+size], dtype=int).reshape(count, 3)
    assert np.all(cells[:, 0] == 2) and min(cells[:, 1:].ravel()) >= 0
    assert max(cells[:, 1:].ravel()) < len(ids)
    native_edges = global_ids[cells[:, 1:]]
    assert sorted(map(tuple, native_edges)) == sorted(map(tuple, markers['wall']))
    cell_offset = words.index('CELL_TYPES')
    assert int(words[cell_offset+1]) == count
    assert set(words[cell_offset+2:cell_offset+2+count]) == {'3'}
    f = {}
    for name in ('Pressure', 'Heat_Flux', 'Temperature', 'Y_Plus'):
        assert values[name].shape == (len(ids), 1) and np.isfinite(values[name]).all()
        f[name] = np.full(len(points), np.nan)
        f[name][global_ids] = values[name][:, 0]
    shear = values['Skin_Friction_Coefficient']
    assert shear.shape == (len(ids), 3) and np.isfinite(shear).all()
    for j, axis in enumerate(('x', 'y')):
        f['Skin_Friction_Coefficient_'+axis] = np.full(len(points), np.nan)
        f['Skin_Friction_Coefficient_'+axis][global_ids] = shear[:, j]
    assert min(f['Pressure'][ids]) > 0 and min(f['Temperature'][ids]) > 0
    wall = wall_values(points, markers, f)
    stations = sample_stations(wall, contract)
    return dict(iteration=int(path.stem.split('_')[-1]),
                p=stations['pressure']['predicted'], q=stations['heat']['predicted'],
                power=integrate_scalar(points, markers['wall'], f['Heat_Flux']),
                wall=wall, coordinateErrorM=float(max(distance)),
                sampling='native-surface-ASCII; exact original mesh integration')


def compare_volume(frame, path, points, markers):
    """同索引体场独立精度较高，逐壁点核对VTK截位及整周积分。"""
    assert frame['iteration'] == int(Path(path).stem.split('_')[-1])
    f = fields(path, True)
    wall = wall_values(points, markers, f)
    assert np.array_equal(frame['wall']['ids'], wall['ids'])
    errors = {key: float(max(abs(frame['wall'][key]-wall[key])/np.maximum(abs(wall[key]), 1.)))
              for key in ('p', 'q', 'fluidT', 'yplus', 'cf')}
    power = integrate_scalar(points, markers['wall'], f['Heat_Flux'])
    errors['power'] = abs(frame['power']-power)/max(abs(power), 1.)
    assert max(errors.values()) <= 6e-7, errors
    return errors
