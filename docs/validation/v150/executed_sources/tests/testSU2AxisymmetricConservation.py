"""用常量轴向流和解析线性温度检查旋转积分、控制体和WLS梯度。"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from su2AxisymmetricAnalysis import median_volumes, wls_gradient, integrate_vector, integrate_scalar


def verify():
    points = np.array([[0., .1], [1., .1], [1., .5], [0., .5]])
    cells = np.array([[0, 1, 2, 3]])
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])
    area, volume = median_volumes(points, cells)
    assert abs(sum(area)-.4) < 1e-14
    assert abs(sum(volume)-np.pi*(.5**2-.1**2)) < 1e-14
    values = 3*points[:, 0]+7*points[:, 1]+11
    gradient = wls_gradient(points, cells, values)[:, 0]
    assert np.max(abs(gradient-np.array([3., 7.]))) < 1e-13
    uniform = integrate_vector(points, edges, np.tile([19., 0.], (4, 1)))
    assert abs(uniform['net']) < 1e-13 and uniform['incoming'] > 0
    linear_heat = integrate_vector(points, edges, -gradient)
    assert abs(linear_heat['net']+2*np.pi*7*.4) < 1e-13
    surface = integrate_scalar(points, edges[:1], np.ones(4))
    assert abs(surface-2*np.pi*.1) < 1e-14
    return dict(passed=True, revolvedVolumeM3=float(sum(volume)),
                linearConductionNetW=linear_heat['net'], uniformNet=uniform['net'])


if __name__ == '__main__':
    import json
    print(json.dumps(verify()))
