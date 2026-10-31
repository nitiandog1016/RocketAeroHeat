"""用真实保存的候选结果绘制方法与法向网格对照，不拟合新的物理解。"""

from pathlib import Path
import json

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from auditShockFittedInterface import current_wall_pressure

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'docs/validation/v103/figures'
ROWS = [
    ('outputs/v102_validation/cap03/curved_cap_evidence.mat', 'Face L2 / forced RH', '#536778'),
    ('outputs/v103_validation/cap01/curved_cap_evidence.mat', 'Nodal J+ / forced RH', '#a86a3d'),
    ('outputs/v103_validation/cap02/curved_cap_evidence.mat', 'Face L2 / actual-face Riemann', '#28718e'),
    ('outputs/v103_validation/cap03/curved_cap_evidence.mat', 'Nodal J+ / actual-face Riemann', '#64866a'),
]


def main():
    """模块一：相同时间/壁面网格四方法；保持原单元压力差，不用重构误差替代。"""
    DEST.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    for source, label, color in ROWS:
        s = loadmat(ROOT / source, simplify_cells=True)['study']
        h = s['history']
        axes[0].plot(h[:, 0] * 1e6, h[:, 15] * 100, label=label, color=color)
        shell = s['cap']['shell']; mesh = shell['mesh']
        boundary = np.asarray(mesh['boundaryFaces'], int) - 1
        walls = np.asarray(shell['wallBoundaryMask'], bool)
        xyz = s['nodes'][np.asarray(mesh['faces'], int)[boundary[walls]] - 1].mean(axis=1)
        theta = np.degrees(np.arctan2(np.linalg.norm(xyz[:, 1:], axis=1), xyz[:, 0]))
        axes[1].scatter(theta, current_wall_pressure(s)[walls] / 1000, s=7, color=color, alpha=.65, label=label)
    axes[0].set(xlabel='Step-end time / us', ylabel='Max original-cell pressure mismatch / %',
                title='198 wall facets, four radial layers')
    axes[1].set(xlabel='Wall facet-centroid polar angle / deg', ylabel='Current-frame wall pressure / kPa',
                title='Same saved time: 9 us')
    # 模块二：只改变法向层数；结果未下降，不冒称网格收敛。
    values, speeds = [], []
    for name in ('cap02', 'cap04', 'cap05'):
        item = json.loads((ROOT / 'docs/validation/v103' / (name + '_audit.json')).read_text(encoding='utf-8'))
        values.append(item['finalShockPressureMismatch'] * 100)
        speeds.append(item['sharedNodeMotion']['maximumUnrepresentedFaceSpeedMS'])
    axes[2].plot([4, 6, 8], values, '-o', color='#28718e')
    axes[2].set(xticks=[4, 6, 8], xlabel='Radial layers (same wall facets)',
                ylabel='Max pressure mismatch / %', title='Radial refinement does not close mismatch')
    for ax in axes:
        ax.grid(True, alpha=.22)
    axes[0].legend(fontsize=7)
    fig.suptitle('M5, alpha 0 deg, R = 0.05 m: unqualified transient candidates', fontsize=12)
    fig.savefig(DEST / 'interface_method_and_radial_comparison.png', dpi=170)
    plt.close(fig)
    print('Saved actual method and radial-sensitivity comparison; no experimental qualification.')


if __name__ == '__main__':
    main()
