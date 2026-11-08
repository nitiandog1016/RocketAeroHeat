"""模块职责：只读最终分区重绘清晰科研图，不重新划域、计算资源或授予热资格。"""
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from eilmerNKSAResourceBudget import ROOT, sha, save


def main():
    base = ROOT / 'outputs/v216_validation'
    out = base / 'plot01'
    out.mkdir(exist_ok=False)
    paths = [ROOT / 'outputs/v144_validation/actual/hybrid_volume.mat',
             base / 'connected01/partition.npz', base / 'connected01/report.json',
             base / 'final01/report.json']
    before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    report = json.loads(paths[2].read_text())
    assert sha(paths[1]) == report['partitionSHA256']
    assert json.loads(paths[3].read_text())['allRegionsConnected']
    data = loadmat(str(paths[0]), variable_names=['sourceVertices', 'sourceFaces'])
    vertices = data['sourceVertices']
    faces = data['sourceFaces'].astype(np.int64) - 1
    with np.load(str(paths[1])) as stored:
        owner = stored['wallOwnerRegion']

    # 原面颜色表示唯一责任，不表示压力/热流；隐藏原图挤在一起的三维刻度。
    colors = ['#446e91', '#d28b49', '#689679', '#9b7aa0', '#647f98', '#ae9b60', '#a87373']
    labels = ['Nose', 'Contraction', 'Short aft', 'Aft 1', 'Aft 2', 'Aft 3', 'Base / wake']
    fig = plt.figure(figsize=(12, 5.3), facecolor='white')
    ax = fig.add_axes([.015, .20, .49, .66], projection='3d')
    ax.set_proj_type('ortho')
    ax.add_collection3d(Poly3DCollection(vertices[faces], facecolors=np.array(colors)[owner],
                                       edgecolors='none', linewidths=0, rasterized=True))
    ax.set_xlim(vertices[:, 0].min(), vertices[:, 0].max())
    ax.set_ylim(-.17, .17)
    ax.set_zlim(-.17, .17)
    ax.set_box_aspect((np.ptp(vertices[:, 0]), .34, .34))
    ax.view_init(elev=16, azim=-68)
    ax.set_axis_off()
    fig.text(.25, .93, 'Original rocket: unique wall responsibility', ha='center', fontsize=12)
    fig.text(.25, .88, '19,752 faces | 1.9646292566 m²', ha='center', fontsize=11)
    fig.text(.25, .26, 'Original geometry and all 32 wall layers retained', ha='center', fontsize=10)
    fig.text(.25, .20, 'X = -2.8283 to +0.0500 m | nose at +X (right)', ha='center', fontsize=10)
    # 条形图复用已审核整数资源账，不把必要载荷当运行峰值。
    bx = fig.add_axes([.59, .17, .35, .65])
    values = [item['payloadLowerBoundGiB'] for item in report['regions']]
    bars = bx.barh(labels, values, color=colors, height=.65)
    bx.invert_yaxis()
    bx.set_xlim(0, 8.6)
    bx.axvline(8, color='#555555', linestyle='--', linewidth=1)
    bx.text(8.08, 6.5, 'WSL 8 GiB', rotation=90, ha='left', va='bottom', fontsize=9)
    for bar, value in zip(bars, values):
        bx.text(value + .07, bar.get_y() + bar.get_height()/2, '%.3f' % value, va='center', fontsize=10)
    bx.set_xlabel('Necessary solver payload (GiB)')
    bx.set_title('Full original domain included\nRun peak not established', fontsize=11, pad=12)
    bx.spines['top'].set_visible(False)
    bx.spines['right'].set_visible(False)
    bx.grid(axis='x', alpha=.15)
    bx.set_axisbelow(True)
    fig.text(.5, .035, 'Partition design only. No new solved pressure or heat flux.',
             ha='center', fontsize=10, color='#444444')
    fig.savefig(str(out/'regions.png'), dpi=150)
    plt.close(fig)
    assert all(sha(ROOT/key) == value for key, value in before.items())
    save(out/'report.json', dict(passed=True, sourceSHA256=before, originalFigurePreserved=True,
                                figureSHA256=sha(out/'regions.png'), numericalReevaluation=False))


if __name__ == '__main__':
    raise ValueError('v216离线设计已经结束；请按NEXT_DECISION登记新任务，当前无构建/查询/网格/CFD准入。')
