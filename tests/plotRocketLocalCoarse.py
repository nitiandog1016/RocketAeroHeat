# -*- coding: utf-8 -*-
"""绘制真实原面/控制体定位、局部像投影与整列存储预算；图中数值不是热流精度。"""
from pathlib import Path
import hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]; case = sys.argv[1]
assert case.replace('_', '').isalnum()
folder = ROOT / 'outputs/v83_validation' / case
revision = int(sys.argv[2]) if len(sys.argv) > 2 else 1; assert revision in (1, 2)
suffix = '' if revision == 1 else str(revision)
source = folder / 'local_coarse_evidence.mat'; target = folder / ('local_coarse_selection' + suffix + '.png')
assert not target.exists()
s = loadmat(source, simplify_cells=True); c = s['coarse']; r = s['report']
assert r['kind'] == 'real'
cut = c['space']['geometry']; geom = cut['inputGeometry']
V = geom['V']; F = geom['F'].astype(int) - 1; nf = len(F)
cells = np.atleast_1d(s['cellIDs']).astype(int) - 1
energy = np.sum(s['remaining'] ** 2, axis=1); total = energy.sum()
face_energy = energy.reshape((-1, nf)).sum(axis=0) / total
centres = cut['cells']['centroid']; layers = cells // nf + 1
history = np.atleast_2d(c['selectionHistory'])
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
fig = plt.figure(figsize=(12, 8), layout='constrained')
ax = fig.add_subplot(221, projection='3d')
normalization = LogNorm(vmin=1e-7, vmax=max(face_energy))
surface = Poly3DCollection(V[F], linewidths=0, rasterized=True)
surface.set_facecolor(plt.cm.viridis(normalization(np.maximum(face_energy, 1e-7))))
ax.add_collection3d(surface); ax.auto_scale_xyz(V[:, 0], V[:, 1], V[:, 2])
ax.set_box_aspect(np.maximum(np.ptp(V, axis=0), 1e-4)); ax.view_init(elev=20, azim=-70)
ax.set(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)', title='Original rocket: remainder squared-norm share per facet')
if revision >= 2:
    ax.set_axis_off()
    ax.text2D(.05, .14, f'Original proportions; body length {np.ptp(V[:, 0]):.3f} m; nose +X', transform=ax.transAxes)
bar = fig.colorbar(plt.cm.ScalarMappable(norm=normalization, cmap='viridis'), ax=ax, shrink=.6, pad=.02)
bar.set_label('Fraction; summed over the original normal CV column' if revision == 1 else 'Facet remainder fraction; not heat-flux error')
ax = fig.add_subplot(222)
colors = ax.scatter(centres[cells, 0], centres[cells, 1], c=layers, cmap='plasma', s=20 + 200 * np.sqrt(energy[cells] / total), alpha=.85)
ax.set(xlabel='Original CV centroid X (m)', ylabel='Original CV centroid Y (m)', title=f'{len(cells)} original seed CVs; size shows remainder contribution')
fig.colorbar(colors, ax=ax, label='Original normal layer')
ax.grid(alpha=.2)
ax = fig.add_subplot(223)
rhs_norm = np.linalg.norm(s['rhs']); baseline = r['oldMixedProjectedRelativeDefect']
ax.axhline(r['originalRemainingRelativeDefect'], color='.5', linestyle=':', label='Saved original-NS direction remainder')
ax.axhline(baseline, color='#61767c', linestyle='--', label='Existing 473-column projection')
ax.plot(history[:, 0], history[:, -1] * np.linalg.norm(s['remaining']) / rhs_norm, color='#c66b20', label='New local enriched projection')
ax.set(xlabel='Selected local conservative columns', ylabel='Predicted remainder / original RHS norm', title='Current frozen linear model; no flow update')
ax.grid(alpha=.2); ax.legend(fontsize=8)
ax = fig.add_subplot(224)
ax.plot(history[:, 0], history[:, 5] / 2**20, color='#326782', label='Stored final sparse image')
ax.axhline(c['maximumImageBytes'] / 2**20, color='#9c3b36', linestyle='--', label='Original 256 MiB limit')
ax.set(xlabel='Selected local conservative columns', ylabel='Sparse image storage (MiB)', title=f'{len(np.atleast_1d(c["removedMixedColumns"]))} whole primitive columns released; all CVs retained')
ax.grid(alpha=.2); ax.legend(fontsize=8)
fig.suptitle('M5 / 30 km / AoA 20 deg / supplied Tw 300 K — local linear-solver research', fontsize=13)
fig.savefig(target, dpi=160); plt.close(fig)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

out = dict(status='generated-scientific-local-CV-location-projection-storage-figure', sourceSHA256=sha(source),
    figureSHA256=sha(target), scriptSHA256=sha(__file__), revision=revision, physicalSolution=False, newProductionSolvedAreaM2=0)
(folder / ('local_coarse_figure' + suffix + '.json')).write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps(out))
