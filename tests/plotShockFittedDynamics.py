"""从实际MAT和独立审核生成定量图；不重算流场，不把参考极限当实验分布。"""

from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'docs/validation/v102/figures'


def main():
    """模块一：解析体平均误差随实际三维网格变化，不拟合或宣称收敛阶。"""
    DEST.mkdir(exist_ok=True)
    data = json.loads((ROOT/'docs/validation/v102/dynamics_final_audit.json').read_text())
    rows = data['tubes']
    nx = [12,24,48]
    fig, ax = plt.subplots(figsize=(6.5,3.8),constrained_layout=True)
    for field, label in ((0,'Density / initial left density'),(4,'Energy / initial left energy')):
        ax.loglog(nx,[r['conservedCellAverageL1ByPhysicalScale'][field] for r in rows],'-o',label=label)
    ax.set_xticks(nx);ax.set_xticklabels(nx)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set(xlabel='Axial cells Nx',ylabel='Volume-weighted relative L1 error',
           title='Sod-type 3D calculation vs exact tetrahedron averages')
    ax.grid(True,which='both',alpha=.25);ax.legend(fontsize=9)
    fig.savefig(DEST/'sod_cell_average_error.png',dpi=170);plt.close(fig)

    # 模块二：不同网格/时间的球鼻原面壁压与实际残差；驻止极限只在轴线上画参考点。
    fig, axes = plt.subplots(1,2,figsize=(11,3.8),constrained_layout=True)
    for folder,label,color in (('cap02','84 faces, 9 us','#28718e'),
                               ('cap03','198 faces, 9 us','#be7e37'),
                               ('cap04','198 faces, 20 us','#567b68')):
        study = loadmat(ROOT/'outputs/v102_validation'/folder/'curved_cap_evidence.mat',
                        simplify_cells=True)['study']
        shell = study['cap']['shell']; mesh = shell['mesh']
        boundary = np.asarray(mesh['boundaryFaces'],int)-1
        walls = np.asarray(shell['wallBoundaryMask'],bool)
        faces = np.asarray(mesh['faces'],int)[boundary[walls]]-1
        center = np.asarray(study['nodes'])[faces].mean(axis=1)
        theta = np.degrees(np.arctan2(np.linalg.norm(center[:,1:],axis=1),center[:,0]))
        pressure = np.asarray(study['lastStep']['wallPressurePa'])[walls]
        axes[0].scatter(theta,pressure/1000,s=8,alpha=.7,label=label,color=color)
        history = np.asarray(study['history'])
        axes[1].semilogy(history[:,0]*1e6,history[:,2],label=label,color=color)
    axes[0].scatter([0],[39.08711980786282],marker='*',s=90,color='#333333',
                    label='Normal-shock stagnation limit (reference only)')
    axes[0].set(xlabel='Facet-centroid polar angle / deg',ylabel='Transient wall pressure / kPa',
                title='M5, alpha 0 deg, radius 0.05 m')
    axes[1].axhline(5e-4,ls='--',color='#777777',label='State-rate gate (other gates also required)')
    axes[1].set(xlabel='Physical time / us',ylabel='Dimensionless state-change rate',
                title='Not steady; shock compatibility not qualified')
    for ax in axes:
        ax.grid(True,alpha=.25);ax.legend(fontsize=7.5)
    fig.savefig(DEST/'cap_pressure_and_residual.png',dpi=170);plt.close(fig)
    print('Saved two quantitative figures from existing evidence.')


if __name__=='__main__':
    main()
