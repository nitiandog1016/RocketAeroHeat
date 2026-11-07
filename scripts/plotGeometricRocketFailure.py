"""只用既有原面诊断CSV绘图，不重复求解/审核；颜色明确属于拒绝的异常场。"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def main():
    root=Path(__file__).resolve().parents[1];base=root/'outputs/v196_validation';out=base/'audit01'
    surface=json.loads((base/'rocket01/surface.json').read_text());v=np.array(surface['vertices']);f=np.array(surface['facesOneBased'])-1
    wall=np.loadtxt(out/'rejected_wall_diagnostic.csv',delimiter=',',skiprows=1)
    residual=np.loadtxt(out/'normalized_residuals.csv',delimiter=',',skiprows=1)
    assert len(wall)==len(f)==19752
    fig=plt.figure(figsize=(12,4.6));ax=fig.add_subplot(121,projection='3d')
    norm=matplotlib.colors.Normalize(200,2200);colors=plt.cm.cividis(norm(wall[:,7]));colors[wall[:,8]==0]=[.85,.16,.2,1]
    ax.add_collection3d(Poly3DCollection(v[f],facecolors=colors,linewidth=0))
    lo,hi=v.min(0),v.max(0);ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]))
    ax.set_box_aspect(hi-lo);ax.view_init(20,-70);ax.set_axis_off()
    ax.set_title('Rejected exception field on original rocket\nRed: wall-neighbor state outside declared domain',fontsize=11,y=.90)
    ax.text2D(.04,.12,'Original geometry; nose +X →\nAll values rejected for engineering use',transform=ax.transAxes,fontsize=9)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='cividis'),ax=ax,shrink=.7,pad=0,label='Raw fluid T (K); not wall T')
    ax=fig.add_subplot(122)
    for i,label in enumerate(['mass','x momentum','y momentum','z momentum','energy']):ax.semilogy(residual[:,0],residual[:,i+1],label=label)
    ax.axhline(1e-8,color='gray',ls='--',label='Required 1e-8');ax.set(xlabel='Reported step',ylabel='Dimensionless residual',title='Finite time steps; failed development')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'rejected_development_final.png',dpi=180,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':main()
