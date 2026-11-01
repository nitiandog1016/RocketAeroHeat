"""只读已核验的预览，修正科研图布局；不重算质量或启动GUI/网格器。"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    base=ROOT/'outputs/v134_validation';out=base/'nose_quality_layout02.png';assert not out.exists()
    check=json.loads((base/'quality_check.json').read_text());assert check['passed']
    source=base/'workbench02/quality_preview.mat'
    assert sha(source)==check['sourceSHA256'][str(source)]
    d=loadmat(str(source),squeeze_me=True,struct_as_record=False);p=d['payload'];info=d['info']
    V=p.vertices;F=p.faces.astype(int);n=int(info.visibleCellFaces)
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    fig=plt.figure(figsize=(9,6.6));ax=fig.add_subplot(111,projection='3d')
    fig.subplots_adjust(left=.04,right=.79,bottom=.14,top=.85)
    ax.add_collection3d(Poly3DCollection(V[F[n:]],facecolors='#d9dfe2',edgecolors='none',alpha=.45))
    coll=Poly3DCollection(V[F[:n]],edgecolors='#273e50',linewidths=.4,cmap='viridis_r')
    coll.set_array(info.displayValues[:n]);coll.set_clim(0,.1);ax.add_collection3d(coll)
    lo=V.min(axis=0);hi=V.max(axis=0)
    ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_zlim(lo[2],hi[2]);ax.set_box_aspect(hi-lo)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(3));axis.set_tick_params(labelsize=8)
    ax.set_xlabel('X / m',labelpad=9);ax.set_ylabel('Y / m',labelpad=9);ax.set_zlabel('Z / m',labelpad=9);ax.view_init(22,-54)
    fig.suptitle('鼻部外层 · 17个真实低质量四面体\n灰色为原壁面；仅网格诊断，未增加气动热资格',fontsize=12,y=.96)
    cb=fig.colorbar(coll,cax=fig.add_axes([.87,.25,.024,.48]));cb.set_label('SICN / 1（越小越差）')
    fig.savefig(str(out),dpi=170,bbox_inches='tight',pad_inches=.22);plt.close(fig)
    (base/'plot_layout02.json').write_text(json.dumps(dict(sourceSHA256=sha(source),scriptSHA256=sha(Path(__file__).resolve()),
        imageSHA256=sha(out),numericalRecalculation=False,replacesClippedInitialPlot=True,nativeScreenshot=False),indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
