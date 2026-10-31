"""输出可单独查阅的残差定位图；颜色表示原生积分残差，不能解释为热流误差。"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
from diagnoseSU2Checkpoint import read_volume

ROOT=Path(__file__).resolve().parents[1]


def main():
    folder=ROOT/'outputs/v110_validation';r=json.loads((folder/'spatial_report_final.json').read_text())
    xyz,f=read_volume(folder/'probe02/volume_flow.vtk')
    surface=loadmat(str(ROOT/'outputs/v108_validation/normal01/normal_volume.mat'),squeeze_me=True,struct_as_record=False)['surface']
    magnitude=abs(f['Residual_Energy']);indices=np.argsort(magnitude)[-max(1,len(xyz)//100):]
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.5),gridspec_kw={'width_ratios':[1.6,1]})
    ax=axes[0];body=surface.V;profile=np.unique(np.round(np.column_stack((body[:,0],np.linalg.norm(body[:,1:],axis=1))),10),axis=0)
    ax.plot(profile[:,0],profile[:,1],color='#343e46',lw=1.2,label='Original body meridian')
    artist=ax.scatter(xyz[indices,0],np.linalg.norm(xyz[indices,1:],axis=1),c=np.log10(magnitude[indices]),cmap='viridis',s=17,alpha=.85)
    ax.axvline(0,color='#777',lw=.8,ls='--');ax.set(xlabel='X (m), nose toward +X',ylabel='Radius from X axis (m)',title='Top 1% of native energy residual magnitude')
    ax.legend(loc='upper right',fontsize=8);fig.colorbar(artist,ax=ax,label='log10 |integrated energy residual|')
    layers=r['layerResidualSquareFraction'];axes[1].bar([x['layer'] for x in layers],[100*x['residualSquareFraction']['rmsRhoE'] for x in layers],color='#376b84')
    axes[1].axvline(24.5,color='#777',lw=.8,ls='--');axes[1].set(xlabel='Layer index: wall=0, far boundary=40',ylabel='Contribution to sum of squared residual (%)',title='Native residual norm by layer')
    for ax in axes:ax.grid(alpha=.18);ax.set_axisbelow(True)
    fig.suptitle('M5 / alpha 20 deg / Tw 300 K: one diagnostic sweep from saved iteration 3750',fontsize=11)
    fig.text(.5,.015,'Integrated finite-volume imbalance, not volume-normalized error or heat-flux accuracy. Unconverged research state.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.045,1,.94]);path=ROOT/'docs/validation/v110/residual_location.png';assert not path.exists();fig.savefig(str(path),dpi=160);plt.close(fig)
    print(str(path))


if __name__=='__main__':main()
