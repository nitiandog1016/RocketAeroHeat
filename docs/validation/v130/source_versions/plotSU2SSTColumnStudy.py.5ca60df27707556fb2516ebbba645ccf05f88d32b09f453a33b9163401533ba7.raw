"""绘制固定对照历史和原节点恢复数量；不是物理资格或新增面积证明。"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    source=ROOT/'outputs/v130_validation/column_comparison.json';r=read(source)
    out=ROOT/'docs/validation/v130/column_comparison.png';assert not out.exists()
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for name,ls in [('control','-'),('column','--')]:
        h=r['variants'][name]['history']
        for key,color in [('rms[RhoE]','#315c84'),('rms[w]','#a34632')]:
            axes[0,0].plot(h[key],ls=ls,color=color,label=name+' '+key)
        axes[0,1].plot([v/1000 for v in h['HF']],ls=ls,label=name)
        s=r['variants'][name]['snapshots']
        axes[1,0].plot([50,99],[s[str(i)]['recoveredOriginalNodesAbove100Floor'] for i in [50,99]],marker='o',ls=ls,label=name+' recovered')
        axes[1,0].plot([50,99],[len(s[str(i)]['newFloorNodeIDs']) for i in [50,99]],marker='x',ls=ls,label=name+' new floors')
    for name in ['source','control_99','column_99']:
        p=r['columnProfiles'][name]
        axes[1,1].semilogy([v*1e6 for v in p['distanceAlongColumnM']],p['omega'],marker='.',label=name)
    axes[0,0].axhline(-8,color='k',ls=':',lw=.8)
    axes[0,0].set(xlabel='New local iteration',ylabel='Native log RMS')
    axes[0,1].set(xlabel='New local iteration',ylabel='Native wall heat power [kW]')
    axes[1,0].axhline(146,color='k',ls=':',lw=.8)
    axes[1,0].set(xlabel='Saved local iteration',ylabel='Node count',xticks=[50,99])
    axes[1,1].set(xlabel='Distance along wall column 2 [um]',ylabel='Omega [1/s]')
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Complete 20-degree rocket: near-wall initialization study\nSource local 149; restart local 0-99; no physical qualification',fontsize=11)
    fig.savefig(str(out),dpi=150);plt.close(fig)
    (out.parent/'plot_provenance.json').write_text(json.dumps(dict(source=str(source),sourceSHA256=sha(source),
        scriptSHA256=sha(Path(__file__).resolve()),imageSHA256=sha(out),physicsQualified=False),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
