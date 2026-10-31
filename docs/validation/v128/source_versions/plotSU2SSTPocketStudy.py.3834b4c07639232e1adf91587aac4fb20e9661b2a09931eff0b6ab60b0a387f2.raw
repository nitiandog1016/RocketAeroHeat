"""用已保存空间证据与原始短程历史展示尾缘ω停滞；不掩盖未收敛状态。"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    source=ROOT/'outputs/v128_validation/pocket_comparison.json';r=read(source)
    dest=ROOT/'docs/validation/v128';p=dest/'pocket_comparison.png';assert not p.exists()
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),constrained_layout=True)
    labels={'control':'Unchanged initial state','repaired':'Positive omega initial guess'}
    colors={'control':'#28628b','repaired':'#b16b21'}
    for name in ['control','repaired']:
        h=r['variants'][name]['history'];color=colors[name];profile=r['columnProfiles'][name]
        axes[0].plot(np.arange(10),h['rms[w]'],label=labels[name],color=color,marker='.',ms=4)
        axes[1].plot(np.arange(10),h['rms[RhoE]'],label=labels[name],color=color,marker='.',ms=4)
        d=np.array(profile['distanceAlongColumnM'])*1e6
        axes[2].semilogy(d,profile['omega'],label=labels[name]+' / final',color=color,marker='.',ms=4)
    axes[0].set_ylabel('Native log10 RMS(omega RHS)');axes[1].set_ylabel('Native log10 RMS(energy RHS)')
    for a in axes[:2]:a.set_xlabel('New local iteration (source: 149)')
    axes[2].axhline(1e-4,color='#777777',ls=':',label='Native numerical floor')
    axes[2].set_xlabel('Distance along tail column 2 / micrometre');axes[2].set_ylabel('Omega / 1/s')
    for a in axes:
        a.grid(True,alpha=.2);a.legend(fontsize=7)
        for side in ['top','right']:a.spines[side].set_visible(False)
    result='passed' if r['allScreeningGatesPassed'] else 'FAILED'
    fig.suptitle('Full rocket, M5 / 20 deg / given Tw 300 K — short initial-state study\nScreening '+result+'; unconverged; newly qualified area = 0',fontsize=12)
    fig.savefig(str(p),dpi=160);plt.close(fig)
    (dest/'plot_provenance.json').write_text(json.dumps(dict(sourceSHA256={str(source):sha(source),str(ROOT/'scripts/plotSU2SSTPocketStudy.py'):sha(ROOT/'scripts/plotSU2SSTPocketStudy.py')},
        imageSHA256=sha(p),physicsQualified=False),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
