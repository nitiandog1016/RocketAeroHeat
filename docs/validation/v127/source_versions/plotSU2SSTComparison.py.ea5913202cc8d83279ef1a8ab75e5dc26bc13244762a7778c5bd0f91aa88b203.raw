"""绘制已终态的CFL、残差和热功率原始历史，帮助辨别数值推进与物理稳定性。"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepareSU2CheckpointContinuation import ROOT,read,sha


def main():
    base=ROOT/'outputs/v127_validation';dest=ROOT/'docs/validation/v127';dest.mkdir(parents=True,exist_ok=True)
    comparison=read(base/'comparison.json');assert not (dest/'comparison.png').exists()
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True,constrained_layout=True)
    colors={'control':'#28628b','cap05':'#c07727'}
    labels={'control':'Control: CFL cap 5','cap05':'Candidate: CFL cap 0.5'}
    for name in ['control','cap05']:
        r=read(base/name/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0
        with (base/name/'history.csv').open() as f:
            reader=csv.DictReader(f,skipinitialspace=True);rows=[{k.strip():v for k,v in row.items()} for row in reader]
        col=lambda key:np.array([float(row[key]) for row in rows])
        x=col('Inner_Iter');color=colors[name]
        axes[0].semilogy(x,col('Max CFL'),color=color,label=labels[name])
        axes[0].semilogy(x,col('Avg CFL'),color=color,ls=':',alpha=.8)
        axes[1].plot(x,col('rms[Rho]')-col('rms[Rho]')[0],color=color,label=labels[name]+' / density')
        axes[1].plot(x,col('rms[RhoE]')-col('rms[RhoE]')[0],color=color,ls='--',label=labels[name]+' / energy')
        axes[2].plot(x,col('HF')/1000,color=color,label=labels[name])
    axes[0].set_ylabel('Local CFL');axes[0].set_title('Solid: maximum CFL; dotted: mean CFL',fontsize=10)
    axes[1].set_ylabel('Change in log10 residual\nrelative to own local step 0')
    axes[1].axhline(0,color='#777777',lw=.7)
    axes[2].set_ylabel('Native wall heat power / kW');axes[2].set_xlabel('Local iteration after common SST checkpoint 49')
    for ax in axes:
        ax.grid(True,alpha=.2);ax.legend(fontsize=8,loc='best')
        for side in ['top','right']:ax.spines[side].set_visible(False)
    fig.suptitle('Full rocket, M5 / 20 deg / given Tw 300 K\nSST numerical comparison; unconverged research fields',fontsize=12)
    fig.savefig(str(dest/'comparison.png'),dpi=160);plt.close(fig)
    p=dest/'plot_provenance.json';p.write_text(json.dumps(dict(
        sourceSHA256={str(x):sha(x) for x in [Path(__file__).resolve(),base/'control/history.csv',base/'cap05/history.csv',base/'comparison.json']},
        imageSHA256=sha(dest/'comparison.png'),physicsQualified=False,newQualifiedArea=0),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
