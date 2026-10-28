"""绘制同预算方向、真实推进及限制界面非线性温度，不代替浏览器验收。"""
from pathlib import Path
import os,json
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from riemannReservoirReference import ROOT,exact,sha
OUT=ROOT/'outputs/v86_validation';folder=OUT/'comparison01'


def main():
    report=json.loads((folder/'temperature_projection_comparison.json').read_text('utf8'))
    s=loadmat(folder/'temperature_projection_directions.mat',simplify_cells=True,variable_names=['rows','free','projected'])
    rows=s['rows'];gas=rows['gas'];scale=rows['scale'];normal=np.asarray(rows['normal']).reshape(-1,3)[0]
    original=np.row_stack([rows['QL'],rows['QR']]);cells=[int(np.atleast_1d(rows['ownerCV'])[0])-1,int(np.atleast_1d(rows['neighbourCV'])[0])-1]
    alphas=[];merits=[];initial=[];eta=[]
    for branch in report['branches']:
        eta.append(branch['probe']['originalNSDirectionDefect'] if branch['probe'] else np.nan)
        m=branch['main']
        if m:
            r=loadmat(Path(branch['resultFolder'])/'wall_velocity_krylov_evidence.mat',simplify_cells=True,variable_names=['result'])['result']
            trace=r['trace'];trace=[trace] if isinstance(trace,dict) else list(np.atleast_1d(trace))
            alphas.append(float(trace[-1]['alpha']) if trace else 0.);merits.append(m['rmsMerit']);initial.append(m['initialRmsMerit'])
        else:alphas.append(0.);merits.append(np.nan);initial.append(np.nan)
    fig,ax=plt.subplots(1,3,figsize=(14,4.5),layout='constrained');names=['Free','Temperature\nprojected'];x=np.arange(2);colors=['#586d80','#267aaa']
    c=[b['linearReport']['relativeCResidual'] for b in report['branches']]
    ax[0].bar(x-.17,c,.32,label='Approximate C');ax[0].bar(x+.17,eta,.32,label='Fresh original NS')
    ax[0].set_yscale('log');ax[0].set_xticks(x,names);ax[0].set_ylabel('Relative direction defect')
    ax[0].axhline(.1,color='#b44d39',ls='--',lw=1,label='Original NS target 0.1')
    ax[0].axhline(.001,color='#777777',ls=':',lw=1,label='Internal C target 0.001');ax[0].legend(fontsize=8,loc='best')
    ax[0].set_title('Same state / 128 LSQR steps each',fontsize=11)
    improvement=[100*(a-b)/a if np.isfinite(a) else 0. for a,b in zip(initial,merits)]
    ax[1].bar(x,improvement,color=colors);ax[1].set_xticks(x,names);ax[1].set_ylabel('Actual RMS reduction (%)')
    ymax=max([.01]+improvement);ax[1].set_ylim(0,1.35*ymax)
    for i in range(2):ax[1].text(i,improvement[i]+.035*ymax,f'{improvement[i]:.5g}%\nalpha={alphas[i]:.7g}',ha='center',fontsize=9)
    ax[1].set_title('Original NS accepted update only',fontsize=11)
    maximum=max([1/32768]+alphas)*2;sample_alphas=np.linspace(0,maximum,101)
    def primitive(q):
        rho=q[0];v=q[1:4]/rho;temp=(q[4]-.5*rho*(v@v))/(rho*gas['R']/(gas['gamma']-1));return np.r_[rho,v,temp]
    for i,key in enumerate(('free','projected')):
        direction=s[key][cells]*scale;temperatures=[]
        for alpha in sample_alphas:
            q=original+alpha*direction
            if np.any(q[:,0]<=0):temperatures.append(np.nan);continue
            l,r=map(primitive,q)
            if min(l[4],r[4])<=0:temperatures.append(np.nan);continue
            try:sample,_,_=exact(l,r,normal,gas);temperatures.append(sample[4])
            except (AssertionError,ValueError):temperatures.append(np.nan)
        ax[2].plot(sample_alphas,temperatures,color=colors[i],label=names[i].replace('\n',' '))
        if alphas[i]>0:
            sample,_,_=exact(*[primitive(q) for q in original+alphas[i]*direction],normal,gas)
            ax[2].plot(alphas[i],sample[4],'o',color=colors[i])
    ax[2].axhline(2000,color='#b44d39',ls='--',lw=1,label='Unchanged 2000 K limit');ax[2].set_xlabel('Step alpha along saved direction')
    ax[2].set_ylabel('Sample temperature at interior face 19207 (K)');ax[2].legend(fontsize=8,loc='best')
    ax[2].ticklabel_format(axis='x',style='sci',scilimits=(0,0));ax[2].set_title('Public Riemann arithmetic reference',fontsize=11)
    for a in ax:a.grid(True,axis='y',alpha=.2)
    fig.suptitle('M5 / 30 km / alpha=20 deg / given Tw=300 K / UNCONVERGED RESEARCH\nDirection restriction is not a physical temperature boundary; qualified area added = 0',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(OUT/f'temperature_projection_comparison.{ext}',dpi=160)
    plt.close(fig)
    receipt=dict(status='plotted-actual-same-budget-directions-updates-and-local-interface-temperature',
                 acceptedAlpha=alphas,rmsReductionPercent=improvement,sourceFileSHA256=sha(folder/'temperature_projection_directions.mat'),
                 newProductionSolvedAreaM2=0,browserVisualOrWebGLRuntimeVerified=False,scriptSHA256=sha(__file__))
    (OUT/'scientific_plot.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8');print(json.dumps(receipt))


if __name__=='__main__':main()
