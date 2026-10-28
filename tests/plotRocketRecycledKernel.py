# -*- coding: utf-8 -*-
"""同时画Arnoldi步和完整A调用的代价，防止只按较少迭代步宣称增强加速。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
revision=sys.argv[2] if len(sys.argv)>2 else '';assert not revision or revision.isdigit()
folder=VAL/case;source=folder/'recycled_kernel_evidence.mat';target=folder/('recycled_kernel_cost_comparison'+revision+'.png');assert not target.exists()
s=loadmat(source,simplify_cells=True);fig,axes=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
rows=lambda x:[x] if isinstance(x,dict) else list(np.atleast_1d(x))
for key,label in [('plain','No retained directions'),('recycled','SVD-retained directions')]:
    r=s[key];c=rows(r['cycles']);norm=np.linalg.norm(s['b']);steps=[0];actions=[0];residual=[1.]
    for entry in c:
        steps.append(steps[-1]+entry['arnoldiSteps']);actions.append(actions[-1]+entry['modelApplications']);residual.append(entry['actualResidualNorm']/norm)
    axes[0].semilogy(steps,residual,'o-',label=label);axes[1].semilogy(actions,residual,'o-',label=label)
for ax in axes[:2]:
    ax.axhline(1e-9,color='.4',ls='--',lw=1);ax.grid(alpha=.2);ax.set_ylabel('Dense-matrix relative residual');ax.legend(fontsize=8)
axes[0].set(xlabel='Arnoldi steps',title='80-variable algebra test: step count')
axes[1].set(xlabel='All current A applications',title='Include retained-image builds and probes')
r=s['unitReport'];c=rows(r['cycles']);axes[2].semilogy([0]+[e['totalArnoldiSteps'] for e in c],[1.]+[e['actualResidualNorm']/np.linalg.norm(s['unitB']) for e in c],'o-',color='#8b5e34')
axes[2].set(xlabel='Arnoldi steps',ylabel='Current unit-Jacobian relative residual',title='True 32-CV Jacobian\nRandom RHS; not flow convergence');axes[2].grid(alpha=.2)
fig.suptitle('Numerical validation only: fewer steps can still require more model work',fontsize=13);fig.savefig(target,dpi=160);plt.close(fig)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

out=dict(status='generated-source-bound-recycled-kernel-work-and-step-comparison',sourceSHA256=sha(source),scientificPNG=target.name,
    figureSHA256=sha(target),scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0)
(folder/('recycled_kernel_plot'+revision+'.json')).write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
