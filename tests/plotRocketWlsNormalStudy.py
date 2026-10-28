# -*- coding: utf-8 -*-
"""比较同Q的旧增强重启与新五带预条件试算，明确不同预算和真实NS门槛。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,parent=sys.argv[1:3]
assert case.replace('_','').isalnum() and parent.replace('_','').isalnum() and case!=parent
target=VAL/case/'wls_normal_direction_plot.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def rows(value):return [value] if isinstance(value,dict) else list(np.atleast_1d(value))
fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');bindings={};details=[];base=None
for which,source in enumerate([parent,case]):
    file=VAL/source/'wall_velocity_krylov_evidence.mat';bindings[str(file.relative_to(ROOT))]=sha(file)
    r=loadmat(file,simplify_cells=True,variable_names=['result'])['result'];a=rows(r['linearLog'])[0]
    if base is None:base=r['initialQ'].copy()
    else:np.testing.assert_array_equal(base,r['initialQ'])
    rhs=r['initialRate']/r['residualScale'];history=np.atleast_1d(a['gmresResidualHistory'])/np.linalg.norm(rhs)
    label=('Prior graph + SVD retention' if which==0 else 'WLS normal five-band + GMRES')+'\nArnoldi budget '+str(int(r['options']['maximumKrylovIterations']))
    color=['#8d683e','#426a8c'][which];ax.semilogy(np.arange(len(history)),history,color=color,label=label)
    ax.scatter([len(history)-1],[a['fullSystemRelativeDefect']],color=color,marker='x',s=70,zorder=4)
    details.append(dict(case=source,maximumArnoldiBudget=int(r['options']['maximumKrylovIterations']),
        actualOriginalNSDefect=float(a['fullSystemRelativeDefect']),recordedModelDefect=float(history[-1]),
        acceptedUpdates=int(r['acceptedNonlinearIterations']),wlsNormalFiveBandEnabled='wlsNormalCorrection' in r,
        selectedModelActions=int(r['selectedAnalyticApplications']),preconditionerModelActions=int(r['coarsePreconditionerModelApplications']),
        coarseReuseVerificationActions=int(r['coarseReuseVerificationModelApplications']),actualNS=int(r['residualEvaluations'])))
ax.axhline(.1,ls='--',color='#555555',label='Original NS target: 0.1')
ax.set(xlabel='Recorded residual index (different restart policies)',ylabel='Defect / original RHS norm',
       title='Same Q / M5 / 30 km / alpha 20 deg / supplied Tw 300 K\nLines: current model; crosses: original NS total-direction rechecks')
ax.grid(alpha=.2);ax.legend(loc='center left',bbox_to_anchor=(1.02,.5),fontsize=9)
fig.suptitle('UNCONVERGED RESEARCH / successive initial guesses and different budgets / no controlled speedup claim',fontsize=10)
image=VAL/case/'wls_normal_direction_comparison.png';assert not image.exists();fig.savefig(image,dpi=160);plt.close(fig)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='generated-source-bound-same-Q-WLS-normal-direction-comparison',processID=os.getpid(),details=details,
    controlledSpeedupClaimed=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__),
    figure=dict(path=str(image.relative_to(ROOT)),sha256=sha(image)))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
