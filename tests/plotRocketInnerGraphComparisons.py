# -*- coding: utf-8 -*-
"""科研比较：近似求逆精度、完整稳态模型与单方向预测严格分列。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';bcase,ccase=sys.argv[1:3]
revision=sys.argv[3] if len(sys.argv)>3 else '';assert not revision or revision.isdigit()
assert all(x.replace('_','').isalnum() for x in (bcase,ccase))
folder=VAL/ccase;target=folder/('inner_steady_graph_figures'+revision+'.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path):
    path=Path(path);bindings[str(path.relative_to(ROOT))]=sha(path);return json.loads(path.read_text(encoding='utf8'))
def rows(x):return [x] if isinstance(x,dict) else x
b=read(VAL/bcase/'inner_graph_real_summary.json');c=read(folder/'steady_graph_real_summary.json')
entries=list(rows(b['summaries']));labels=['B: one scan','B: inner solve'];colors=['#7894aa','#426a8c']
rejected=c['status'].startswith('rejected')
if not rejected:
    entries+=rows(c['summaries']);labels+=['C: one scan','C: inner solve'];colors+=['#b28b69','#926135']
fig,axes=plt.subplots(1,3,figsize=(13,4.6),layout='constrained');x=np.arange(len(entries))
approx=[entry.get('relativeBResidual',entry.get('relativeCResidual')) for entry in entries]
full=[entry['afterCoarseRelativeResidual'] for entry in entries];predicted=[entry['predictedTotalDefectAfterBestScalar'] for entry in entries]
for ax,values,title,ylabel in zip(axes,[approx,full,predicted],
        ['Approximate system solve','Full A after coarse correction','Best scalar MODEL prediction'],
        ['Residual / current remaining norm','Residual / current remaining norm','Predicted residual / original RHS norm']):
    ax.bar(x,values,color=colors,width=.65);ax.set(title=title,xticks=x,xticklabels=labels,ylabel=ylabel)
    ax.tick_params(axis='x',labelrotation=28,labelsize=8);ax.grid(axis='y',alpha=.2)
    for i,value in enumerate(values):ax.text(i,value,'%.4g'%value,ha='center',va='bottom',fontsize=8)
axes[0].set_yscale('log');axes[0].axhline(1e-3,color='#555555',ls='--',lw=1,label='Inner target 1e-3');axes[0].legend(fontsize=8)
axes[1].set_yscale('log');axes[1].axhline(1,color='#555555',ls=':',lw=1)
axes[1].set_title('Full A after coarse correction\nDotted line: unchanged (1)',fontsize=10)
axes[2].axhline(.1,color='#555555',ls='--',lw=1,label='Required original NS target');axes[2].legend(loc='lower center',fontsize=8)
fig.suptitle('M5 / 30 km / alpha 20 deg / same Q / UNCONVERGED RESEARCH\nB includes I/tau; C removes it for all five fields. Predictions are not original NS acceptance.',fontsize=11)
if rejected:fig.text(.5,.005,'C factor rejected: '+c['errorIdentifier']+'; no C direction or flow result',ha='center',fontsize=9)
image=folder/('inner_steady_graph_comparison'+revision+'.png');assert not image.exists();fig.savefig(image,dpi=160);plt.close(fig)
small=read(VAL/'innergraphunit01/dense_schur_references.json');results=small['candidates'];fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
labels2=['All shifted','Density unshifted','Energy unshifted','Density+energy\nunshifted','All unshifted']
axes[0].bar(np.arange(len(results)),[entry['steadyAPreconditionedCondition'] for entry in results],color='#426a8c');axes[0].set(yscale='log',ylabel='2-norm condition of A times exact C inverse',title='Entire small system / 32 CV / 160 unknowns')
axes[1].bar(np.arange(len(results)),[entry['bestScalarModelRelativeDefect'] for entry in results],color='#926135');axes[1].set(ylabel='Best scalar model defect / RHS norm',title='Small-domain reference only')
for ax in axes:ax.set(xticks=np.arange(len(results)),xticklabels=labels2);ax.tick_params(axis='x',labelrotation=25,labelsize=8);ax.grid(axis='y',alpha=.2)
fig.suptitle('Dense Schur/time-diagonal references / NOT a real-rocket or physical validation\nPartial removal is not assumed better; no physical equation or mesh changed.',fontsize=11)
smallimage=folder/('dense_schur_reference_comparison'+revision+'.png');assert not smallimage.exists();fig.savefig(smallimage,dpi=160);plt.close(fig)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='generated-source-bound-inner-and-steady-graph-scientific-comparisons',processID=os.getpid(),
    figures=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in (image,smallimage)],CfactorRejected=rejected,
    newModelApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
