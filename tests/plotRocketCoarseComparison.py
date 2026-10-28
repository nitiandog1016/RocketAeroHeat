# -*- coding: utf-8 -*-
"""从原始研究摘要及实际GMRES历史绘制粗校正对照；不是物理精度或浏览器验收。"""
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
folder=VAL/case;target=folder/('coarse_comparison'+revision+'.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
sources=[VAL/'coarsepost01/coarse_post_summary.json',VAL/'coarselsrocket01/least_squares_coarse_summary.json',
         VAL/'krylov08/wall_velocity_krylov_evidence.mat',folder/'wall_velocity_krylov_evidence.mat']
bindings={str(p.relative_to(ROOT)):sha(p) for p in sources}
post=json.loads(sources[0].read_text(encoding='utf8'));ls=json.loads(sources[1].read_text(encoding='utf8'))
fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained');ax=axes[0]
values=[ls['initialFineModelRelativeResidual']]+[x[0] for x in post['history']]+[ls['correctedFineModelRelativeResidual']]
names=['Graph scan','Galerkin','+1 post','+2 posts','+3 posts','+4 posts','Fine-L2 coarse']
bars=ax.bar(np.arange(len(values)),values,color=['#8a969e']+['#bb7458']*5+['#346b8c'])
ax.set_yscale('log');ax.set_xticks(np.arange(len(values)),names,rotation=25,ha='right');ax.set_ylabel('Fine-model remaining norm / original RHS norm')
ax.set_title('Same Q and RHS: coarse algebraic success is insufficient');ax.grid(True,axis='y',alpha=.2)
for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value*1.08,f'{value:.3g}',ha='center',fontsize=8)
ax.set_ylim(min(values)*.65,max(values)*2)
ax=axes[1];details=[]
for source,label in zip(sources[2:],['v83 graph / krylov08','v83 graph + fine-L2 coarse / '+case]):
    s=loadmat(source,simplify_cells=True,variable_names=['result']);r=s['result'];logs=r['linearLog'];logs=[logs] if isinstance(logs,dict) else list(logs);a=logs[0]
    hist=np.atleast_1d(a['gmresResidualHistory']);rhs=np.linalg.norm(r['initialRate']/r['residualScale']);ax.semilogy(np.arange(len(hist)),hist/rhs,label=label)
    ax.scatter([len(hist)-1],[a['fullSystemRelativeDefect']],marker='x',s=45)
    details.append(dict(case=source.parent.name,modelFinalRelativeResidual=float(hist[-1]/rhs),actualNSRelativeDefect=float(a['fullSystemRelativeDefect']),
                        acceptedIterations=int(r['acceptedNonlinearIterations']),actualNSCalls=int(r['residualEvaluations'])))
    del s,r
ax.axhline(.1,color='#444444',linestyle='--',label='Original NS target 0.1');ax.set(xlabel='GMRES recorded iteration',ylabel='Relative linear residual',title='Model history; crosses show original NS recheck')
ax.grid(True,alpha=.2);ax.legend(fontsize=8,loc='center right')
fig.suptitle('M5 / 30 km / alpha 20 deg / prescribed Tw 300 K / same unchanged real endpoint\nNumerical direction study only; no added qualified surface area',fontsize=11)
png=folder/('coarse_correction_comparison'+revision+'.png');fig.savefig(png,dpi=160);plt.close(fig)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
report=dict(status='generated-source-bound-coarse-model-and-original-NS-direction-comparison',processID=os.getpid(),details=details,
    fineModelResidualRatios=values,scientificPNG=png.name,inputSHA256=bindings,scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0)
target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))
