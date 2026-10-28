# -*- coding: utf-8 -*-
"""用保存的实际拒绝摘要绘制科研诊断图；不把内部线性剩余当热流误差。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1] if len(sys.argv)>1 else 'block_ilu_figures01';assert name.replace('_','').isalnum();folder=VAL/name;assert not folder.exists();folder.mkdir()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
paths=[VAL/'block_linear_research_summary.json',VAL/'blockilufailurediag01/block_ilu_failure_diagnostic_summary.json',VAL/'blockfillreal01/block_fill1_real_summary.json',VAL/'blockilureal01_independent_real_C_and_factors.json',VAL/'blockfillreal01_independent_real_fill1_factors.json']
summary,zero,one,a0,a1=[json.loads(p.read_text(encoding='utf-8-sig')) for p in paths];studies=summary['studies'];labels=['ILU0 / default','ILU0 / scaled','ILU1 / scaled']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');colors=['#355f8d','#238a8d','#56a868']
ax=axes[0,0];memory=np.array([379080388,one['enrichmentReport']['persistentFactorStorageUpperBoundBytes']])/1048576
ax.bar(['ILU0','ILU1'],memory,color=[colors[0],colors[2]],width=.5);ax.axhline(768,color='0.4',linestyle='--',label='Core-array budget: 768 MiB');ax.set_ylim(0,900);ax.set_ylabel('Factor core arrays / MiB');ax.set_title('Storage bound (not process peak)');ax.legend(frameon=False)
for i,v in enumerate(memory):ax.text(i,v+12,f'{v:.1f}',ha='center')
ax=axes[0,1];errors=[a0['fullLUApproximationRelativeDifference'],a1['fullLUApproximationRelativeDifference']];ax.bar(['ILU0','ILU1'],errors,color=[colors[0],colors[2]],width=.5);ax.set_ylabel(r'$\|C-LU\|_F/\|C\|_F$');ax.set_title('More fill did not improve the full approximation');ax.set_ylim(0,7.5)
for i,v in enumerate(errors):ax.text(i,v+.12,f'{v:.3f}',ha='center')
ax=axes[1,0];x=np.arange(3);w=.33
ax.bar(x-w/2,[s['basisOrthogonalityError'] for s in studies],width=w,label='Basis orthogonality',color=colors[0]);ax.bar(x+w/2,[s['modelRelationRelativeDifference'] for s in studies],width=w,label='Combined image relation',color='#cb7c41');ax.axhline(1e-8,color='0.35',linestyle='--',label='Original gate');ax.set_yscale('log');ax.set_ylim(1e-13,1e3);ax.set_xticks(x,labels,rotation=12);ax.set_ylabel('Dimensionless error');ax.set_title('Every candidate was rejected');ax.legend(frameon=False,fontsize=9,loc='upper right')
ax=axes[1,1];remaining=[100*s['innerRelativeCResidual'] for s in studies];ax.bar(labels,remaining,color=colors,width=.55);ax.set_ylim(0,115);ax.set_ylabel('Internal RHS remaining / %');ax.set_title('At least 99.88% of the internal RHS remained');ax.tick_params(axis='x',labelrotation=12)
for i,v in enumerate(remaining):ax.text(i,v+2,f'{v:.4f}%',ha='center')
fig.suptitle('Same-state five-conservation linear studies\n0 new NS evaluations / 0 physical updates / 0 new qualified area',fontsize=13)
png=folder/'block_ilu_rejected_studies.png';pdf=folder/'block_ilu_rejected_studies.pdf';fig.savefig(png,dpi=160);fig.savefig(pdf);plt.close(fig)
out=dict(status='generated-scientific-diagnostic-figures-from-source-bound-rejected-studies',processID=os.getpid(),visuallyInspected=False,
         newActualNSEvaluations=0,physicalSolution=False,inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in paths},
         outputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (png,pdf)},scriptSHA256=sha(__file__))
(folder/'figure_manifest.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
