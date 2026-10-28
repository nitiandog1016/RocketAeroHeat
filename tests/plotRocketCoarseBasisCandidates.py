# -*- coding: utf-8 -*-
"""比较保存真实剩余上的粗基预测与像存储；不冒称原NS接受或物理精度。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
folder=VAL/case;target=folder/'coarse_basis_plot.json';assert not target.exists()
files=[folder/'coarse_candidates_summary.json',VAL/'layercoarserocket01/coarse_layer_summary.json']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();bindings={str(p.relative_to(ROOT)):sha(p) for p in files}
s=json.loads(files[0].read_text(encoding='utf8'));parent=json.loads(files[1].read_text(encoding='utf8'));assert s['kind']=='real'
p=s['primitive'];q=s['spatial'];values=[p['parentProjectedRelativeDefect'],p['candidateProjectedRelativeDefect'],q['candidateProjectedRelativeDefect']]
labels=['Conservative parent\n'+str(s['parentCoarseUnknowns'])+' unknowns','Primitive coupled\n'+str(p['coarseUnknowns'])+' unknowns','Targeted nose spatial\n'+str(q['coarseUnknowns'])+' unknowns']
fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');colors=['#8b969e','#ba8463','#41748f']
ax=axes[0];bars=ax.bar(np.arange(3),values,color=colors);ax.set_xticks(np.arange(3),labels);ax.set_ylim(0,max(values)*1.17);ax.set_ylabel('Predicted remaining norm / original right-hand-side norm');ax.set_title('Same saved original NS linear remaining')
for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value+max(values)*.02,f'{value:.6f}',ha='center',fontsize=9)
ax.grid(True,axis='y',alpha=.2)
sizes=np.asarray([parent['imageBytes'],p['imageBytes'],q['imageBytes']])/1048576;ax=axes[1];bars=ax.bar(np.arange(3),sizes,color=colors);ax.set_xticks(np.arange(3),labels);ax.set_ylim(0,280)
ax.axhline(256,color='#555555',ls='--',label='Per-image bound: 256 MiB');ax.set_ylabel('Sparse coarse image storage (MiB)');ax.set_title('Memory and model construction costs differ')
for bar,value in zip(bars,sizes):ax.text(bar.get_x()+bar.get_width()/2,value+4,f'{value:.2f}',ha='center',fontsize=9)
ax.grid(True,axis='y',alpha=.2);ax.legend(loc='lower right',fontsize=8)
fig.suptitle('M5 / 30 km / alpha 20 deg / supplied Tw 300 K / fixed physical state\nSelected-model projections only; no accepted flow update and no added qualified area',fontsize=11)
png=folder/'coarse_basis_candidates.png';fig.savefig(png,dpi=160);plt.close(fig)
out=dict(status='generated-source-bound-real-coarse-basis-prediction-and-memory-comparison',processID=os.getpid(),scientificPNG=png.name,
    inputSHA256=bindings,scriptSHA256=sha(Path(__file__)),physicalSolution=False,newProductionSolvedAreaM2=0)
assert all(sha(ROOT/p)==h for p,h in bindings.items());target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
