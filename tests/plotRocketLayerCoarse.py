# -*- coding: utf-8 -*-
"""显示粗基丰富前后的模型投影与原NS剩余所在层；预测不冒称实际流场改进。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
folder=VAL/case;target=folder/'layer_coarse_plot.json';assert not target.exists()
files=[folder/'coarse_layer_summary.json',VAL/'krylov09/coarse_direction_diagnostic.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in files};s=json.loads(files[0].read_text(encoding='utf8'));d=json.loads(files[1].read_text(encoding='utf8'))
assert s['kind']=='real';fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
shares=np.asarray(d['layerRemainingSquaredNormFractions'])*100;ax=axes[0];ax.bar(np.arange(1,len(shares)+1),shares,color='#527b95')
ax.set(xlabel='Wall-normal layer number (original CVs)',ylabel='Share of original NS linear remaining squared norm (%)',
    title='Saved failed direction: layer 3 contributes most');ax.set_xticks(np.arange(1,len(shares)+1));ax.grid(True,axis='y',alpha=.2)
values=[s['failedRemainingNormRelativeToOriginalRHS'],s['parentPredictedRemainingRelativeToOriginalRHS'],s['enrichedPredictedRemainingRelativeToOriginalRHS']]
labels=['Saved original NS\nremaining','Existing coarse\nmodel projection','Layer-enriched coarse\nmodel projection'];ax=axes[1];bars=ax.bar(np.arange(3),values,color=['#8b949b','#6c91a9','#2d6588'])
ax.set_xticks(np.arange(3),labels);ax.set(ylabel='Remaining norm / original right-hand-side norm',title='Last two bars are model predictions only',ylim=(0,1.18*max(values)))
for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+.015*max(values),f'{v:.6f}',ha='center',fontsize=9)
ax.axhline(.1,color='#555555',ls='--',label='Original NS acceptance target 0.1');ax.grid(True,axis='y',alpha=.2);ax.legend(loc='center right',fontsize=8)
fig.suptitle('M5 / 30 km / alpha 20 deg / supplied Tw 300 K / unchanged physical state\nNested connected volume coarse spaces; no new accepted flow update',fontsize=11)
png=folder/'layer_coarse_projection.png';fig.savefig(png,dpi=160);plt.close(fig)
out=dict(status='generated-source-bound-layer-coarse-projection-figure',processID=os.getpid(),scientificPNG=png.name,inputSHA256=bindings,
    scriptSHA256=sha(Path(__file__)),physicalSolution=False,newProductionSolvedAreaM2=0)
assert all(sha(ROOT/p)==h for p,h in bindings.items());target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
