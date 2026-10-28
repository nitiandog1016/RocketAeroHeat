# -*- coding: utf-8 -*-
"""按原始混合筛选历史绘制投影收益与实际稀疏像字节，分清模型预测与NS验收。"""
from pathlib import Path
import hashlib,json,os,sys,gc
from scipy.io import loadmat
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
folder=VAL/case;target=folder/'mixed_selection_plot.json';assert not target.exists()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=folder/'mixed_coarse_evidence.mat';digest=sha(file);s=loadmat(file,simplify_cells=True,variable_names=['coarse','report'])
history=s['coarse']['selectionHistory'].copy();report=s['report'];budget=s['coarse']['maximumImageBytes'];del s;gc.collect()
fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
axes[0].plot(history[:,0],history[:,6]*report['originalRemainingRelativeDefect'],'o-',ms=3,color='#386e8b')
axes[0].axhline(report['baseProjectedRelativeDefect'],ls='--',color='#777777',label='Conservative base projection')
axes[0].set(xlabel='Selected additional primitive columns',ylabel='Predicted remaining norm / original RHS norm',title='Fixed saved original NS remaining; model projection only')
axes[0].grid(True,alpha=.2);axes[0].legend(fontsize=8)
axes[1].plot(history[:,0],history[:,5]/1048576,'o-',ms=3,color='#386e8b');axes[1].axhline(budget/1048576,ls='--',color='#555555',label='Original image bound 256 MiB')
axes[1].set(xlabel='Selected additional primitive columns',ylabel='Sparse image storage (MiB)',title='Whole columns retained; no small-value truncation')
axes[1].grid(True,alpha=.2);axes[1].legend(loc='lower right',fontsize=8)
fig.suptitle('M5 / 30 km / alpha 20 deg / prescribed Tw 300 K / unchanged physical state\nMixed fixed column lifts; no added qualified surface area',fontsize=11)
png=folder/'mixed_coarse_selection.png';fig.savefig(png,dpi=160);plt.close(fig);assert sha(file)==digest
out=dict(status='generated-source-bound-mixed-selection-gain-and-image-memory-figure',processID=os.getpid(),scientificPNG=png.name,
    inputSHA256={str(file.relative_to(ROOT)):digest},scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0)
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
