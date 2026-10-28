# -*- coding: utf-8 -*-
"""从已审计的固定方向JSON绘制前后误差比较；不重新求解或替换物理误差。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=VAL/'wlsrocket01'
source=folder/'wls_direction_real_summary.json';prior=VAL/'graphdirection01/graph_direction_summary.json'
report=json.loads(source.read_text(encoding='utf8'));old=json.loads(prior.read_text(encoding='utf8'));assert report['sourceCase']==old['sourceCase']
revision=sys.argv[1] if len(sys.argv)>1 else '';assert not revision or revision.isdigit()
datafile=folder/'wls_direction_real_evidence.mat';graphfile=VAL/'graphdirection01/graph_direction_evidence.mat'
data=loadmat(datafile,simplify_cells=True,variable_names=['allEvidence','scale']);item=data['allEvidence'];ref=item['refinements'][-1]
full=(ref['ratePlus']-ref['rateMinus'])/(2*ref['h'])/data['scale'];graph=loadmat(graphfile,simplify_cells=True,variable_names=['graphAction'])['graphAction']
graph_errors=np.linalg.norm(graph-full,axis=0)/np.linalg.norm(full,axis=0)
analytic_errors=np.linalg.norm(item['euler']+item['viscous']-full,axis=0)/np.linalg.norm(full,axis=0)
np.testing.assert_allclose(analytic_errors,report['equationL2RelativeDifferences'],rtol=2e-9,atol=0)
figure=folder/('wls_direction_comparison'+revision+'.png');target=folder/('wls_direction_figure'+revision+'.json');assert not figure.exists() and not target.exists()
fig,ax=plt.subplots(figsize=(12,4.4),layout='constrained');x=np.arange(5);width=.36
ax.bar(x-width/2,graph_errors,width,color='#8a969f',label='Two-point lateral graph')
ax.bar(x+width/2,analytic_errors,width,color='#426a8c',label='Full WLS + selected boundary direction')
ax.set(yscale='log',ylim=(1e-9,2),xticks=x,xticklabels=['Mass','X momentum','Y momentum','Z momentum','Energy'],
    ylabel='Relative action error (L2 per equation)',title='Same saved krylov07 direction / half-step central NS reference')
ax.legend(loc='upper left',bbox_to_anchor=(1.01,1));ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
fig.suptitle('M5 / 30 km / alpha 20 deg / prescribed Tw = 300 K / UNCONVERGED RESEARCH\nDirection consistency is not pressure or heat-flux accuracy; no state update',fontsize=11)
fig.savefig(figure,dpi=160);plt.close(fig)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
out=dict(status='rendered-from-verified-fixed-direction-summary-data',processID=os.getpid(),figure=str(figure.relative_to(ROOT)),figureSHA256=sha(figure),
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (source,prior,datafile,graphfile)},scriptSHA256=sha(__file__),
    graphEquationRelativeErrors=graph_errors.tolist(),analyticEquationRelativeErrors=analytic_errors.tolist(),referenceKind='same-half-step-central-NS',
    newOriginalNSEvaluations=0,newProductionSolvedAreaM2=0)
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
