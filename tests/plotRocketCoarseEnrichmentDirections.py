# -*- coding: utf-8 -*-
"""同一真实端点的90/210维粗空间外层方向对照，区分模型历史与原NS复核。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
previous=sys.argv[2] if len(sys.argv)>2 else 'krylov09';assert previous.replace('_','').isalnum() and previous!=case
revision=sys.argv[3] if len(sys.argv)>3 else '';assert not revision or revision.isdigit()
folder=VAL/case;target=folder/('coarse_enrichment_direction_plot'+revision+'.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
sources=[VAL/previous/'wall_velocity_krylov_evidence.mat',folder/'wall_velocity_krylov_evidence.mat'];bindings={str(p.relative_to(ROOT)):sha(p) for p in sources}
fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained');checks=[];base=None
for file in sources:
    r=loadmat(file,simplify_cells=True,variable_names=['result'])['result'];a=r['linearLog'];a=a if isinstance(a,dict) else list(a)[0]
    if base is None:base=r['initialQ'].copy()
    else:np.testing.assert_array_equal(base,r['initialQ'])
    assert r['options']['maximumKrylovIterations']>0 and r['options']['graphSweeps']==1 and r['options']['initialPseudoDt']==1e-7
    history=np.atleast_1d(a['gmresResidualHistory'])/np.linalg.norm(r['initialRate']/r['residualScale'])
    label=file.parent.name+' / '+str(int(r['coarseUnknowns']))+' coarse / budget '+str(int(r['options']['maximumKrylovIterations']));ax.semilogy(np.arange(len(history)),history,label=label)
    ax.scatter([len(history)-1],[a['fullSystemRelativeDefect']],marker='x',s=55)
    checks.append(dict(case=file.parent.name,coarseUnknowns=int(r['coarseUnknowns']),actualNSDirectionDefect=float(a['fullSystemRelativeDefect']),
        modelDirectionDefect=float(history[-1]),acceptedFlowUpdates=int(r['acceptedNonlinearIterations']),maximumArnoldiBudget=int(r['options']['maximumKrylovIterations'])));del r,a
ax.axhline(.1,color='#555555',ls='--',label='Original NS acceptance target 0.1')
ax.set(xlabel='GMRES recorded iteration',ylabel='Relative linear defect',title='Same real Q / M5 / 30 km / alpha 20 deg / supplied Tw 300 K\nLines: selected model; crosses: original NS combined-direction rechecks')
ax.grid(True,alpha=.2);ax.legend(loc='center left',bbox_to_anchor=(1.02,.5),fontsize=9)
png=folder/('coarse_enrichment_directions'+revision+'.png');fig.savefig(png,dpi=160);plt.close(fig)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='generated-source-bound-same-endpoint-coarse-enrichment-direction-comparison',processID=os.getpid(),details=checks,scientificPNG=png.name,
    inputSHA256=bindings,scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0)
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
