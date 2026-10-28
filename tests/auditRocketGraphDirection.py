# -*- coding: utf-8 -*-
"""复核真实失败方向的质量项/图作用/原NS差别与逐层误差，并绘制科研诊断图。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1]
revision=sys.argv[2] if len(sys.argv)>2 else '';assert not revision or revision.isdigit()
assert case.replace('_','').isalnum();target=VAL/(case+'_independent_graph_direction'+revision+'.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'graph_direction_evidence.mat';s=loadmat(file,simplify_cells=True);source=Path(s['sourceFile'])
bindings={str(p.relative_to(ROOT)):sha(p) for p in (file,source)};r=loadmat(source,simplify_cells=True,variable_names=['result'])['result'];report=s['report']
log=r['linearLog'];log=log if isinstance(log,dict) else log[0];e=log['steadyDirectionEvidence'];op=r['firstOperator']
np.testing.assert_array_equal(e['directionU'],s['directionU']);np.testing.assert_array_equal(e['value'],s['fullAction'])
np.testing.assert_allclose(s['massAction'],s['directionU']/op['pseudoTimeStepSeconds'],rtol=3e-15,atol=0)
np.testing.assert_allclose(s['graphAction'],s['massAction']-s['preconditionerAction'],rtol=3e-15,atol=0)
rhs=r['initialRate']/r['residualScale'];np.testing.assert_array_equal(rhs,s['rhs']);full=s['fullAction'];graph=s['graphAction'];error=graph-full
values=dict(graphActionRelativeToFullNS=float(np.linalg.norm(error)/np.linalg.norm(full)),
    equationActionRelativeErrors=np.linalg.norm(error,axis=0)/np.linalg.norm(full,axis=0),
    graphSteadyLinearDefect=float(np.linalg.norm(-graph-rhs)/np.linalg.norm(rhs)),
    fullNSSteadyLinearDefect=float(np.linalg.norm(-full-rhs)/np.linalg.norm(rhs)))
for key,value in values.items():np.testing.assert_allclose(value,report[key],rtol=3e-12,atol=1e-15)
nf=int(op['sourceFaces']);nl=int(op['layers']);fraction=np.array([np.sum(error[i*nf:(i+1)*nf]**2)/np.sum(error**2) for i in range(nl)])
np.testing.assert_allclose(fraction,report['layerErrorFraction'],rtol=3e-12,atol=1e-30)
order=np.argsort(np.sum(error**2,axis=1))[::-1][:20];np.testing.assert_array_equal(order+1,report['largestErrorCellIDs'])
np.testing.assert_array_equal(order%nf+1,report['largestErrorOriginalFaceIDs'])
history=np.asarray(log['gmresResidualHistory'])/np.linalg.norm(rhs);np.testing.assert_allclose(history,report['normalizedGMRESHistory'],rtol=3e-12,atol=0)
assert r['acceptedNonlinearIterations']==0 and not r['steadyConverged'] and report['newOriginalNSEvaluations']==0
fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained');layers=np.arange(1,nl+1)
axes[0].bar(layers,100*fraction,color='#426a8c');axes[0].set(xlabel='Wall-normal layer',ylabel='Share of squared action error (%)',xticks=layers,
    title='Approximate graph action versus full NS')
axes[0].text(.97,.92,'First two layers: %.3f%%'%(100*sum(fraction[:2])),ha='right',transform=axes[0].transAxes)
axes[1].semilogy(np.arange(len(history)),history,'o-',ms=3,color='#426a8c',label='Actual GMRES model history')
axes[1].axhline(r['options']['linearRelativeTolerance'],ls='--',color='#a65e35',label='Required full NS target: 0.1')
axes[1].set(xlabel='Krylov directions',ylabel='Relative linear residual',ylim=(.08,1.12),title='24 directions: target not met')
axes[1].legend(loc='center left');axes[1].grid(alpha=.2)
fig.suptitle('M5 / 30 km / alpha 20 deg / UNCONVERGED RESEARCH\nSingle saved direction; action error is not pressure or heat-flux accuracy',fontsize=11)
figure=VAL/case/('graph_direction_diagnostic'+revision+'.png');assert not figure.exists();fig.savefig(figure,dpi=160);plt.close(fig)
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-failed-graph-direction-arithmetic-and-layer-diagnosis',case=case,processID=os.getpid(),
    graphActionRelativeToFullNS=values['graphActionRelativeToFullNS'],graphSteadyLinearDefect=values['graphSteadyLinearDefect'],
    fullNSSteadyLinearDefect=values['fullNSSteadyLinearDefect'],firstTwoLayersSquaredErrorFraction=float(sum(fraction[:2])),
    newOriginalNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__),
    figure=dict(path=str(figure.relative_to(ROOT)),sha256=sha(figure)))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
