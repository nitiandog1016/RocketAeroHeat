# -*- coding: utf-8 -*-
"""用已逐MAT核对的工作台候选数据绘制同方向坐标对照，保留真实拒绝和未收敛身份。"""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=VAL/'primitivereservoir01'
target=folder/'trial_coordinate_comparison.png';receipt=folder/'trial_coordinate_figure_receipt.json'
assert not target.exists() and not receipt.exists()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();bindings={};cases=[]
for name,label,color in [('workbench35','Conservative straight path','C0'),('workbench36','Primitive tangent path','C1')]:
    path=VAL/name/'wall-velocity-study-data.json';bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    d=json.loads(path.read_text(encoding='utf8'))['cases'][2];tr=d['solverDiagnostics']['nonlinearTrials']
    if isinstance(tr,dict):tr=[tr]
    assert len(tr)==13 and sum(x['accepted'] for x in tr)==1
    cases.append((d,tr,label,color))
fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
for d,tr,label,color in cases:
    valid=[x for x in tr if x['candidateRmsMerit'] is not None]
    alpha=np.array([x['alpha'] for x in valid]);base=tr[0]['baseRmsMerit'];initial_max=d['solverDiagnostics']['segmentResidualTrend']['initialMaximumResidual']
    for ax,key,denominator in [(axes[0],'candidateRmsMerit',base),(axes[1],'candidateMaximumResidual',initial_max)]:
        values=np.array([x[key] for x in valid])/denominator
        ax.plot(alpha,values,'o-',label=label,color=color)
        accept=[i for i,x in enumerate(valid) if x['accepted']]
        ax.scatter(alpha[accept],values[accept],s=110,facecolors='none',edgecolors=color,linewidths=1.6,zorder=4)
        ax.set_xscale('log',base=2);ax.grid(True,alpha=.25);ax.axhline(1,color='.4',linestyle='--',linewidth=1)
for ax in axes:
    ax.set_xticks([1/4096,1/2048,1/1024],['1/4096','1/2048','1/1024']);ax.set_xlabel('Actual candidate alpha (dimensionless)')
axes[0].set_ylabel('RMS / initial RMS');axes[0].set_title('Original RMS descent test')
axes[1].set_ylabel('Maximum residual / initial maximum');axes[1].set_yscale('log');axes[1].set_title('Largest local residual, same initial state')
axes[0].legend(fontsize=9);fig.suptitle('Same Q20 and identical original-NS direction / UNCONVERGED RESEARCH\nCircles with outer ring: accepted. Dashed line: initial value. Physical time = 0',fontsize=11)
notes=[]
for d,tr,label,_ in cases:
    guards=sum(bool(x['errorIdentifier']) for x in tr)
    notes.append(label+': '+str(guards)+' state-rejected trials; '+str(d['summary']['residualEvaluations'])+' main NS calls')
fig.supxlabel('Only full-NS-evaluated candidates plotted; no value is assigned to failed states.\n'+' | '.join(notes),fontsize=9)
fig.savefig(target,dpi=160);plt.close(fig)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='rendered-actual-same-direction-trial-coordinate-comparison',inputSHA256=bindings,figureSHA256=sha(target),
         scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0,browserVisualOrWebGLRuntimeVerified=False)
receipt.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
