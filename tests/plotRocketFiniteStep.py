"""同一48单元核心、固定Q24邻域、相同六步预算的原始残差科学图。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
s=loadmat(OUT/'patch02/patch_comparison.mat',simplify_cells=True)
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':10,'axes.unicode_minus':False})
fig,axes=plt.subplots(1,3,figsize=(14.5,5),layout='constrained')
for name,label,color in [('baseline','精确Newton局部基线','#7b8794'),('guarded','最小原始变量幅度','#126caa')]:
    rows=s[name]['summary']['history'];its=[x['iteration'] for x in rows]
    for ax,key in zip(axes[:2],['rms','maximum']):ax.semilogy(its,[x[key] for x in rows],'-o',label=label,color=color,lw=2,ms=4)
matched=loadmat(OUT/'patch03/patch_comparison.mat',simplify_cells=True)['baseline']['summary']['history']
for ax,key in zip(axes[:2],['rms','maximum']):ax.semilogy([x['iteration'] for x in matched],[x[key] for x in matched],'-s',label='同0.0999约束：最小守恒变量幅度',color='#8d609e',lw=1.8,ms=4)
for ax,title in zip(axes[:2],['原始 RMS / (1/s)','最大局部残差 / (1/s)']):
    ax.set(title=title,xlabel='实际连续接受步');ax.grid(alpha=.2)
axes[0].legend(frameon=False,loc='upper left',bbox_to_anchor=(0,-.18),fontsize=9)
d=loadmat(OUT/'diagnosis01/finite_step_decomposition.mat',simplify_cells=True)
f=[];v=[]
for row in d['sweep']:
    i=np.flatnonzero(row['faceID']==238752)[0];f.append(row['fractionOfAcceptedStep']);v.append(row['nonlinearFluxRemainder'][i,0])
cross=json.loads((OUT/'mechanism.json').read_text())['crossingFractionOfAcceptedStep']
f.append(cross);v.append(0.);order=np.argsort(f);f=np.asarray(f)[order];v=np.asarray(v)[order]
axes[2].plot(f,v,'-o',color='#b65825',ms=4);axes[2].axvline(cross,ls='--',color='#626d75',label='右激波跨过采样位置')
axes[2].set(title='面 238752：质量通量的一阶偏差',xlabel='Q22→Q24 已接受增量的比例',ylabel='kg / (m²·s)');axes[2].legend(frameon=False);axes[2].grid(alpha=.2)
fig.suptitle('真实几何缩小复现：局部改善明显，但后续仍停滞；不代表整箭收敛',fontsize=14)
fig.savefig(OUT/'finite_step_comparison.png',dpi=160);fig.savefig(OUT/'finite_step_comparison.pdf');plt.close(fig)
