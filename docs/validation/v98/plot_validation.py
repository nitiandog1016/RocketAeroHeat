"""模块用途：从冻结MAT绘制可分享的层内剖面及真实体形残差，图不充当物理资格。"""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
BASE=ROOT/'tests/data/shock_layer_transport'
files=[BASE/name for name in ('parabolic_shock_layer_evidence.mat','conic_shock_layer_evidence.mat','translated_conic_layer_evidence.mat')]
p,c,t=[loadmat(file,simplify_cells=True) for file in files]
font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False})
fig,axes=plt.subplots(1,3,figsize=(12,3.9));colors=['#647787','#176b85']
for j,(label,style) in enumerate([('定比热','--'),('NASA7连续积分','-')]):
 col=p['cases'][j]['column'];n=col['normalDistanceM']*1000
 for ax,key,unit in [(axes[0],'pressure','kPa'),(axes[1],'temperature','K')]:
  factor=.001 if key=='pressure' else 1
  ax.plot(n,col['state'][key]*factor,style,color=colors[j],marker='o',markersize=2.6,label=label)
  ax.set_xlabel('激波内法向距离 n / mm');ax.set_ylabel('层内压力 / kPa' if key=='pressure' else '层内温度 / K');ax.grid(color='#e2e7eb',linewidth=.6);ax.legend(frameon=False,fontsize=8)
axes[0].set_title('同一给定共同激波的验证列');axes[1].set_title('同入口熵与总焓恢复状态')
values=np.array([max(abs(p['bodyFit']['initial']['residualM'])),p['report']['maximumBodyResidualM'],c['report']['maximumBodyResidualM'],t['report']['maximumFitResidualM'],t['report']['maximumHeldOutResidualM']])*1e6
axes[2].bar(np.arange(5),values,color=['#8c9aa5','#708b9c','#4e7d95','#176b85','#966b48'],width=.65)
axes[2].axhline(10,color='#a24c47',linestyle='--',linewidth=1,label='原匹配门槛10 µm');axes[2].set_yscale('log');axes[2].set_ylim(6,1000)
axes[2].set_xticks(np.arange(5));axes[2].set_xticklabels(['初始','4参数','5参数','7参数\n拟合点','7参数\n独立检查'])
for j,value in enumerate(values):axes[2].text(j,value*1.12,'%.1f'%value,ha='center',fontsize=8)
axes[2].set_ylabel('最大真实体形偏差 / µm');axes[2].set_title('拟合点改善不等于整个区域匹配');axes[2].legend(frameon=False,fontsize=8);axes[2].grid(axis='y',color='#e2e7eb',linewidth=.6);axes[2].set_axisbelow(True)
fig.suptitle('v98 内部薄激波层原型 · M5 / 30 km / 攻角20°',fontsize=12)
fig.text(.5,.015,'数值复算与工程原型验证；真实体形匹配未通过，热流、壁温及驻点整体未求。',ha='center',fontsize=9,color='#526879')
fig.tight_layout(rect=(0,.07,1,.94));fig.savefig(HERE/'shock_layer_validation.png',dpi=170);fig.savefig(HERE/'shock_layer_validation.svg');plt.close(fig)
sha=lambda file:hashlib.sha256(file.read_bytes()).hexdigest()
report={'status':'plotted-frozen-numerical-results-with-explicit-unqualified-scope','sourceFiles':[{'path':file.relative_to(ROOT).as_posix(),'sha256':sha(file)} for file in files],'plotProgramSHA256':sha(Path(__file__)),'maximumBodyResidualMicrometers':values.tolist(),'experimentalValidation':False,'newProductionSolvedAreaM2':0}
(HERE/'plot_binding.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
print(json.dumps({'status':report['status'],'maximumBodyResidualMicrometers':values.tolist()}))
