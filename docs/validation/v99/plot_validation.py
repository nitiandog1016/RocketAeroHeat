"""模块用途：绘制真实数值校核与已定位的轴线近似限制，不表示新增合格流场。"""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.io import loadmat
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
evidence=ROOT/'tests/data/shock_surface_mass/general_mass_evidence.mat';s=loadmat(evidence,simplify_cells=True)
axis=json.loads((HERE/'raw/axis_limit_audit.json').read_text(encoding='utf8'))
font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False})
fig,axes=plt.subplots(1,3,figsize=(12,3.8));m=s['general']['massCoordinateMap'];logEta=m['logMassRatio'];mask=logEta>=np.log(.01)
eta=np.exp(logEta[mask]);new=m['sourceCoordinateGridM'][:,mask];K=s['shape']['K'];q=m['coordinate']
old=np.column_stack([expm(K*x/np.trace(K))@q for x in logEta[mask]])
axes[0].plot(eta,np.linalg.norm(new-old,axis=0)*1e6,color='#176b85');axes[0].set_xlabel('同一激波线的质量比 eta');axes[0].set_ylabel('入口横向坐标差 / µm');axes[0].set_title('一般图与旧椭圆映射确有区别');axes[0].grid(color='#e2e7eb')
h=s['steps']*1000;error=abs(s['divergence']-1)
axes[1].loglog(h,error,'o-',color='#176b85');axes[1].set_xlabel('Cartesian差分步长 / mm');axes[1].set_ylabel('|div(a·grad(f)) − 1|');axes[1].set_title('投影连续式差分加密');axes[1].grid(color='#e2e7eb')
values=[c['spuriousAxisTangentialSpeedMS'] for c in axis['cases']];axes[2].bar([0,1],values,color=['#708b9c','#176b85'],width=.6)
axes[2].set_xticks([0,1]);axes[2].set_xticklabels(['定比热','NASA7连续积分']);axes[2].set_ylabel('显式式轴线延拓余速 / (m/s)');axes[2].set_title('现有压力闭合仍有驻点限制');axes[2].axhline(0,color='#a24c47',linestyle='--',label='轴线正则性要求切向速度为0');axes[2].set_ylim(0,80)
for j,v in enumerate(values):axes[2].text(j,v+2,'%.3f'%v,ha='center')
axes[2].legend(frameon=False,fontsize=8,loc='upper right');fig.suptitle('v99 一般共同激波质量映射 · 数值验证与适用性诊断',fontsize=12)
fig.text(.5,.015,'新投影质量映射通过数值复核；驻点整体、真实体形与热流仍未完成，新增合格面积0。',ha='center',color='#526879')
fig.tight_layout(rect=(0,.07,1,.94));fig.savefig(HERE/'general_mass_validation.png',dpi=170);fig.savefig(HERE/'general_mass_validation.svg');plt.close(fig)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(HERE/'plot_binding.json').write_bytes((json.dumps({'evidenceSHA256':sha(evidence),'axisAuditSHA256':sha(HERE/'raw/axis_limit_audit.json'),'plotProgramSHA256':sha(Path(__file__)),'newProductionSolvedAreaM2':0},indent=2)+'\n').encode('utf8'))
print(json.dumps({'status':'plotted-frozen-mass-validation-and-axis-limit-diagnostic'}))
