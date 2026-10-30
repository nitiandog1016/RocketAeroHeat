"""用实际保存结果生成可导出的定量图：驻点改善、守恒加密、稳定性失败及真实体形残差。"""

from pathlib import Path
import argparse
import json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from auditAxisymmetricShockLayerPatch import plain


def plot(evidence,audit_file,folder):
    """模块用途：所有曲线来自原结果或已说明旧解析压力式，不把未解热流补成曲线。"""
    folder.mkdir(parents=True,exist_ok=True)
    s=plain(loadmat(evidence,simplify_cells=True));a=json.loads(audit_file.read_text(encoding='utf-8'))
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.8));colors=['#226b83','#a65c30'];labels=['Constant cp','NASA7']
    rho=float(s['upstream']['density']);U=float(np.linalg.norm(s['upstream']['velocityVector']))
    for mode in range(2):
        p=s['cases'][mode,-1];eta=np.asarray(p['eta']);shock=s['limits'][mode]['shock']
        normal=float(np.dot(shock['velocityVector'],shock['normal']))
        old=float(shock['pressure'])+.5*rho*U*normal*(1-eta*eta)
        axes[0,0].plot(eta,p['pressurePa'][:,0]/1000,color=colors[mode],label=labels[mode]+' coupled')
        axes[0,0].plot(eta,old/1000,color=colors[mode],ls='--',label=labels[mode]+' old explicit')
        rows=[r for r in a['primaryPatches'] if r['gasMode']==p['gas']['thermalPropertyMode']]
        axes[0,1].semilogy([r['radialNodes'] for r in rows],[r['physicalResiduals']['continuity'] for r in rows],'-o',color=colors[mode],label=labels[mode])
        rows=[r for r in a['independentIntegrations'] if r['gasMode']==p['gas']['thermalPropertyMode']]
        axes[1,0].semilogy([r['radialNodes'] for r in rows],
            [max(r['result']['normalizedRegularSpeedDifference'],r['result']['normalizedNormalSpeedDifference']) for r in rows],'-o',color=colors[mode],label=labels[mode])
    ax=axes[0,0];ax.set(xlim=(1,0),xlabel=r'Mass ratio $\eta$ (shock 1 to generated body 0)',ylabel='Axis pressure / kPa',title='(a) Axis closure, same upstream H and source entropy');ax.legend(fontsize=8,ncol=2)
    ax=axes[0,1];ax.axhline(1e-5,color='#333333',ls=':',label='Numerical gate');ax.set(xlabel='Radial collocation nodes',ylabel='Direct physical continuity residual / normalizer',title='(b) Independent physical-coordinate differentiation');ax.legend(fontsize=8)
    ax=axes[1,0];ax.axhline(2e-7,color='#333333',ls=':',label='Fixed comparison gate');ax.set(xlabel='Radial collocation nodes',ylabel='max(error(w)/(U k), error(v)/U)',title='(c) DOP853 vs ode113: high-order stability not passed');ax.set_xticks([13,19]);ax.legend(fontsize=8)
    b=s['realBody'];r=np.linalg.norm(b['coordinates'],axis=1).reshape(3,3)*1000;res=np.asarray(b['residualM']).reshape(3,3)*1e6
    ax=axes[1,1]
    for j,color in enumerate(['#226b83','#a65c30','#59625e']):ax.plot(r[:,j],res[:,j],'-o',color=color,label=r'$\phi$ = %d deg'%(120*j))
    ax.axhline(10,color='#333333',ls=':');ax.axhline(-10,color='#333333',ls=':',label='Original +/-10 um gate')
    ax.set(xlabel='Shock radial coordinate / mm',ylabel='Predicted minus real surface distance / um',title='(d) Fixed candidate on original 19,752 triangles; no fit');ax.legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.22);ax.tick_params(labelsize=9)
    fig.suptitle('M5 / 30 km: coupled axisymmetric shock-layer research',fontsize=13)
    fig.text(.5,.015,'Original rocket alpha = 20 deg; flow axisymmetry is an additional assumption. No new qualified heat-transfer area.',ha='center',fontsize=9,color='#495966')
    fig.tight_layout(rect=(0,.045,1,.96));fig.savefig(folder/'coupled_patch_quantitative.png',dpi=180);fig.savefig(folder/'coupled_patch_quantitative.svg');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('audit',type=Path);p.add_argument('folder',type=Path)
    args=p.parse_args();plot(args.evidence,args.audit,args.folder)
