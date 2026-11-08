"""模块：只读绘制压力—温度候选范围与参照偏差，不生成或暗示新的箭体流场。"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from standardDryAir import pressure_ceiling
ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'outputs/v204_validation/plot02';out.mkdir(exist_ok=False)
    r=json.loads((ROOT/'outputs/v204_validation/standard01/report.json').read_text(encoding='utf-8'))
    design=json.loads((ROOT/'outputs/v204_validation/domain_tests04/report.json').read_text(encoding='utf-8'))
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.4));t=np.linspace(100,2000,1000)
    a=axes[0];a.fill_between(t,.01,pressure_ceiling(t)/1000,color='#dde9f2');a.plot(t,pressure_ceiling(t)/1000,color='#286287',lw=2)
    a.scatter([226.50908361133,1231.35237218,design['isentropicDesign']['T0']],
        [1.1970279007383,100.292294882,design['isentropicDesign']['p0']/1000],color=['#276747','#b96e25','#962d4b'],s=30,zorder=3)
    a.set(xlim=(100,2000),ylim=(.01,1800),yscale='log',xlabel='Temperature [K]',ylabel='Pressure [kPa]',title='Named low-density engineering domain')
    a.text(.03,.04,'Dots: freestream / prior p with inferred T /\nisentropic design point (not a flow solution)',transform=a.transAxes,fontsize=8,color='#333333')
    a=axes[1];colors=['#286287','#b96e25','#276747','#962d4b']
    for key,label,color in zip(['rho','a','cp','cv'],['Density','Sound speed','Cp','Cv'],colors):
        ts=sorted({x['T'] for x in r['rows']});err=[100*max(abs(x['errors'][key]) for x in r['rows'] if x['T']==v) for v in ts]
        a.plot(ts,err,label=label,color=color,lw=1.6)
    a.axhline(1,color='#555555',ls='--',lw=.8,label='Density / sound gate: 1%');a.axhline(2,color='#777777',ls=':',lw=.8,label='Cp / Cv gate: 2%')
    a.set(xlim=(100,2000),ylim=(0,2.1),xlabel='Temperature [K]',ylabel='Maximum sampled absolute difference [%]',title='Ideal table vs full standard-Air EOS');a.legend(fontsize=8,loc='upper right')
    for a in axes:a.grid(True,alpha=.2);a.spines['top'].set_visible(False);a.spines['right'].set_visible(False)
    fig.suptitle('RocketAeroHeat v204 | Property-domain check; no new CFD or heat qualification',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(out/'pressure_domain.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
