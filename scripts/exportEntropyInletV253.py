"""模块职责：真实交界环/剖面/匹配差量接同页研究入口；保留三维原面与所有缺失点。"""
from pathlib import Path
import json,sys,subprocess
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v253_validation';SRC=BASE/'inlet02'
CLOSED=True


def export():
    out=BASE/'workbench01';out.mkdir(exist_ok=False)
    audit=read(BASE/'audit02/report.json');assert audit['passed']
    rows=read(SRC/'inlets.json');a=np.load(str(SRC/'inlet_profiles.npz'));labels={'base':'2°剖面 · Tw300 K','fine':'1°剖面 · Tw300 K','history':'2°剖面 · 上游壁温变化'}
    data=dict(schema='rah-entropy-inlet-research-v253',physicalQualification=False,profileReadyForCFD=False,
        oldHeatOnly=True,originalFaces=19752,originalAreaM2=1.9646292565715553,sourceSHA256=sha(SRC/'inlets.json'),branches={})
    for branch,label in labels.items():
        stations=[]
        for raw in rows:
            if raw['branch']!=branch:continue
            s=dict(raw);i=s['index'];s.update(wallXYZ=a['wallXYZ'][i].tolist(),normal=a['unitNormal'][i].tolist(),tangent=a['unitTangent'][i].tolist(),
                metricM=float(a['metricM'][i]),labelJacobian=float(a['labelJacobian'][i]))
            s['profile']=None
            if s['available']:
                s['profile']=dict(y=a[branch+'_yM'][i].tolist(),temperature=a[branch+'_temperatureK'][i].tolist(),
                    speed=a[branch+'_speedMps'][i].tolist(),density=a[branch+'_densityKgm3'][i].tolist(),rhoE=a[branch+'_physicalConservative'][i,:,4].tolist())
            stations.append(s)
        data['branches'][branch]=dict(label=label,stations=stations)
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'inlet-data.js').write_text('window.ENTROPY_INLET='+payload+';\n',encoding='utf8');(out/'inlets.json').write_text(payload+'\n',encoding='utf8')
    html=(ROOT/'web/entropy-inlet-v253.html').read_text(encoding='utf8').replace('src="entropy-inlet-v253.js"','src="../../../web/entropy-inlet-v253.js"')
    (out/'inlet.html').write_text(html,encoding='utf8')
    old=ROOT/'outputs/v252_validation/workbench02/index.html';html=old.read_text(encoding='utf8')
    html=html.replace('src="case-data.js"','src="../../v252_validation/workbench02/case-data.js"').replace('href="profiles.html"','href="inlet.html"').replace('查看沿程热流与实际边界层剖面 →','查看真实前体入口与变熵相容性 →')
    html=html.replace('<main>','<p style="margin:10px 20px;color:#596775">本页三维热流沿用v252球鼻候选。本轮新结果在入口查看页，尚无新增热流面积。</p><main>')
    (out/'index.html').write_text(html,encoding='utf8')
    # 大型原面JSON/JS直接按只读路径复用，避免再次复制一整份结果。
    restored=read(out/'inlets.json');assert restored==data
    fig(out,a,rows)
    test=subprocess.run(['node',str(ROOT/'tests/testEntropyInletViewerV253.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout
    save(out/'report.json',dict(passed=True,branches=3,inletStations=138,availableProfiles=123,missingProfiles=15,
        fullNormalSamples=73923,sourceDataRoundtripExact=True,physicalQualification=False,newQualifiedAreaM2=0,
        sharedWebGLChanged=False,browserInteractionTested=False,existingSurfaceDataCopied=False,
        reusedSurfaceDataPath='outputs/v252_validation/workbench02/case-data.js',reusedSurfaceDataSHA256=sha(ROOT/'outputs/v252_validation/workbench02/case-data.js'),
        currentInletDataSHA256=sha(out/'inlets.json')))


def fig(out,a,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    phi=np.rad2deg(a['phiRad']);w=a['wallXYZ'];base=[r for r in rows if r['branch']=='base'];valid=np.array([r['available'] for r in base])
    fig,axs=plt.subplots(2,2,figsize=(11,7.3),layout='constrained')
    ax=axs[0,0];ax.plot(np.r_[w[:,1],w[0,1]],np.r_[w[:,2],w[0,2]],color='.8');ax.scatter(w[valid,1],w[valid,2],s=20,label='41 interpolated profiles');ax.scatter(w[~valid,1],w[~valid,2],s=30,color='.6',label='5 missing source profiles');ax.set_aspect('equal');ax.set(xlabel='Y [m]',ylabel='Z [m]',title='Actual 46-edge inlet ring');ax.legend(fontsize=9)
    ax=axs[0,1];ax.plot(phi,a['metricM']/np.linalg.norm(w[:,1:],axis=1));ax.axhline(1,color='.6',ls='--');ax.set(xlabel='Body azimuth [deg]',ylabel='Projected h / ring radius [-]',title='Cross-stream metric is not ring radius')
    for branch,label in [('base','2 deg, Tw300'),('fine','1 deg, Tw300'),('history','upstream Tw pulse')]:
        rr=[r for r in rows if r['branch']==branch]
        for ax,key,title in [(axs[1,0],'edgeTemperatureRelativeCorrection','Edge-temperature correction'),(axs[1,1],'edgeSpeedRelativeCorrection','Edge-speed correction')]:
            y=[100*r[key] if r['available'] else np.nan for r in rr];ax.plot(phi,y,label=label);ax.set(xlabel='Body azimuth [deg]',ylabel='Correction [%]',title=title);ax.legend(fontsize=8)
    for ax in axs.flat:ax.grid(alpha=.2)
    fig.suptitle('v253: source profiles versus mass-matched entropy | joint inlet NOT ready | no new qualified area')
    fig.savefig(out/'entropy_inlet.png',dpi=160);plt.close(fig)
    j=int(np.argmin(abs(phi)));fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for branch,label in [('base','2 deg, Tw300'),('fine','1 deg, Tw300'),('history','upstream Tw pulse')]:
        y=a[branch+'_yM'][j]*1000;axs[0].plot(a[branch+'_temperatureK'][j],y,label=label);axs[1].plot(a[branch+'_speedMps'][j],y,label=label)
    axs[0].set(xlabel='Fluid temperature [K]',ylabel='Wall-normal distance [mm]');axs[1].set(xlabel='Speed [m/s]',ylabel='Wall-normal distance [mm]')
    for ax in axs:ax.legend(fontsize=9);ax.grid(alpha=.2)
    fig.suptitle('Original interpolated profiles, body azimuth %.2f deg | not adjusted to matched entropy'%phi[j]);fig.savefig(out/'inlet_profiles.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v253 fixed inlet export closed; follow NEXT_DECISION')
    export()
