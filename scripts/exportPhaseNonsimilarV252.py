"""模块职责：三分支非相似剖面接原面WebGL，输出有单位的沿程/法向曲线。"""
import json,re,sys,subprocess
from pathlib import Path
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v252_validation';SOURCE=BASE/'matlab01'
BRANCHES={'base':'2°沿程 · Tw300 K','fine':'1°沿程 · Tw300 K','history':'2°沿程 · 上游壁温变化'}
CLOSED=True


def export(stage='export'):
    out=BASE/('workbench02' if stage=='export2' else 'workbench01');out.mkdir(exist_ok=False)
    audit=read(BASE/'audit01/report.json');assert audit['allNumericalGatesPassed']
    model=read(SOURCE/'geometry_viewer.json');chains={k:read(SOURCE/(k+'.json')) for k in BRANCHES}
    fields=model['fields'];ranges={};newkeys=[];profiledata={}
    next(f for f in fields if f['key']=='wallInput')['label']='给定壁温 · 基准工况300 K'
    specs=[('heatFlux','入壁对流热流','W/m²'),('wallShearStress','壁面剪应力','Pa'),('delta','有限边界层厚度','m'),
           ('deltaStar','压缩性位移厚度','m'),('thetaPhysical','压缩性动量厚度','m'),
           ('massThickness','质量厚度','m'),('energyThickness','机械能厚度','m')]
    surfaces={k:read(SOURCE/('surface_'+k+'.json')) for k in BRANCHES}
    for key,label,unit in specs:
        values=[v for h in surfaces.values() for v in h[key] if v is not None]
        common=[min(values),max(values)]
        for branch,caption in BRANCHES.items():
            h=surfaces[branch];name=key+'_'+branch;newkeys.append(name);ranges[name]=common
            fields.append(dict(key=name,label=label+' · '+caption,unit=unit,values=h[key],
                status=[int(v) for v in h['validFace']],methodIndex=[0]*19752,methodNames=[h['method']],scope=h['scope']+'；本版β在0—1范围',
                missingReason=h['reason'],computedStateLabel='未验证非相似工程候选',availability='available',
                stateNote='上游剖面沿程传递；输入为旋转球鼻近似，未授真实三维/实验资格。不是物理时间。',
                comparisonVariable=key,branch=branch))
    # 各壁温历史显式给定；不能把中段变化隐去，仅显示终点300 K。
    for branch,caption in BRANCHES.items():
        h=surfaces[branch]
        fields.append(dict(key='wallInput_'+branch,label='给定壁温 · '+caption,unit='K',values=h['wallTemperature'],
            status=[2 if v else 0 for v in h['validFace']],methodIndex=[0]*19752,methodNames=['输入壁温历史；不是求解固体温度'],
            scope=h['scope'],missingReason=h['reason'],availability='available'))
    hb=surfaces['base'];hp=surfaces['history'];difference=[b-a if a is not None and b is not None else None for a,b in zip(hb['heatFlux'],hp['heatFlux'])]
    fields.append(dict(key='historyHeatDifference',label='上游壁温变化造成的热流差',unit='W/m²',values=difference,
        status=[int(v) for v in hb['validFace']],methodIndex=[0]*19752,methodNames=['上游壁温变化分支减全程300 K分支；同一非相似方法'],
        scope=hb['scope'],missingReason=hb['reason'],computedStateLabel='候选分支差值',availability='available'))
    model.update(defaultField='heatFlux_base',comparisonColorRanges=ranges,
        caseNote='v252 同相态Levy–Lees非相似层流。全原面仍在，仅1306球帽面为研究候选；鼻外/强三维区域未解。新合格面积0。MATLAB完整结果存在，原父监督终态缺失；独立复算与积分检查通过。',
        resultCondition=dict(caption='M5 · α20° · 原19752面 · 旋转球帽工程近似 · 3种沿程/壁温分支 · 非物理时间'),
        researchIdentity=dict(physicalQualification=False,productionQualified=False,newQualifiedAreaM2=0,
            scope='rotated-spherical-research',heatSummary=hb['summary'],branchSummaries={k:h['summary'] for k,h in surfaces.items()},
            chainSHA256={k:c['sourceSHA256'] for k,c in chains.items()},auditSHA256=sha(BASE/'audit01/report.json'),
            originalMatlabSupervisorTerminalKnown=False))
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'rocket_workbench.json').write_text(payload+'\n',encoding='utf8')
    (out/'case-data.js').write_text('window.ROCKET_CASE='+payload+';\n',encoding='utf8')
    html=(ROOT/'web/index.html').read_text(encoding='utf8')
    html=re.sub(r'((?:src|href)=")((?:viewer|section|export)\.(?:js|css)(?:\?[^\"]*)?)(")',r'\1../../../web/\2\3',html)
    html=html.replace('<main>','<p style="margin:10px 20px"><a href="profiles.html">查看沿程热流与实际边界层剖面 →</a></p><main>')
    (out/'index.html').write_text(html,encoding='utf8')
    for branch,c in chains.items():
        stations=[]
        for r in c['profiles']:
            # 页面降采样只用于绘图；完整601点/自适应解仍在原始结果中保存。
            p=r['profile'];idx=np.linspace(0,len(p['eta'])-1,101,dtype=int)
            stations.append(dict(angleDeg=float(np.rad2deg(r['station']['thetaRad'])),arcM=r['station']['arcM'],
                heatFlux=r['heatFlux'],wallShearStress=r['wallShearStress'],wallTemperature=r['wallTemperature'],
                edgeTemperature=r['edge']['temperatureK'],edgeSpeed=r['edge']['speedMps'],
                y=[p['physicalY'][i] for i in idx],temperature=[p['temperature'][i] for i in idx],
                speedRatio=[p['states'][1][i] for i in idx],sourceChainSHA256=c['sourceSHA256']))
        profiledata[branch]=dict(label=BRANCHES[branch],stations=stations,method=c['method'],scope=c['scope'])
    save(out/'profiles.json',profiledata)
    (out/'profiles-data.js').write_text('window.PHASE_PROFILES='+json.dumps(profiledata,ensure_ascii=False,separators=(',',':'))+';\n',encoding='utf8')
    html=(ROOT/'web/phase-profiles-v252.html').read_text(encoding='utf8').replace('src="phase-profiles-v252.js"','src="../../../web/phase-profiles-v252.js"')
    (out/'profiles.html').write_text(html,encoding='utf8')
    returned=read(out/'rocket_workbench.json');assert returned['vertices']==model['vertices'] and returned['faces']==model['faces']
    for f in returned['fields']:assert f==next(g for g in model['fields'] if f['key']==g['key'])
    plot(out,chains,model,surfaces)
    test=subprocess.run(['node',str(ROOT/'tests/testPhaseNonsimilarViewerV252.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
    (out/'frontend_test.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout
    save(out/'report.json',dict(passed=True,fields=len(fields),originalFaces=19752,newFields=newkeys,
        exactSurfaceValuesChecked=sum(len(f['values']) for f in fields),candidateFaces=1306,physicalQualification=False,
        newQualifiedAreaM2=0,sharedWebGLChanged=False,browserInteractionTested=False,
        profileBranchCounts={k:len(v['stations']) for k,v in profiledata.items()},
        frontendSHA256=sha(ROOT/'web/viewer.js'),dataSHA256=sha(out/'rocket_workbench.json')))


def plot(out,chains,model,surfaces):
    """模块四：科研图使用固定单位/共同色标，未解面灰色，原真实尺度保持。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    labels={'base':'2 deg, Tw=300 K','fine':'1 deg, Tw=300 K','history':'2 deg, upstream wall-T pulse'}
    for branch,c in chains.items():
        deg=np.rad2deg([r['station']['thetaRad'] for r in c['profiles']]);q=np.array([r['heatFlux'] for r in c['profiles']])/1e6
        axs[0,0].plot(deg,q,label=labels[branch]);axs[0,1].plot(deg,c['wallTemperatureHistory'],label=labels[branch])
        j=np.argmin(abs(deg-40));p=c['profiles'][j]['profile'];y=np.array(p['physicalY'])*1e3
        axs[1,0].plot(np.array(p['states'])[1],y,label=labels[branch]);axs[1,1].plot(p['temperature'],y,label=labels[branch])
    for ax,title,xlabel,ylabel in [(axs[0,0],'Streamwise heating','Surface polar angle [deg]','Heat into wall [MW/m2]'),
        (axs[0,1],'Prescribed wall-temperature history','Surface polar angle [deg]','Given Tw [K]'),
        (axs[1,0],'Velocity profile at 40 deg','u / Ue [-]','Distance from wall [mm]'),
        (axs[1,1],'Temperature profile at 40 deg','Fluid temperature [K]','Distance from wall [mm]')]:
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('v252: non-similar laminar engineering candidate | no experimental / 3D qualification')
    fig.savefig(out/'profiles_and_history.png',dpi=160);plt.close(fig)
    V=np.array(model['vertices']);F=np.array(model['faces']);valid=np.array(surfaces['base']['validFace'],bool)
    q=np.array(surfaces['base']['heatFlux'],float)/1e6;norm=Normalize(np.nanmin(q),np.nanmax(q))
    colors=np.tile([.76,.78,.8,1.],(len(F),1));colors[valid]=plt.get_cmap('viridis')(norm(q[valid]))
    fig=plt.figure(figsize=(12,5));fig.subplots_adjust(left=.07,right=.86,top=.84,bottom=.22,wspace=.35)
    for j,zoom in enumerate((False,True)):
        ax=fig.add_subplot(1,2,j+1);mask=np.min(V[F,0],axis=1)>.025 if zoom else np.ones(len(F),bool)
        tri=V[F[mask]];dims=[1,2] if zoom else [0,1];depth=0 if zoom else 2;order=np.argsort(tri.mean(axis=1)[:,depth])
        proj=tri[:,:,dims];ax.add_collection(PolyCollection(proj[order],facecolors=colors[mask][order],edgecolors='none',rasterized=True))
        ax.autoscale_view();ax.set_aspect('equal');ax.set(xlabel='Y [m]' if zoom else 'X [m]',ylabel='Z [m]' if zoom else 'Y [m]',title='Original whole rocket: XY view' if not zoom else 'Nose detail: YZ view')
        ax.tick_params(labelsize=9);ax.grid(alpha=.15)
    cax=fig.add_axes([.90,.28,.015,.43]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax,label='Heat into wall [MW/m2]')
    fig.suptitle('v252 | original 19752 faces | alpha=20 deg | prescribed Tw=300 K')
    fig.text(.06,.04,'Grey = unsolved here. Coloured = spherical-cap engineering candidate (0.1791% area). Qualified new area = 0.',fontsize=10)
    fig.savefig(out/'rocket_nonsimilar.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v252 nonsimilar export closed; follow NEXT_DECISION')
    export(sys.argv[1])
