"""模块职责：从原生连续解目标直接接1306原面，保留旧失败对照及逐变量公式。"""
import json,subprocess
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v255_validation';CLOSED=False


def export():
    audit=read(BASE/'audit01/report.json');assert audit['passed']
    out=BASE/'workbench01';out.mkdir(exist_ok=False);src=BASE/audit['matlabFolder']
    source=read(BASE/'prepare01/source.json');geo=read(BASE/'prepare01/geometry.json')
    raw=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'));old=read(ROOT/'outputs/v254_validation/matlab02/base.json')
    labels={'base':'正则积分 · 主容差','fine':'正则积分 · 容差细化','seed':'正则积分 · 半种子控制'}
    chains={key:read(src/f'{key}.json') for key in labels};base=chains['base'];fine=chains['fine']
    ids=np.array(geo['coneAdmittedFaceIDs']);candidates=np.array(source['surfaceFaceIDs']);positions=np.flatnonzero(np.isin(candidates,ids));assert np.array_equal(candidates[positions],ids)
    target=np.array(source['surfaceTargetIndicesZeroBased'])[positions];j=np.array(base['targetIndices'])[target]-1;jf=np.array(fine['targetIndices'])[target]-1
    method='CR-112180质量匹配＋CR-111922层流积分＋CR-111921热相关式；正则变量与ode113'
    scope='拟合旋转球帽0—60°/Tw300K/小横流工程候选；三维动量/实际剖面/实验资格未闭合；驻点解析极限与有限种子区分'
    specs=[('entropyIntegralHeat','入壁热流 · 驻点正则积分','W/m²','heatFluxWm2',base,j,method),
        ('entropyRefinedHeat','入壁热流 · 固定容差细化','W/m²','heatFluxWm2',fine,jf,method),
        ('entropyIntegralMassThickness','质量厚度 · 正则积分','m','massThicknessM',base,j,'CR-111922 θm壁面性质积分＋CR-112180 δm=5.27θm；不是实际u/T剖面'),
        ('entropyIntegralTe','外缘温度 · 联合质量熵恢复','K','edgeTemperatureK',base,j,'同p及匹配激波熵；同相态Cp解析熵差反解；不是求解固体壁温'),
        ('entropyIntegralUe','外缘速度 · 同总焓恢复','m/s','edgeSpeedMps',base,j,'同总焓U=sqrt(2Δh)；相同Cp表稳定焓差'),
        ('entropyIntegralPressure','压力 · 修正牛顿球帽近似','Pa','pressurePa',base,j,'修正牛顿p=p∞+(p0−p∞)cos²θ；拟合旋转球帽工程近似'),
        ('entropyIntegralBeta','压力梯度代理参数','1','betaPressureProxy',base,j,'CR-111921式146；−p_s/(ρU)压力梯度代理；不是真实变熵速度导数')]
    fields=[]
    for key,label,unit,name,c,index,formula in specs:
        fields.append(dict(key=key,label=label,unit=unit,values=np.array(c[name])[index].tolist(),method=formula))
    heat=np.array(fields[0]['values']);angles=np.array(base['angleRad'])[j];oldQ=np.interp(angles,old['angleRad'],old['heatFluxWm2'])
    fields.append(dict(key='entropyHeatDifference',label='热流差 · 新解减旧v254诊断值',unit='W/m²',values=(heat-oldQ).tolist(),
        method='本轮原生连续解减旧v254线性映射诊断值；旧起步未收敛，不作为精确基准',computedStateLabel='版本差值诊断'))
    area=float(raw['originalFaceArea'][ids-1].sum());power=float(heat@raw['originalFaceArea'][ids-1])
    common=[min(fields[0]['values']+fields[1]['values']),max(fields[0]['values']+fields[1]['values'])]
    data=dict(schema='rah-regular-integral-v255',originalFaces=19752,originalAreaM2=1.9646292565715553,faceIDs=ids.tolist(),surfaceFields=fields,
        colorRanges={k:common for k in ['entropyIntegralHeat','entropyRefinedHeat']},method=method,scope=scope,branches={},
        withheldFaceIDs=geo['coneRejectedFaceIDs'],withheldReasons={str(i):'整面顶点超60°凸锥范围' for i in geo['coneRejectedFaceIDs']},
        caseNote='v255 驻点正则化已通过全采样数值对照，1306原球帽面为工程候选；其他18446面本方法未知。原19752面及完整分母保持。无实际剖面/三维或实验资格，新合格面积0。',
        identity=dict(physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,actualProfileSolved=False,solvedSolidWallTemperatureK=None,
            wallTemperatureGivenK=300,candidateFaces=len(ids),candidateAreaM2=area,signedCandidatePowerW=power,fullSampleNumericalGatesPassed=True,
            finiteSeedRad=1e-4,zeroAngleAnalytic=True,positiveBelowSeedSolved=False,recoveredPriorFaceIDs=audit['recoveredPreviouslyWithheldFaceIDs'],
            sourceSHA256=sha(BASE/'prepare01/source.json'),matlabSourceSHA256=sha(ROOT/'solveRegularizedEntropyIntegral.m'),actualMatlabTerminalKnown=True),
        legacyAlong=dict(angleDeg=np.rad2deg(old['angleRad']).tolist(),heatFluxWm2=old['heatFluxWm2']))
    for branch,c in chains.items():
        records=[]
        for i,index in enumerate(c['targetIndices'][:46]):
            k=index-1;row=dict(index=i,upstreamFaceID=int(raw['upstreamFaceIDs'][i]),downstreamFaceID=int(raw['downstreamFaceIDs'][i]),
                phiDeg=float(np.rad2deg(raw['phiRad'][i])),thetaDeg=float(np.rad2deg(c['angleRad'][k])),wallXYZ=raw['wallXYZ'][i].tolist(),
                actualMetricM=float(raw['metricM'][i]),fittedMetricM=float(c['metricM'][k]*raw['labelJacobian'][i]),givenWallTemperatureK=300,
                solvedSolidWallTemperatureK=None,actualProfile=None,profileReadyForCFD=False,physicalQualification=False,oldProfileAvailable=bool(np.isfinite(raw['base_yM'][i,0])))
            for key in ['pressurePa','heatFluxWm2','edgeTemperatureK','edgeSpeedMps','massThicknessM','betaPressureProxy','alongMomentumDefectPam']:row[key]=c[key][k]
            records.append(row)
        data['branches'][branch]=dict(label=labels[branch],stations=records,sourceSHA256=sha(src/f'{branch}.json'),
            along=dict(angleDeg=np.rad2deg(c['angleRad']).tolist(),heatFluxWm2=c['heatFluxWm2']))
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'coupled.json').write_text(payload+'\n',encoding='utf8');(out/'coupled-data.js').write_text('window.COUPLED_ENTROPY='+payload+';\n',encoding='utf8')
    html=(ROOT/'outputs/v252_validation/workbench02/index.html').read_text(encoding='utf8')
    html=html.replace('<p style="margin:10px 20px"><a href="profiles.html">查看沿程热流与实际边界层剖面 →</a></p>',
        '<details style="margin:12px 20px"><summary>展开本页：正则积分、46处入口与旧起步失败对照</summary><iframe title="正则积分曲线与来源" src="integral.html" style="width:100%;height:760px;border:1px solid #dbe1e7"></iframe></details>')
    html=html.replace('<script src="case-data.js"></script>','<script src="../../v252_validation/workbench02/case-data.js"></script><script src="coupled-data.js"></script><script src="../../../web/regularized-entropy-v255.js"></script>')
    (out/'index.html').write_text(html,encoding='utf8')
    html=(ROOT/'web/regularized-entropy-v255.html').read_text(encoding='utf8').replace('src="regularized-entropy-v255.js"','src="../../../web/regularized-entropy-v255.js"');(out/'integral.html').write_text(html,encoding='utf8')
    assert read(out/'coupled.json')==data
    plot(out,data,chains,raw,old)
    proc=subprocess.run(['node',str(ROOT/'tests/testRegularizedEntropyViewerV255.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(proc.stdout);assert proc.returncode==0,proc.stdout
    save(out/'report.json',dict(passed=True,originalFaces=19752,originalAreaM2=1.9646292565715553,candidateFaces=len(ids),candidateAreaM2=area,
        signedCandidatePowerW=power,heatRangeWm2=[float(heat.min()),float(heat.max())],restoredPreviousFaces=len(audit['recoveredPreviouslyWithheldFaceIDs']),
        methodsExplicitPerField=True,rawTargetMappingExact=True,fields=len(fields),inletRecords=138,physicalQualification=False,newQualifiedAreaM2=0,
        profileReadyForCFD=False,browserInteractionTested=False,sharedWebGLChanged=False,oldLargeDataCopied=False,
        oldDataSHA256=sha(ROOT/'outputs/v252_validation/workbench02/case-data.js'),dataSHA256=sha(out/'coupled.json')))


def plot(out,data,chains,raw,old):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for key,label in [('base','Regularized ode113'),('fine','Fixed tighter tolerance')]:
        r=chains[key];angle=np.rad2deg(r['angleRad']);q=np.array(r['heatFluxWm2'])/1000
        axs[0].plot(angle,q,label=label);near=angle<=5;axs[1].plot(angle[near],q[near],label=label)
        stations=data['branches'][key]['stations'];axs[2].plot([s['phiDeg'] for s in stations],[s['heatFluxWm2']/1000 for s in stations],label=label)
    t=np.rad2deg(old['angleRad']);near=t<=5;axs[1].plot(t[near],np.array(old['heatFluxWm2'])[near]/1000,'--',label='v254: unconverged startup')
    for ax,title,x in zip(axs,['Full sampled polar range','Startup: raw old defect retained','Actual 46-position ring'],['Polar angle [deg]','Polar angle [deg]','Body azimuth [deg]']):
        ax.set(title=title,xlabel=x,ylabel='Heat into wall [kW/m2]');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('v255 | M5, alpha=20 deg, given Tw=300 K | numerical comparison, NOT physical qualification')
    fig.savefig(out/'regularized_integral.png',dpi=160);plt.close(fig)
    V=raw['originalVertices'];F=raw['originalFaces'];ids=np.array(data['faceIDs'])-1;q=np.array(data['surfaceFields'][0]['values'])/1000;norm=Normalize(q.min(),q.max())
    colors=np.tile([.76,.78,.8,1.],(len(F),1));colors[ids]=plt.get_cmap('viridis')(norm(q))
    fig=plt.figure(figsize=(12,5));fig.subplots_adjust(left=.07,right=.86,top=.84,bottom=.22,wspace=.35)
    for i,zoom in enumerate((False,True)):
        ax=fig.add_subplot(1,2,i+1);mask=np.min(V[F,0],axis=1)>.025 if zoom else np.ones(len(F),bool);tri=V[F[mask]];dims=[1,2] if zoom else [0,1];order=np.argsort(tri.mean(axis=1)[:,0 if zoom else 2])
        ax.add_collection(PolyCollection(tri[:,:,dims][order],facecolors=colors[mask][order],edgecolors='none',rasterized=True));ax.autoscale_view();ax.set_aspect('equal')
        ax.set(xlabel='Y [m]' if zoom else 'X [m]',ylabel='Z [m]' if zoom else 'Y [m]',title='Nose detail: YZ' if zoom else 'Original full rocket: XY');ax.grid(alpha=.15)
    cax=fig.add_axes([.91,.28,.015,.43]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax,label='Heat into wall [kW/m2]')
    fig.suptitle('v255 | original 19752 faces | alpha=20 deg | given Tw=300 K')
    fig.text(.06,.04,'Grey = unsolved here. Coloured = 1306 spherical-cap engineering candidates (0.1791% area). New qualified area = 0.',fontsize=10)
    fig.savefig(out/'rocket_regularized_integral.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v255 fixed workbench export closed; follow NEXT_DECISION')
    export()
