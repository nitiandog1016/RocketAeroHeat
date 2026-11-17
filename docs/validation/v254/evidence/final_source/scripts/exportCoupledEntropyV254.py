"""模块职责：积分候选稀疏映射回原面，同页入口曲线与可核对科研图；不复制旧大数据。"""
import json,subprocess,sys
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v254_validation';CLOSED=True


def export(stage):
    assert read(BASE/'audit01/report.json')['passed']
    assert stage=='export4','旧出口已封存，只允许整面筛查及分开拒绝原因的新出口。'
    out=BASE/'workbench04';out.mkdir(exist_ok=False)
    src=read(BASE/'prepare01/source.json');a=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'))
    old=read(ROOT/'outputs/v252_validation/matlab01/surface_base.json')
    screen=read(BASE/'surface_audit01/report.json');assert screen['passed']
    valid=np.isin(np.arange(1,19753),screen['retainedFaceIDs']);ids=np.flatnonzero(valid)
    # 旧anglesRad是31个沿程站，不是19752个面角；从原面质心与风轴直接计算。
    meta=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json')
    withheld=np.array(screen['withheldFaceIDs']);tri=a['originalVertices'][a['originalFaces'][withheld-1]]-np.array(meta['sphereCenter'])
    ctr=tri.mean(axis=1);cn=ctr/np.linalg.norm(ctr,axis=1)[:,None];vn=tri/np.linalg.norm(tri,axis=2)[:,:,None]
    theta=np.arccos(cn@(-np.array(meta['condition']['flowDirection'])))
    angularRadius=np.max(np.arctan2(np.linalg.norm(np.cross(vn,cn[:,None,:]),axis=2),np.einsum('ijk,ik->ij',vn,cn)),axis=1)
    startMask=theta-angularRadius<np.deg2rad(screen['withheldStartIntervalDeg'][1]);endMask=theta+angularRadius>np.pi/3
    assert np.all(startMask|endMask)
    reasons={str(int(fid)):('整面角界触及驻点起步未收敛区' if startMask[i] else '保守整面角界超出60°供体范围') for i,fid in enumerate(withheld)}
    radius=a['originalVertices'][a['originalFaces'][ids]].mean(axis=1)-np.array(meta['sphereCenter'])
    angles=np.arccos((radius@(-np.array(meta['condition']['flowDirection'])))/np.linalg.norm(radius,axis=1))
    stations=np.array(old['anglesRad']);lo=np.array(old['lowerStation'])[ids]-1;hi=np.array(old['upperStation'])[ids]-1
    recovered=stations[lo]+np.array(old['upperWeight'])[ids]*(stations[hi]-stations[lo])
    mappingAngleError=float(np.max(abs(angles-recovered)));assert mappingAngleError<1e-12
    assert len(ids)==screen['retainedFaces'] and angles.min()>=0 and angles.max()<np.pi/3
    assert len(a['originalFaces'])==19752 and abs(sum(a['originalFaceArea'])-1.9646292565715553)<2e-13
    labels={'base':'联合变熵 · 0.5°','fine':'联合变熵 · 0.25°','reference':'固定正激波熵 · 0.5°'}
    chains={k:read(BASE/'matlab02'/f'{k}.json') for k in labels};base=chains['base'];reference=chains['reference']
    method='CR-112180质量匹配＋CR-111922壁面性质积分＋CR-111921层流热相关式；隐式梯形'
    scope='拟合旋转球帽0—60°/Tw300K/小横流工程候选；实际三维动量、完整剖面及实验资格未闭合；β为压力梯度代理'
    fields=[]
    specs=[('entropyIntegralHeat','入壁热流 · 变熵积分','W/m²','heatFluxWm2',base),
           ('entropyFixedHeat','入壁热流 · 同公式固定熵','W/m²','heatFluxWm2',reference),
           ('entropyIntegralMassThickness','质量厚度 · 变熵积分','m','massThicknessM',base),
           ('entropyIntegralTe','外缘温度 · 变熵积分','K','edgeTemperatureK',base),
           ('entropyIntegralUe','外缘速度 · 变熵积分','m/s','edgeSpeedMps',base),
           ('entropyIntegralPressure','压力 · 旋转球帽近似','Pa','pressurePa',base),
           ('entropyIntegralBeta','压力梯度代理参数','1','betaPressureProxy',base)]
    for key,label,unit,name,c in specs:
        values=np.interp(angles,c['angleRad'],c[name]);fields.append(dict(key=key,label=label,unit=unit,values=values.tolist()))
    difference=np.array(fields[0]['values'])-fields[1]['values']
    fields.append(dict(key='entropyHeatDifference',label='热流差 · 联合变熵减固定熵',unit='W/m²',values=difference.tolist()))
    heat=np.array(fields[0]['values']);power=float(heat@a['originalFaceArea'][ids]);area=float(a['originalFaceArea'][ids].sum())
    ranges={k:[min(fields[0]['values']+fields[1]['values']),max(fields[0]['values']+fields[1]['values'])] for k in ('entropyIntegralHeat','entropyFixedHeat')}
    data=dict(schema='rah-coupled-integral-v254',originalFaces=19752,originalAreaM2=1.9646292565715553,faceIDs=(ids+1).tolist(),
        surfaceFields=fields,colorRanges=ranges,method=method,scope=scope,branches={},withheldFaceIDs=screen['withheldFaceIDs'],withheldReasons=reasons,
        caseNote=f"v254 变熵积分工程候选；原19752面及面积完整保留。{len(ids)}球帽面通过本轮数值筛查；{int(startMask.sum())}面触及未收敛起步区，{int(endMask.sum())}面保守角界超60°；其余区域未知。实际剖面/三维动量/实验资格未就绪，新合格面积0。",
        identity=dict(physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,actualProfileSolved=False,
            solvedSolidWallTemperatureK=None,wallTemperatureGivenK=300,candidateFaces=len(ids),rawCandidateFaces=1306,candidateAreaM2=area,signedCandidatePowerW=power,
            fullCapNumericalGatePassed=False,withheldStartIntervalDeg=screen['withheldStartIntervalDeg'],
            betaRangeExperimentallyValidated=False,sourceSHA256=src['sourceSHA256'],matlabSourceSHA256=sha(ROOT/'solvePhaseAwareEntropyIntegral.m'),
            actualMatlabTerminalKnown=True,oldV252ParentTerminalKnown=False))
    for branch,c in chains.items():
        records=[]
        for i,index in enumerate(c['targetIndices'][:46]):
            j=index-1
            row=dict(index=i,upstreamFaceID=int(a['upstreamFaceIDs'][i]),downstreamFaceID=int(a['downstreamFaceIDs'][i]),
                phiDeg=float(np.rad2deg(a['phiRad'][i])),thetaDeg=float(np.rad2deg(c['angleRad'][j])),wallXYZ=a['wallXYZ'][i].tolist(),
                actualMetricM=float(a['metricM'][i]),fittedMetricM=float(c['metricM'][j]*a['labelJacobian'][i]),
                givenWallTemperatureK=300,solvedSolidWallTemperatureK=None,actualProfile=None,profileReadyForCFD=False,physicalQualification=False,
                oldProfileAvailable=bool(np.isfinite(a['base_yM'][i,0])))
            for name in ('pressurePa','heatFluxWm2','edgeTemperatureK','edgeSpeedMps','massThicknessM','betaPressureProxy','alongMomentumDefectPam'):
                row[name]=c[name][j]
            records.append(row)
        data['branches'][branch]=dict(label=labels[branch],stations=records,sourceSHA256=sha(BASE/'matlab02'/f'{branch}.json'),
            along=dict(angleDeg=np.rad2deg(c['angleRad']).tolist(),heatFluxWm2=c['heatFluxWm2']))
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'coupled.json').write_text(payload+'\n',encoding='utf8');(out/'coupled-data.js').write_text('window.COUPLED_ENTROPY='+payload+';\n',encoding='utf8')
    html=(ROOT/'outputs/v252_validation/workbench02/index.html').read_text(encoding='utf8')
    html=html.replace('<p style="margin:10px 20px"><a href="profiles.html">查看沿程热流与实际边界层剖面 →</a></p>',
        '<details style="margin:12px 20px"><summary>展开本页：46处入口与变熵积分曲线</summary><iframe title="变熵积分曲线与来源" src="integral.html" style="width:100%;height:760px;border:1px solid #dbe1e7"></iframe></details>')
    html=html.replace('<script src="case-data.js"></script>',
        '<script src="../../v252_validation/workbench02/case-data.js"></script><script src="coupled-data.js"></script><script src="../../../web/coupled-entropy-v254.js"></script>')
    (out/'index.html').write_text(html,encoding='utf8')
    html=(ROOT/'web/coupled-entropy-v254.html').read_text(encoding='utf8').replace('src="coupled-entropy-v254.js"','src="../../../web/coupled-entropy-v254.js"')
    (out/'integral.html').write_text(html,encoding='utf8')
    assert read(out/'coupled.json')==data
    for f,(_,_,_,key,c) in zip(fields,specs):assert np.array_equal(np.array(f['values']),np.interp(angles,c['angleRad'],c[key]))
    plot(out,data,chains,a,valid)
    test=subprocess.run(['node',str(ROOT/'tests/testCoupledEntropyViewerV254.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout
    save(out/'report.json',dict(passed=True,originalFaces=19752,originalAreaM2=1.9646292565715553,candidateFaces=len(ids),candidateAreaM2=area,
        candidateAreaPercent=100*area/1.9646292565715553,signedCandidatePowerW=power,heatRangeWm2=[float(heat.min()),float(heat.max())],
        rawCandidateFaces=1306,withheldFaces=screen['withheldFaces'],withheldAreaM2=screen['withheldAreaM2'],fullCapNumericalGatePassed=False,
        startupWithheldFaces=int(startMask.sum()),outerBoundWithheldFaces=int(endMask.sum()),
        heatDifferenceRangeWm2=[float(difference.min()),float(difference.max())],fields=len(fields),inletBranches=3,totalIntegralPositions=138,
        actualProfilesFilled=0,physicalQualification=False,newQualifiedAreaM2=0,sharedWebGLChanged=False,browserInteractionTested=False,mappingAngleAbsoluteErrorRad=mappingAngleError,
        oldLargeDataCopied=False,overlaySHA256=sha(out/'coupled.json'),oldSurfaceDataSHA256=sha(ROOT/'outputs/v252_validation/workbench02/case-data.js')))


def plot(out,data,chains,a,valid):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for key,label in [('base','Coupled entropy, 0.5 deg'),('fine','Coupled entropy, 0.25 deg'),('reference','Fixed entropy, 0.5 deg')]:
        c=chains[key];angle=np.rad2deg(c['angleRad']);phi=[s['phiDeg'] for s in data['branches'][key]['stations']]
        axs[0,0].plot(angle,np.array(c['heatFluxWm2'])/1000,label=label)
        for ax,name,factor in [(axs[0,1],'heatFluxWm2',.001),(axs[1,0],'massThicknessM',1000),(axs[1,1],'edgeSpeedMps',1)]:
            ax.plot(phi,[s[name]*factor for s in data['branches'][key]['stations']],label=label)
    for ax,title,x,y in [(axs[0,0],'Fitted spherical cap','Polar angle [deg]','Heat into wall [kW/m2]'),
        (axs[0,1],'Actual 46-position ring','Body azimuth [deg]','Heat into wall [kW/m2]'),
        (axs[1,0],'Integral thickness; not a full profile','Body azimuth [deg]','Mass thickness [mm]'),
        (axs[1,1],'Entropy feedback on edge speed','Body azimuth [deg]','Edge speed [m/s]')]:
        ax.set(title=title,xlabel=x,ylabel=y);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axs[0,0].axvspan(*data['identity']['withheldStartIntervalDeg'],color='.8',alpha=.45,label='Unconverged start (retained raw)')
    axs[0,0].legend(fontsize=7)
    fig.suptitle('v254 | M5, alpha=20 deg, given Tw=300 K | engineering candidate, NOT CFD inlet')
    fig.savefig(out/'coupled_integral.png',dpi=160);plt.close(fig)
    V=a['originalVertices'];F=a['originalFaces'];q=np.array(data['surfaceFields'][0]['values'])/1000;norm=Normalize(q.min(),q.max())
    colors=np.tile([.76,.78,.8,1.],(len(F),1));colors[valid]=plt.get_cmap('viridis')(norm(q))
    fig=plt.figure(figsize=(12,5));fig.subplots_adjust(left=.07,right=.86,top=.84,bottom=.22,wspace=.35)
    for i,zoom in enumerate((False,True)):
        ax=fig.add_subplot(1,2,i+1);mask=np.min(V[F,0],axis=1)>.025 if zoom else np.ones(len(F),bool)
        tri=V[F[mask]];dims=[1,2] if zoom else [0,1];order=np.argsort(tri.mean(axis=1)[:,0 if zoom else 2])
        ax.add_collection(PolyCollection(tri[:,:,dims][order],facecolors=colors[mask][order],edgecolors='none',rasterized=True))
        ax.autoscale_view();ax.set_aspect('equal');ax.set(xlabel='Y [m]' if zoom else 'X [m]',ylabel='Z [m]' if zoom else 'Y [m]',title='Nose detail: YZ' if zoom else 'Original full rocket: XY');ax.grid(alpha=.15)
    cax=fig.add_axes([.91,.28,.015,.43]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax,label='Heat into wall [kW/m2]')
    fig.suptitle('v254 | original 19752 faces | alpha=20 deg | given Tw=300 K')
    fig.text(.06,.04,'Grey = unsolved / startup withheld. Coloured = screened integral candidate (%.4f%% area). Qualified new area = 0.'%(100*data['identity']['candidateAreaM2']/data['originalAreaM2']),fontsize=10)
    fig.savefig(out/'rocket_coupled_integral.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v254 fixed integral export closed; follow NEXT_DECISION')
    export(sys.argv[1])
