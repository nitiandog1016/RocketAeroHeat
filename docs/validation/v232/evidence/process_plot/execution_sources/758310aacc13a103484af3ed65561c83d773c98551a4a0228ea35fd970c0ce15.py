"""模块职责：审核真实800更新保存状态与阶段身份，导出原面科研比较，不继续流动。"""
from pathlib import Path
import csv,json,math,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runPhaseEulerDevelopmentV232 import ROOT,BASE,DOC,CASE,BUILD,BASELINE,read,save,sha
import auditPhaseEulerRocketV230 as old

def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);run=read(BASE/'process_run/report.json')
    assert run['passed'] and run['inputsStable'] and not run['liveGroupAfterCleanup'] and run['stopReason'] is None
    end=read(CASE/'native_completion.json');assert end['normalReturn'] and end['updates']==800
    assert read(CASE/'prefix_verified.json')['first30NativeAndPhysicalLinesEqual']
    prefix={}
    for n in ('history.csv','physical_history.csv'):
        before=(BASELINE/'rocket01'/n).read_bytes();now=b''.join((CASE/n).read_bytes().splitlines(keepends=True)[:31])
        assert now==before,n;prefix[n]=sha(BASELINE/'rocket01'/n)
    for n in ('state_000000.csv','state_000001.csv'):assert sha(CASE/n)==sha(BASELINE/'rocket01'/n)
    steps=[int(p.stem.split('_')[1]) for p in CASE.glob('state_*.csv')];steps.sort();assert steps==[0,1]+list(range(50,801,50))
    h=np.loadtxt(str(CASE/'physical_history.csv'),delimiter=',',skiprows=1);assert h.shape==(800,22)
    assert np.array_equal(h[:,0],np.arange(1,801)) and np.all(h[:,1]==18027) and np.all(h[:,2:4]==0)
    r=np.loadtxt(str(CASE/'native_rejections.csv'),delimiter=',',skiprows=1);assert r.shape==(800,2)
    assert np.array_equal(r[:,0],h[:,0]) and np.all(r[:,1]==r[:,1].astype(int)) and np.all((r[:,1]>=0)&(r[:,1]<=18027))
    with (CASE/'history.csv').open(encoding='utf8') as f:nh=[{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(f,skipinitialspace=True)]
    assert len(nh)==800
    native4=np.array([[row[k] for k in ('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]')] for row in nh])
    errnative=float(np.max(abs(native4-np.log10(h[:,4:8]))));assert errnative<1e-8
    baseline=read(BASELINE/'audit02/report.json');scale=baseline['boundaryBalanceScale']
    balance=np.abs(h[:,[9,13]]-h[:,[14,18]])/np.array([scale['massKgs'],scale['energyW']]);assert float(balance.max())<1e-8
    # 壁压合力由更新后primitive读取，边界通量/RHS来自更新所用阶段：不可作为同阶段动量收支。
    cross_stage=np.linalg.norm(h[:,10:13]-h[:,15:18]-h[:,19:22],axis=1)/np.maximum(1,np.linalg.norm(h[:,19:22],axis=1))
    old.CASE=CASE;xyz=old.coords();g=np.load(str(BASELINE/'audit02/geometry.npz'))
    area=g['area'];ids=g['volumeNodeIDs'];tri=g['faces'];assert len(tri)==19752
    assert abs(float(area.sum())-1.9646292565715553)<1e-12
    selected=[50,400,800];saved=[];previous_p=None;last_pair=None;latest=None
    for step in steps[2:]:
        path=CASE/('state_%06d.csv'%step);field=old.fields(path,xyz);a=field['raw'];assert np.isfinite(a).all()
        physical=-a[:,10:15];physical[:,4]-=300000*physical[:,0]
        rms=np.sqrt(np.mean(physical**2,axis=0));rms_error=float(np.max(abs(rms-h[step-1,4:9])/np.maximum(1,abs(h[step-1,4:9]))));assert rms_error<2e-9
        p=field['p'][ids[tri]].mean(axis=1)
        change=None if previous_p is None else float(np.sum(area*abs(p-previous_p))/np.sum(area*abs(previous_p)))
        previous_p=p;last_pair=change
        sums=np.array([math.fsum(physical[:,j]) for j in range(5)])
        magnitude=np.sum(abs(physical),axis=0);magnitude[4]=np.sum(abs(a[:,14])+300000*abs(a[:,10]))
        eps=np.finfo(float).eps/2;ops=4*len(a)+8;bound=ops*eps/(1-ops*eps)*magnitude
        assert np.all(abs(sums-h[step-1,9:14])<=bound)
        stat=dict(update=step,sourceSHA256=sha(path),rmsRelativeError=rms_error,areaPressureChangeFromPrevious50=change,
          range={k:[float(field[k].min()),float(field[k].max())] for k in ('p','T','H','S')},allNodesJointDomainAccepted=True,
          sumWithinExplicitRoundoffBound=True,lastPreprocessingRejectedPoints=int(r[step-1,1]))
        saved.append(stat)
        if step in selected:
            data=np.column_stack([np.arange(1,19753),area]+[field[k][ids[tri]].mean(axis=1) for k in ('p','T','H','S')]+
              [field['v'][ids[tri],j].mean(axis=1) for j in range(3)])
            np.savetxt(str(out/('original_faces_%06d.csv'%step)),data,delimiter=',',comments='',fmt='%.17g',
              header='faceID,areaM2,pressurePa,fluidTemperatureK,totalEnthalpyJkg,entropyJkgK,Ux,Uy,Uz')
        if step==800:latest=field
    # 原始末场低温节点仅定位为候选；原生总数不能反推真实回退节点ID。
    order=np.argsort(latest['T'])[:12];low=[]
    for i in order:low.append(dict(pointID=int(i),xyz=xyz[i].tolist(),p=float(latest['p'][i]),T=float(latest['T'][i])))
    save(out/'low_temperature_candidates.json',dict(actualRejectedIDsKnown=False,notUniqueCauseEvidence=True,nodes=low))
    first=int(r[np.flatnonzero(r[:,1])[0],0]) if np.any(r[:,1]) else None
    save(out/'report.json',dict(passed=True,normalCompletion=True,completedUpdates=800,savedUpdates=steps,selectedUpdates=selected,
      prefix30ByteIdentical=True,prefixSHA256=prefix,initialAndUpdate1ByteIdentical=True,allSavedNewStatesJointDomainAccepted=True,
      savedFrameChecks=saved,historyNativeFirstFourLogError=errnative,massEnergyBoundaryAssemblyMaxScaled=float(balance.max()),
      boundaryBalanceScale=scale,steadyBoundaryNetScaled=(abs(h[-1,[14,18]])/np.array([scale['massKgs'],scale['energyW']])).tolist(),
      lastLogPhysicalRMS=np.log10(h[-1,4:9]).tolist(),allFiveBelowMinus8=bool(np.all(np.log10(h[-1,4:9])<=-8)),
      last50PressureAreaRelativeChange=last_pair,last50PressureStable=bool(last_pair<=.01),
      rejectionFirstUpdate=first,rejectionMaxLastPass=int(r[:,1].max()),rejectionLastPass=int(r[-1,1]),
      rejectionCountRole='last native preprocessing, not all substages or per-node identity',
      lastCrossStageMomentumDifference=float(cross_stage[-1]),momentumSameStageVerified=False,
      wallPressureForceStage='post-update primitive; differs from pre-update residual and characteristic boundary cache',
      spatialFarfieldVerified=False,scientificAcceptance=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
      originalFaces=19752,areaM2=float(area.sum()),actualRejectedNodeIDsKnown=False))
    print(json.dumps({k:read(out/'report.json')[k] for k in ('passed','completedUpdates','lastLogPhysicalRMS','last50PressureAreaRelativeChange','rejectionFirstUpdate','rejectionLastPass')}))

def export():
    out=BASE/'workbench01';out.mkdir(exist_ok=False);report=read(BASE/'audit01/report.json');assert report['passed']
    updates=report['selectedUpdates'];data={i:np.loadtxt(str(BASE/'audit01'/('original_faces_%06d.csv'%i)),delimiter=',',skiprows=1) for i in updates}
    model=read(BASELINE/'workbench01/rocket_euler_research.json');template=model['fields'][0];fields=[];ranges={};n=19752
    for key,label,unit,col,states in [('pressure','压力','Pa',2,updates),('fluidTemperature','Euler壁节点流体温度','K',3,updates),
      ('totalEnthalpy','总焓（hs298.15=0）','J/kg',4,[800]),('entropy','同模型比熵','J/(kg·K)',5,[800]),
      ('velocityX','速度 X','m/s',6,[800]),('velocityY','速度 Y','m/s',7,[800]),('velocityZ','速度 Z','m/s',8,[800])]:
        lo=min(float(data[s][:,col].min()) for s in states);hi=max(float(data[s][:,col].max()) for s in states)
        for step in states:
            f=dict(template);name=key+'_'+str(step)
            f.update(key=name,label=label+' · 更新%d'%step,unit=unit,values=data[step][:,col].tolist(),status=[1]*n,methodIndex=[0]*n,
              methodNames=['SU2 8.6 相态Air Euler/general Roe/implicit Euler＋原生域拒绝回退；更新后U反算，原壁三角顶点均值'],
              scope='M5/30km/α20°；800更新正常结束但未收敛；末预处理8点回退，非外缘/热解/固体壁温。',
              successfulUpdates=step,comparisonVariable=key,computedStateLabel='未收敛研究数值',availability='available',
              stateNote='更新编号不是物理时间；800步正常退出不授物理资格。')
            fields.append(f);ranges[name]=[lo,hi]
    fields.extend(model['fields'][-2:])
    model.update(fields=fields,comparisonColorRanges=ranges,defaultField='pressure_800',
      resultCondition=dict(caption='M5 · 30 km · α20° · 更新50/400/800 · 未收敛研究场'),
      caseNote='800更新正常完成；前30条与原计算相同。8个末预处理拒绝点不等于已知其面元位置。热流及固体壁温未知。',
      researchIdentity=dict(productionQualified=False,physicalQualification=False,model='phase-euler-native-primitive-rejection-v1',
        savedUpdates=report['savedUpdates'],selectedUpdates=updates,completedUpdates=800,normalCompletion=True,
        originalFaces=n,heatFluxSolved=False,lastPreprocessingRejectedPoints=8,
        sourceSHA256={str(i):sha(BASE/'audit01'/('original_faces_%06d.csv'%i)) for i in updates}))
    target=out/'rocket_euler_research.json';target.write_text(json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n',encoding='utf8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    vertices=np.asarray(model['vertices']);tri=vertices[np.array(model['faces'],int)]
    right=np.array([1.,0.,.1]);right/=np.linalg.norm(right);up=np.array([-.025,1.,.25]);up-=right*(up@right);up/=np.linalg.norm(up)
    xy=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(3,2,figsize=(13,8),facecolor='white')
    for c,(col,cmap,label,scale) in enumerate([(2,'viridis','Pressure [kPa]',.001),(3,'inferno','Euler wall-node fluid T [K]',1.)]):
        norm=Normalize(min(data[s][:,col].min() for s in updates)*scale,max(data[s][:,col].max() for s in updates)*scale)
        for row,step in enumerate(updates):
            ax=axs[row,c];ax.add_collection(PolyCollection(xy[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(data[step][order,col]*scale))))
            ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off();ax.margins(.03);ax.set_title('Saved update %d'%step,fontsize=10)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axs[:,c].tolist(),shrink=.75,label=label)
    fig.suptitle('RocketAeroHeat | original 19,752 faces | M5 / alpha 20 deg\n800 updates completed normally; not converged; common color limits',fontsize=12)
    fig.text(.08,.035,'Research only. Last native preprocessing rejected 8 points. Heat flux and solid temperature remain unknown.',fontsize=10)
    fig.subplots_adjust(top=.88,bottom=.10,left=.03,right=.80,hspace=.2,wspace=.1)
    fig.savefig(str(out/'rocket_euler_comparison.png'),dpi=160);plt.close(fig)
    result=subprocess.run(['node',str(ROOT/'tests/testPhaseEulerViewerV232.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout
    save(out/'report.json',dict(passed=True,fields=len(fields),originalFaces=n,selectedUpdates=updates,
      exportSHA256=sha(target),imageSHA256=sha(out/'rocket_euler_comparison.png'),dataRoundTripExact=True,
      browserInteractionTested=False,sharedWebChanged=False,newHeatAreaM2=0,newQualifiedAreaM2=0))
    print(json.dumps(read(out/'report.json')))

def main(stage):
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes;guard.size_bytes=lambda _:original(BASE)
    inputs=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/runPhaseEulerDevelopmentV232.py',
      ROOT/'scripts/auditPhaseEulerRocketV230.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
      ROOT/'scripts/su2PhaseAirFields.py',DOC/'PLAN.md',BASE/'process_run/report.json',BASELINE/'audit02/geometry.npz',BASELINE/'audit02/report.json']
    inputs += [p for p in CASE.iterdir() if p.is_file()]
    exporting=stage in ('export','export2','plot')
    if exporting:inputs += [ROOT/'web/viewer.js',ROOT/'tests/testPhaseEulerViewerV232.js',BASE/'audit01/report.json']
    # 全BASE容量包装已覆盖审核目录；不能再使用会二次累加audit目录的旧kind。
    seconds=120
    if exporting:seconds-=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind'] in ('audit','ui'))
    assert seconds>0
    task='task_plot' if stage=='plot' else 'task_export' if exporting else 'task_audit'
    return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),task],seconds,512,inputs,BASE,'ui' if exporting else 'offline')

def plot():
    # 仅重排科学图；固定三阶段原面数据与JSON，不重算物性或覆盖首图证据。
    out=BASE/'plot02';out.mkdir(exist_ok=False);updates=[50,400,800]
    data={i:np.loadtxt(str(BASE/'audit01'/('original_faces_%06d.csv'%i)),delimiter=',',skiprows=1) for i in updates}
    g=np.load(str(BASELINE/'audit02/geometry.npz'));tri=g['vertices'][g['faces']]
    right=np.array([1.,0.,.1]);right/=np.linalg.norm(right);up=np.array([-.025,1.,.25]);up-=right*(up@right);up/=np.linalg.norm(up)
    xy=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize,LogNorm
    fig,axs=plt.subplots(3,2,figsize=(13,7),facecolor='white')
    fig.subplots_adjust(top=.85,bottom=.22,left=.03,right=.97,hspace=.20,wspace=.12)
    ranges={}
    for c,(col,cmap,label,scale) in enumerate([(2,'viridis','Pressure [kPa] (log scale)',.001),(3,'inferno','Euler wall-node fluid T [K]',1.)]):
        lo=min(data[s][:,col].min() for s in updates)*scale;hi=max(data[s][:,col].max() for s in updates)*scale
        assert lo>0;norm=LogNorm(lo,hi) if c==0 else Normalize(lo,hi);ranges[str(col)]=[float(lo),float(hi)]
        for row,step in enumerate(updates):
            ax=axs[row,c];ax.add_collection(PolyCollection(xy[order],linewidths=0,edgecolors='none',antialiaseds=False,
              facecolors=plt.get_cmap(cmap)(norm(data[step][order,col]*scale))))
            ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off();ax.margins(.03);ax.set_title('Saved update %d'%step,fontsize=10)
        cax=fig.add_axes([.07+c*.49,.14,.36,.02]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',label=label)
    fig.suptitle('RocketAeroHeat | original 19,752 faces | M5 / alpha 20 deg\n800 updates completed normally; not converged; common color limits',fontsize=12)
    fig.text(.05,.035,'Research only. Last native preprocessing rejected 8 points. Heat flux and solid temperature remain unknown.',fontsize=10)
    target=out/'rocket_euler_comparison.png';fig.savefig(str(target),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,imageSHA256=sha(target),dataUnchanged=True,pressureScale='logarithmic',
      temperatureScale='linear',commonRanges=ranges,sourceSHA256={str(i):sha(BASE/'audit01'/('original_faces_%06d.csv'%i)) for i in updates},
      JSONSHA256=sha(BASE/'workbench01/rocket_euler_research.json'),firstPlotPreserved=True,browserInteractionTested=False))

if __name__=='__main__':
    if sys.argv[1]=='task_audit':audit()
    elif sys.argv[1]=='task_export':export()
    elif sys.argv[1]=='task_plot':plot()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
