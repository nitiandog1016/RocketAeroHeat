"""模块职责：只读审核失败候选、实际保存场及原面科研出口，不重启流动。"""
from pathlib import Path
import csv,json,math,subprocess,sys
import numpy as np
from scipy.spatial import cKDTree
import runSU2RocketPhaseV222 as guard
from runPhaseEulerRocketV230 import ROOT,BASE,CASE,BUILD,DOC,read,save,sha
from su2PhaseAirFields import inverse_energy
from phaseAwareViscousAir import PhaseAwareViscousAir

def inspect():
    out=BASE/'crash01';out.mkdir(exist_ok=False)
    tool=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/objdump.exe'
    argv=[str(tool),'-d','--start-address=0x140b1a9f0','--stop-address=0x140b1aa10',str(BUILD/'rocket.exe')]
    r=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15);assert r.returncode==0
    (out/'fault_disassembly.txt').write_bytes(r.stdout)
    text=r.stdout.decode('utf8');assert '__trap_terminate_handler' in text and 'ud2' in text
    report=Path('C:/ProgramData/Microsoft/Windows/WER/ReportArchive/AppCrash_rocket.exe_821494c91f5feaa518c38950b0a1035633a7fe_dd86ac11_0916db5c-3d25-47ab-acc4-8970e8483d02/Report.wer')
    data=report.read_text(encoding='utf-16');assert 'rocket.exe' in data and 'c000001d' in data.lower() and '0000000000b1aa00' in data.lower()
    # 仅发布本程序故障签名，系统报告其他机器信息不复制到仓库。
    save(out/'report.json',dict(passed=True,executableSHA256=sha(BUILD/'rocket.exe'),werSourceSHA256=sha(report),
      windowsExceptionCode='0xc000001d',imageBase='0x140000000',faultModuleOffset='0xb1aa00',
      symbol='__gnu_cxx::__trap_terminate_handler()',instruction='ud2',
      terminateEstablished=True,originalExceptionTypeKnown=False,underlyingPhysicalOrSoftwareCauseKnown=False,
      notEvidenceOfUnsupportedCPUInstruction=True,CFDReplays=0,disassemblyCommand=argv))
    print(json.dumps(dict(passed=True,terminateEstablished=True,causeKnown=False,replays=0)))

def coords():
    lines=(CASE/'mesh.su2').read_text(encoding='utf8').splitlines();j=next(i for i,s in enumerate(lines) if s.startswith('NPOIN='))
    count=int(lines[j].split('=')[1]);v=np.empty((count,3))
    for line in lines[j+1:j+1+count]:
        w=line.split();v[int(w[3])]=list(map(float,w[:3]))
    return v

def fields(path,xyz):
    with path.open(encoding='utf8') as f:header=next(csv.reader(f))
    assert header==['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy','Volume']+['RHS%d'%i for i in range(5)]
    a=np.loadtxt(str(path),delimiter=',',skiprows=1);assert a.shape==(18027,15)
    ids=a[:,0].astype(int);assert np.array_equal(np.sort(ids),np.arange(18027)) and np.array_equal(ids,a[:,0])
    assert np.array_equal(a[:,1:4],xyz[ids]);a=a[np.argsort(ids)];assert np.isfinite(a[:,:10]).all()
    rho=a[:,4];assert np.all(rho>0) and np.all(a[:,9]>0)
    v=a[:,5:8]/rho[:,None];e=a[:,8]/rho-.5*np.sum(v*v,axis=1)-300000
    air=PhaseAwareViscousAir();T=inverse_energy(air,e);p=rho*air.R*T;q=air.evaluate(p,T)
    H=np.asarray(q['enthalpyJkg'])+.5*np.sum(v*v,axis=1);S=np.asarray(q['entropyJkgK'])
    return dict(raw=a,p=p,T=T,H=H,S=S,v=v)

def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);run=read(BASE/'process_run/report.json')
    assert not run['passed'] and run['exitCode']==3221225501 and run['stopReason'] is None
    assert run['inputsStable'] and not run['liveGroupAfterCleanup']
    assert read(BASE/'crash01/report.json')['terminateEstablished']
    xyz=coords();frames={i:fields(CASE/('state_%06d.csv'%i),xyz) for i in (0,1)}
    assert np.isnan(frames[0]['raw'][:,10:]).all() and np.isfinite(frames[1]['raw']).all()
    air=PhaseAwareViscousAir();initial=air.evaluate(1197.0279007383137,226.50908361133003)
    v=5*float(initial['soundSpeedMps'])*np.array([-math.cos(math.radians(20)),0,math.sin(math.radians(20))])
    rho=float(initial['densityKgm3']);U=np.r_[rho,rho*v,rho*(float(initial['internalEnergyJkg'])+300000+.5*(v@v))]
    initial_error=float(np.max(abs(frames[0]['raw'][:,4:9]-U)/np.maximum(1,abs(U))));assert initial_error<2e-9
    h=np.loadtxt(str(CASE/'physical_history.csv'),delimiter=',',skiprows=1);assert h.shape==(30,22)
    assert np.array_equal(h[:,0],np.arange(1,31)) and np.all(h[:,1]==18027) and np.all(h[:,2:4]==0)
    physical=-frames[1]['raw'][:,10:15];physical[:,4]-=300000*physical[:,0]
    rms=np.sqrt(np.mean(physical**2,axis=0));sums=physical.sum(axis=0)
    err=lambda a,b:float(np.max(abs(a-b)/np.maximum(1,abs(b))))
    rms_error=err(rms,h[0,4:9]);sum_error=err(sums,h[0,9:14]);assert rms_error<2e-9 and sum_error<2e-9
    balance=np.abs(h[:,[9,13]]-h[:,[14,18]])/np.maximum(1.,np.maximum(abs(h[:,[9,13]]),abs(h[:,[14,18]])))
    assert float(balance.max())<1e-8
    with (CASE/'history.csv').open(encoding='utf8') as f:native=[{k.strip():float(v) for k,v in r.items()} for r in csv.DictReader(f,skipinitialspace=True)]
    assert len(native)==30
    native4=np.array([[r[k] for k in ('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]')] for r in native])
    native_error=float(np.max(abs(native4-np.log10(h[:,4:8]))));assert native_error<1e-8
    # 全域动量差额仅是推定约束反力，与压力合力独立比较，不能凭恒等式授资格。
    reaction=h[:,10:13]-h[:,15:18];wall=h[:,19:22]
    reaction_error=np.linalg.norm(reaction-wall,axis=1)/np.maximum(1,np.linalg.norm(wall,axis=1))
    surf=read(ROOT/'outputs/v194_validation/rocket01/surface.json');vertices=np.array(surf['vertices']);tri=np.array(surf['facesOneBased'],dtype=int)-1
    distance,ids=cKDTree(xyz).query(vertices);assert np.all(distance==0) and len(set(ids.tolist()))==len(ids)
    assert np.array_equal(surf['faceIDs'],np.arange(1,19753))
    area=np.linalg.norm(np.cross(vertices[tri[:,1]]-vertices[tri[:,0]],vertices[tri[:,2]]-vertices[tri[:,0]]),axis=1)/2
    assert abs(float(sum(area))-1.9646292565715553)<1e-12
    stats={}
    for step,f in frames.items():
        data=np.column_stack([np.arange(1,19753),area]+[f[k][ids[tri]].mean(axis=1) for k in ('p','T','H','S')]+[f['v'][ids[tri],j].mean(axis=1) for j in range(3)])
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),data,delimiter=',',header='faceID,areaM2,pressurePa,fluidTemperatureK,totalEnthalpyJkg,entropyJkgK,Ux,Uy,Uz',comments='',fmt='%.17g')
        stats[str(step)]={k:[float(f[k].min()),float(f[k].max())] for k in ('p','T','H','S')}
    np.savez_compressed(str(out/'geometry.npz'),vertices=vertices,faces=tri,area=area,volumeNodeIDs=ids)
    save(out/'report.json',dict(passed=True,scientificAcceptance=False,normalCompletion=False,completedUpdates=30,
      savedFullStateUpdates=[0,1],latestCompleteFrameIsNotUpdate30=True,initialError=initial_error,
      firstPhysicalRMSRelativeError=rms_error,firstPhysicalSumScaledError=sum_error,historyNativeFirstFourLogError=native_error,
      massEnergyBoundaryBalanceMaxScaled=float(balance.max()),lastLogPhysicalRMS=np.log10(h[-1,4:9]).tolist(),
      allFiveBelowMinus8=bool(np.all(np.log10(h[-1,4:9])<=-8)),
      momentumReactionVsPressureRelative=[float(reaction_error.min()),float(reaction_error.max()),float(reaction_error[-1])],
      fullMomentumBalanceAccepted=False,spatialFarfieldAccepted=False,steadyStableAccepted=False,
      originalFaces=19752,areaM2=float(sum(area)),ranges=stats,newHeatAreaM2=0,newQualifiedAreaM2=0,
      underlyingFailureCauseKnown=False,CFDReplays=0))
    print(json.dumps(read(out/'report.json')))

def export():
    out=BASE/'workbench01';out.mkdir(exist_ok=False);audit=read(BASE/'audit01/report.json');assert audit['passed']
    g=np.load(str(BASE/'audit01/geometry.npz'));n=19752;area=g['area'];data={i:np.loadtxt(str(BASE/'audit01'/('original_faces_%06d.csv'%i)),delimiter=',',skiprows=1) for i in (0,1)}
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json');template=model['fields'][0];fieldsout=[];ranges={}
    for key,label,unit,col,steps in [('pressure','压力','Pa',2,(0,1)),('fluidTemperature','Euler壁节点流体温度','K',3,(0,1)),
      ('totalEnthalpy','总焓（hs298.15=0）','J/kg',4,(1,)),('entropy','同模型比熵','J/(kg·K)',5,(1,)),
      ('velocityX','速度 X','m/s',6,(1,)),('velocityY','速度 Y','m/s',7,(1,)),('velocityZ','速度 Z','m/s',8,(1,))]:
        lo=min(float(data[i][:,col].min()) for i in steps);hi=max(float(data[i][:,col].max()) for i in steps)
        if lo==hi:lo-=1e-12;hi+=1e-12
        for step in steps:
            f=dict(template);name=key+'_'+str(step)
            f.update(key=name,label=label+(' · 初始场' if step==0 else ' · 第1次更新'),unit=unit,values=data[step][:,col].tolist(),
              status=[2 if step==0 else 1]*n,methodIndex=[0]*n,
              methodNames=['SU2 8.6 · 相态Air Euler/general Roe/implicit Euler；原壁三角顶点算术均值；'+('相容给定均匀初态' if step==0 else '原生第一次更新后守恒态反算')],
              scope='M5/30km/α20°；未收敛，第30次更新后异常停止；此完整场仅第%d次更新；非边界层外缘/固体壁温/热解。'%step,
              computedStateLabel='未收敛研究数值',stateNote='初0是给定初态；更新编号不是物理时间；不具备热资格。',availability='available',successfulUpdates=step,
              comparisonVariable=key)
            fieldsout.append(f);ranges[name]=[lo,hi]
    for key,label,unit in [('heatFlux','入壁对流热流','W/m²'),('wallSolved','求解固体壁温','K')]:
        f=dict(template);f.update(key=key,label=label,unit=unit,values=[None]*n,status=[0]*n,methodIndex=[0]*n,
          methodNames=['未求解：无黏外流不提供壁热流或固体温度'],scope='需相容工程边界层/局部黏性反馈及独立固体热模型。',availability='unavailable',
          computedStateLabel='未求解',stateNote='Euler占位零不能作为热流结果。');fieldsout.append(f)
    model.update(vertices=g['vertices'].tolist(),faces=g['faces'].tolist(),faceArea=area.tolist(),totalArea=float(sum(area)),
      region=['原完整箭体']*n,regions=[dict(name='原完整箭体',area=float(sum(area)),faces=n)],
      missingReason=['无黏外流未求解壁热流或固体温度']*n,
      heatingGaps=[dict(reason='工程热层/黏性反馈未接入',areaFraction=1.,faces=n)],fields=fieldsout,
      defaultField='pressure_1',comparisonColorRanges=ranges,resultCondition=dict(caption='M5 · 30 km · α20° · 初0/更新1完整研究场；运行在30步后异常停止'),
      caseNote='只能查看已保存的初始和第1次更新场，不能把它们当作第30步末场或收敛结果。无物理时间；所有热流保持未知。',
      researchIdentity=dict(productionQualified=False,physicalQualification=False,model='rah-su2-phase-air-euler-eplus300k-v1',
        savedUpdates=[0,1],completedUpdates=30,normalCompletion=False,originalFaces=19752,heatFluxSolved=False,
        sourceSHA256={str(i):sha(BASE/'audit01'/('original_faces_%06d.csv'%i)) for i in (0,1)}))
    target=out/'rocket_euler_research.json';target.write_text(json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n',encoding='utf8')
    restored=read(target);assert restored['faces']==g['faces'].tolist()
    for f in restored['fields'][:9]:
        col={'pressure':2,'fluidTemperature':3,'totalEnthalpy':4,'entropy':5,'velocityX':6,'velocityY':7,'velocityZ':8}[f['comparisonVariable']]
        assert np.array_equal(f['values'],data[f['successfulUpdates']][:,col])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    tri=g['vertices'][g['faces']];right=np.array([1.,0.,.1]);right/=np.linalg.norm(right)
    up=np.array([-.025,1.,.25]);up-=right*(up@right);up/=np.linalg.norm(up)
    projected=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(2,1,figsize=(12,6),facecolor='white')
    for ax,col,cmap,label,scale in zip(axs,(2,3),('viridis','inferno'),('Pressure [kPa]','Euler wall-node fluid T [K]'),(.001,1)):
        vals=data[1][:,col]*scale;norm=Normalize(float(vals.min()),float(vals.max()))
        ax.add_collection(PolyCollection(projected[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(vals[order]))))
        ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off();ax.margins(.03)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,shrink=.8,label=label)
    fig.suptitle('RocketAeroHeat | original 19,752 faces | M5 / alpha 20 deg\nComplete saved update 1 only; run stopped abnormally after update 30',fontsize=12)
    fig.text(.08,.025,'Research only: not converged, not a boundary-layer edge. Heat flux and solid temperature remain unknown.',fontsize=10)
    fig.subplots_adjust(top=.86,bottom=.10,left=.04,right=.95,hspace=.10);fig.savefig(str(out/'rocket_euler_research.png'),dpi=160);plt.close(fig)
    test=subprocess.run(['node',str(ROOT/'tests/testPhaseEulerViewerV230.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
    (out/'frontend.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout
    save(out/'report.json',dict(passed=True,fields=len(fieldsout),originalFaces=n,originalAreaM2=float(sum(area)),savedUpdates=[0,1],
      exportSHA256=sha(target),imageSHA256=sha(out/'rocket_euler_research.png'),dataRoundTripExact=True,
      browserInteractionTested=False,sharedWebChanged=False,newHeatAreaM2=0,newQualifiedAreaM2=0))
    print(json.dumps(read(out/'report.json')))

def main(stage):
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes;guard.size_bytes=lambda _:original(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseEulerRocketV230.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',BASE/'process_run/report.json']
    locks += [p for p in CASE.iterdir() if p.is_file()]
    if stage=='export':locks += [ROOT/'web/viewer.js',ROOT/'tests/testPhaseEulerViewerV230.js']
    return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],120,512,locks,BASE,'offline')

if __name__=='__main__':
    if sys.argv[1].startswith('task_'):globals()[sys.argv[1][5:]]()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
