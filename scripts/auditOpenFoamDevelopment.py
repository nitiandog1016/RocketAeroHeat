"""同步发展控制收尾：真实方程残差、壁面参照和科研曲线，资格分开报告。"""
from pathlib import Path
import argparse,json,re,csv,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from openFoamConservation import geometry,flux_and_residual,reference_scales
from openFoamFields import audit_state
from openFoamAirThermo import evaluate

BASE=ROOT/'outputs/v199_validation'


def log_progress(text,inspect_gaps=False):
    """模块一：区分已宣布时刻、日志完整步和实际保存场，避免多算中断步。"""
    announced=[];completed=[];pending=None;orphan=[];completion_count=0
    for line_number,line in enumerate(text.splitlines(),1):
        match=re.match(r'^Time = ([\deE.+-]+)s',line)
        if match:
            if pending is not None:raise ValueError('上一宣布步缺少完成日志')
            pending=float(match.group(1));announced.append(pending)
        elif line.startswith('ExecutionTime = '):
            completion_count+=1
            if pending is None:
                if not inspect_gaps:raise ValueError('完成日志没有对应时刻')
                orphan.append(line_number)
            else:completed.append(pending);pending=None
    if not announced or any(b<=a for a,b in zip(announced,announced[1:])):
        raise ValueError('缺少或倒序物理时间')
    return dict(announcedSteps=len(announced),completedLogSteps=len(completed),
        completionMarkerCount=completion_count,unpairedCompletionLines=orphan,logTimePairingComplete=not orphan,
        lastAnnouncedTimeSeconds=announced[-1],lastCompletedLogTimeSeconds=completed[-1] if completed else None)


def main(inspect_interrupted=False,output_name=None):
    start=time.monotonic();case=BASE/'development02'
    name=output_name or ('development_inspection03' if inspect_interrupted else 'development_audit01')
    if not re.fullmatch(r'[a-zA-Z0-9_]+',name):raise ValueError('只接受本轮目录名')
    out=BASE/name;out.mkdir(exist_ok=False)
    for name in ['auditOpenFoamDevelopment.py','openFoamConservation.py','openFoamFields.py','openFoamAirThermo.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    run=json.loads((BASE/'development_run02/report.json').read_text())
    if not inspect_interrupted:assert run['passed'],'原生控制没有正常完成'
    process=run['process'];assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup']
    log=(BASE/'development_run02/native/stdout.txt').read_text(encoding='utf-8')
    progress=log_progress(log,inspect_gaps=inspect_interrupted)
    assert progress['announcedSteps']<=16000 and 'out of temperature range' not in log
    if not inspect_interrupted:
        assert abs(progress['lastCompletedLogTimeSeconds']-.003)<1e-13
        assert progress['announcedSteps']==progress['completedLogSteps']
    else:assert process['stopReason']=='supervisor-error' and process['nativeExitCode']==-15 and 'read_text' in process['supervisorTraceback']
    refs=json.loads((ROOT/'outputs/v198_validation/independent01/reference.json').read_text(encoding='utf-8'))
    g=geometry(case);vol=g['volume'];ni=len(g['neighbor']);scale,speed=reference_scales()
    wall=next(p for p in g['patches'] if p['name']=='wall');ids=np.arange(wall['startFace'],wall['startFace']+wall['nFaces'])
    owner=g['owner'][ids];area=g['area'][ids];center=g['points'][g['faces'][ids]].mean(axis=1);s=-center[:,0]
    identity=json.loads((case/'wall_identity.json').read_text());assert identity['originalFaceIDs']==list(range(1,241))
    assert np.array_equal(np.sort(g['faces'][ids],axis=1),np.sort(identity['originalConnections'],axis=1))
    frames=[];walls=[]
    for name,row in sorted(run['saved'].items(),key=lambda x:float(x[0])):
        for field,h in row['fieldsSHA256'].items():assert sha(case/name/field)==h
        f,flux,residual=flux_and_residual(case/name,g);state=audit_state(f)
        norms=np.sum(abs(residual)*vol[:,None],axis=0)/vol.sum()/scale/speed
        balance=flux[ni:].sum(axis=0)/(scale*speed) # 全域积分残差，以单位参考面积归一；非瞬时守恒缺陷。
        h0=evaluate(f['p'],f['T'])['enthalpyJkg']+np.einsum('ij,ij->i',f['U'],f['U'])/2
        frames.append(dict(timeSeconds=float(name),state=state,residualVolumeL1Scaled=norms.tolist(),
            residualVolumeL1Log10=np.log10(np.maximum(norms,1e-300)).tolist(),integratedBoundaryFluxScaled=balance.tolist(),
            totalEnthalpyRelativeMax=float(abs(h0/refs['states'][0]['totalEnthalpyJkg']-1).max()),
            totalEnthalpyRelativeVolumeMean=float(np.average(abs(h0/refs['states'][0]['totalEnthalpyJkg']-1),weights=vol))))
        walls.append(np.column_stack([f['p'][owner],f['T'][owner],f['U'][owner]]))
    assert len(frames)==(4 if inspect_interrupted else 6)
    last,previous=walls[-1],walls[-2]
    change={key:float(np.sum(area*abs(last[:,i]-previous[:,i]))/np.sum(area*abs(last[:,i]))) for i,key in enumerate(['pressure','temperature'])}
    comparisons=[]
    for region,lo,hi,ref in [('compression',.35,.5,refs['states'][1]),('expansion',.85,.95,refs['states'][2])]:
        mask=(s>=lo)&(s<=hi);avg=np.average(last[mask],weights=area[mask],axis=0)
        comparisons.append(dict(region=region,streamwiseRangeM=[lo,hi],faces=int(mask.sum()),pressurePa=float(avg[0]),
            temperatureK=float(avg[1]),referencePressurePa=ref['pressurePa'],referenceTemperatureK=ref['temperatureK'],
            pressureRelativeError=float(avg[0]/ref['pressurePa']-1),temperatureRelativeError=float(avg[1]/ref['temperatureK']-1)))
    improvements=np.array(frames[0]['residualVolumeL1Log10'])-np.array(frames[-1]['residualVolumeL1Log10'])
    observed_screening=max(change.values())<=.05 and all(max(abs(x['pressureRelativeError']),abs(x['temperatureRelativeError']))<=.1 for x in comparisons) and min(improvements[[0,4]])>=.3
    screening=run['passed'] and len(frames)==6 and observed_screening
    converged=max(frames[-1]['residualVolumeL1Scaled'])<=1e-8
    result=dict(auditPassed=True,profile=run['profile'],**progress,
        process=process,frames=frames,lastWallAreaL1Change=change,comparisons=comparisons,
        residualLogImprovement=improvements.tolist(),developmentScreeningPassed=bool(screening),
        savedFrameNumericalThresholdPassed=bool(converged),numericallyConverged=False if inspect_interrupted else bool(converged),
        observedIncompleteWindowCriteriaPassed=bool(observed_screening),planAcceptance=run['passed'],
        inspectedInterrupted=inspect_interrupted,lastCompleteSavedTimeSeconds=frames[-1]['timeSeconds'],
        originalNumericalThreshold=1e-8,physicalQualification=False,newQualifiedHeatAreaM2=0,
        originalControlWallExact=True,originalRocketRun=False,elapsedSeconds=time.monotonic()-start)
    with (out/'wall.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['faceID','s_m','y_m','z_m','area_m2','p_Pa','Tfluid_K','Ux_mps','Uy_mps','Uz_mps'])
        for i in range(len(area)):writer.writerow([i+1,s[i],*center[i,1:],area[i],*last[i]])
    save(out/'report.json',result)
    # 模块二：标准科研图只显示定量结果；温度为无黏流体壁邻，非壁温或热流。
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(12,3.4));order=np.argsort(s)
    for i,(key,label) in enumerate([('pressurePa','Pressure (Pa)'),('temperatureK','Adjacent-fluid temperature (K)')]):
        axs[i].scatter(s[order],last[order,i],s=8,label='Native wall owner')
        for row in comparisons:axs[i].hlines(row['reference'+('PressurePa' if i==0 else 'TemperatureK')],*row['streamwiseRangeM'],colors='#b45430',linewidth=2)
        axs[i].set(xlabel='Streamwise s (m)',ylabel=label);axs[i].grid(alpha=.2)
    for j,label in [(0,'Mass'),(4,'Total energy')]:
        axs[2].plot([x['timeSeconds']*1000 for x in frames],[x['residualVolumeL1Log10'][j] for x in frames],'-o',label=label)
    axs[2].axhline(-8,color='gray',linestyle='--',linewidth=1);axs[2].set(xlabel='Physical development time (ms)',ylabel='log10 normalized residual')
    axs[2].legend(frameon=False);axs[2].grid(alpha=.2)
    title='OpenFOAM 13 / NASA7 / synchronized Euler + constant face states'
    if inspect_interrupted:title+='\nInterrupted run | last complete saved field: 2 ms | fluid temperature, no wall heat'
    fig.suptitle(title)
    fig.tight_layout();fig.savefig(out/'control.png',dpi=180);plt.close(fig)
    print(json.dumps({k:result[k] for k in ['auditPassed','announcedSteps','completedLogSteps','lastWallAreaL1Change','comparisons','residualLogImprovement','developmentScreeningPassed','numericallyConverged','elapsedSeconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--inspect-interrupted',action='store_true')
    parser.add_argument('--output-name');args=parser.parse_args()
    main(args.inspect_interrupted,args.output_name)
