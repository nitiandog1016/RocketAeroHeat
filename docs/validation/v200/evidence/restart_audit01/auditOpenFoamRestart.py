"""只读验收正常恢复及六保存时刻；保留旧失败、数值与物理资格的区别。"""
from pathlib import Path
import csv,json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from openFoamConservation import geometry,flux_and_residual,reference_scales
from openFoamFields import audit_state
from openFoamAirThermo import evaluate
from auditOpenFoamDevelopment import log_progress
from runOpenFoamRestart import validate_source

BASE=ROOT/'outputs/v200_validation'


def main():
    start=time.monotonic();out=BASE/'restart_audit01';out.mkdir(exist_ok=False);case=BASE/'restart01'
    for name in ['auditOpenFoamRestart.py','openFoamConservation.py','auditOpenFoamDevelopment.py','runOpenFoamRestart.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    run=json.loads((BASE/'restart_run01/report.json').read_text());assert run['passed']
    p=run['process'];assert p['nativeExitCode']==0 and p['stopReason'] is None
    assert all(p['inputStable'].values()) and not p['liveGroupAfterCleanup'] and p['logExportPassed']
    assert not p['observationIORetries'] and p['liveLogFilesystem'] not in ['9p','drvfs']
    validate_source(run['source']);validate_source(run['source'],case)
    log=(BASE/'restart_run01/native/stdout.txt').read_text();progress=log_progress(log)
    assert progress['announcedSteps']==progress['completedLogSteps']<=7000 and progress['lastCompletedLogTimeSeconds']==.003
    oldpath=ROOT/'outputs/v199_validation/development_inspection03/report.json'
    assert sha(oldpath)==sha(ROOT/'docs/validation/v199/evidence/development_inspection03/report.json')
    old=json.loads(oldpath.read_text());assert not old['planAcceptance']
    g=geometry(case);vol=g['volume'];ni=len(g['neighbor']);scale,speed=reference_scales()
    identity=json.loads((case/'wall_identity.json').read_text());ids=np.array(identity['foamFaceIDsZeroBased']);owner=g['owner'][ids]
    assert identity['originalFaceIDs']==list(range(1,241))
    assert np.array_equal(np.sort(g['faces'][ids],axis=1),np.sort(identity['originalConnections'],axis=1))
    area=g['area'][ids];center=g['points'][g['faces'][ids]].mean(axis=1);s=-center[:,0]
    refs=json.loads((ROOT/'outputs/v198_validation/independent01/reference.json').read_text());frames=list(old['frames']);walls=[]
    for name in ['0.0025','0.003']:
        for field,h in run['saved'][name]['fieldsSHA256'].items():assert sha(case/name/field)==h
        f,flux,residual=flux_and_residual(case/name,g);norms=np.sum(abs(residual)*vol[:,None],axis=0)/vol.sum()/scale/speed
        h0=evaluate(f['p'],f['T'])['enthalpyJkg']+np.einsum('ij,ij->i',f['U'],f['U'])/2
        frames.append(dict(timeSeconds=float(name),state=audit_state(f),residualVolumeL1Scaled=norms.tolist(),
            residualVolumeL1Log10=np.log10(np.maximum(norms,1e-300)).tolist(),
            integratedBoundaryFluxScaled=(flux[ni:].sum(axis=0)/scale/speed).tolist(),
            totalEnthalpyRelativeMax=float(abs(h0/refs['states'][0]['totalEnthalpyJkg']-1).max())))
        walls.append(np.column_stack([f['p'][owner],f['T'][owner],f['U'][owner]]))
    finaltime=(case/'0.003/uniform/time').read_text();index=int(re.search(r'\bindex\s+(\d+);',finaltime).group(1))
    assert index==7889+progress['completedLogSteps']
    first=float(re.search(r'^Time = ([\deE.+-]+)s',log,re.M).group(1));assert .002<first<.002001
    change={key:float(np.sum(area*abs(walls[-1][:,i]-walls[-2][:,i]))/np.sum(area*abs(walls[-1][:,i]))) for i,key in enumerate(['pressure','temperature'])}
    comparisons=[]
    for region,lo,hi,ref in [('compression',.35,.5,refs['states'][1]),('expansion',.85,.95,refs['states'][2])]:
        mask=(s>=lo)&(s<=hi);avg=np.average(walls[-1][mask],weights=area[mask],axis=0)
        comparisons.append(dict(region=region,pressurePa=float(avg[0]),temperatureK=float(avg[1]),
            referencePressurePa=ref['pressurePa'],referenceTemperatureK=ref['temperatureK'],streamwiseRangeM=[lo,hi],
            pressureRelativeError=float(avg[0]/ref['pressurePa']-1),temperatureRelativeError=float(avg[1]/ref['temperatureK']-1)))
    improvement=np.array(frames[0]['residualVolumeL1Log10'])-np.array(frames[-1]['residualVolumeL1Log10'])
    converged=max(frames[-1]['residualVolumeL1Scaled'])<=1e-8
    screening=max(change.values())<=.05 and min(improvement[[0,4]])>=.3 and all(max(abs(x['pressureRelativeError']),abs(x['temperatureRelativeError']))<=.1 for x in comparisons)
    assert [r['timeSeconds'] for r in frames]==[.0005,.001,.0015,.002,.0025,.003]
    result=dict(passed=bool(converged and screening),normalRestart=True,oldRunPlanAcceptance=False,newRestartPlanAcceptance=True,
        combinedSavedTimesComplete=True,continuousBitwiseEquivalent=False,progress=progress,firstAnnouncedTimeSeconds=first,
        sourceTimeIndex=7889,finalTimeIndex=index,frames=frames,lastWallAreaL1Change=change,comparisons=comparisons,
        residualLogImprovement=improvement.tolist(),numericallyConverged=bool(converged),developmentScreeningPassed=bool(screening),
        originalControlWallExact=True,physicalQualification=False,newQualifiedHeatAreaM2=0,sourceReportSHA256=sha(oldpath),
        newReportSHA256=sha(BASE/'restart_run01/report.json'),elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    with (out/'wall.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['faceID','s_m','y_m','z_m','area_m2','p_Pa','Tfluid_K','Ux_mps','Uy_mps','Uz_mps'])
        for i in range(len(ids)):writer.writerow([i+1,s[i],*center[i,1:],area[i],*walls[-1][i]])
    # 模块二：展示真实控制压力/流体温度及恢复前后残差，不显示不存在的壁热。
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.6));order=np.argsort(s)
    for i,label in enumerate(['Pressure (Pa)','Adjacent-fluid temperature (K)']):
        axes[i].scatter(s[order],walls[-1][order,i],s=7)
        for row in comparisons:axes[i].hlines(row['referencePressurePa' if i==0 else 'referenceTemperatureK'],*row['streamwiseRangeM'],color='#b45430',lw=2)
        axes[i].set(xlabel='Streamwise s (m)',ylabel=label);axes[i].grid(alpha=.2)
    for j,label in [(0,'Mass'),(4,'Total energy')]:axes[2].plot([f['timeSeconds']*1000 for f in frames],[f['residualVolumeL1Log10'][j] for f in frames],'-o',label=label)
    axes[2].axhline(-8,color='gray',ls='--');axes[2].axvline(2,color='#999999',ls=':',label='Restart source')
    axes[2].set(xlabel='Physical development time (ms)',ylabel='log10 normalized residual');axes[2].legend(frameon=False);axes[2].grid(alpha=.2)
    fig.suptitle('OpenFOAM 13 / NASA7 / synchronized Euler + constant face states\nNormal saved-field restart to 3 ms | fluid temperature, no wall heat')
    fig.tight_layout();fig.savefig(out/'control.png',dpi=180);plt.close(fig)
    print(json.dumps({k:result[k] for k in ['passed','progress','lastWallAreaL1Change','comparisons','numericallyConverged','elapsedSeconds']}))


if __name__=='__main__':main()
