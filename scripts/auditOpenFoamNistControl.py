"""低温控制只读审计：同模型真实残差、准确壁面参照及科研图。"""
from pathlib import Path
import csv,json,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from nistIntegratedAir import evaluate,R,MODEL
from openFoamFields import audit_state
from openFoamConservation import geometry,flux_and_residual
from auditOpenFoamDevelopment import log_progress

BASE=ROOT/'outputs/v202_validation'


def main():
    start=time.monotonic();out=BASE/'control_audit01';out.mkdir(exist_ok=False);case=BASE/'control01'
    for name in ['auditOpenFoamNistControl.py','openFoamConservation.py','openFoamFields.py','nistIntegratedAir.py','auditOpenFoamDevelopment.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    run=json.loads((BASE/'control_run01/report.json').read_text(encoding='utf-8'));process=run['process']
    assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup'] and process['logExportPassed']
    refs=json.loads((BASE/'reference02/report.json').read_text(encoding='utf-8'));assert refs['passed'] and refs['gasModel']==MODEL
    log=(BASE/'control_run01/native/stdout.txt').read_text(encoding='utf-8');progress=log_progress(log)
    g=geometry(case);vol=g['volume'];ni=len(g['neighbor']);gas=evaluate(1000.,180.)
    rho=float(gas['densityKgm3']);speed=3*float(gas['soundSpeedMps'])
    scale=np.array([rho,*([rho*speed]*3),rho*(abs(float(gas['internalEnergyJkg']))+speed**2)])
    identity=json.loads((case/'wall_identity.json').read_text());ids=np.array(identity['foamFaceIDsZeroBased']);owner=g['owner'][ids]
    assert identity['originalFaceIDs']==list(range(1,241))
    assert np.array_equal(np.sort(g['faces'][ids],axis=1),np.sort(identity['originalConnections'],axis=1))
    area=g['area'][ids];center=g['points'][g['faces'][ids]].mean(axis=1);s=-center[:,0]
    frames=[];walls=[]
    for name,row in sorted(run['saved'].items(),key=lambda kv:float(kv[0])):
        for field,h in row['fieldsSHA256'].items():assert sha(case/name/field)==h
        f,flux,residual=flux_and_residual(case/name,g,evaluate,R)
        norms=np.sum(abs(residual)*vol[:,None],axis=0)/vol.sum()/scale/speed
        h0=evaluate(f['p'],f['T'])['enthalpyJkg']+np.einsum('ij,ij->i',f['U'],f['U'])/2
        frames.append(dict(timeSeconds=float(name),state=audit_state(f,evaluate),
            residualVolumeL1Scaled=norms.tolist(),residualVolumeL1Log10=np.log10(np.maximum(norms,1e-300)).tolist(),
            integratedBoundaryFluxScaled=(flux[ni:].sum(axis=0)/scale/speed).tolist(),
            totalEnthalpyRelativeMax=float(abs(h0/refs['states'][0]['totalEnthalpyJkg']-1).max())))
        walls.append(np.column_stack([f['p'][owner],f['T'][owner],f['U'][owner]]))
    complete=len(frames)==8 and np.allclose([f['timeSeconds'] for f in frames],np.arange(1,9)*.0005,rtol=0,atol=1e-13)
    changes=None;comparisons=[];improvement=None;converged=False;screening=False
    if len(walls)>=2:
        changes={key:float(np.sum(area*abs(walls[-1][:,i]-walls[-2][:,i]))/np.sum(area*abs(walls[-1][:,i]))) for i,key in enumerate(['pressure','temperature'])}
        for region,lo,hi,ref in [('compression',.35,.5,refs['states'][1]),('expansion',.85,.95,refs['states'][2])]:
            mask=(s>=lo)&(s<=hi);avg=np.average(walls[-1][mask],weights=area[mask],axis=0)
            comparisons.append(dict(region=region,streamwiseRangeM=[lo,hi],faces=int(mask.sum()),pressurePa=float(avg[0]),temperatureK=float(avg[1]),
                referencePressurePa=ref['pressurePa'],referenceTemperatureK=ref['temperatureK'],
                pressureRelativeError=float(avg[0]/ref['pressurePa']-1),temperatureRelativeError=float(avg[1]/ref['temperatureK']-1)))
        improvement=np.array(frames[0]['residualVolumeL1Log10'])-np.array(frames[-1]['residualVolumeL1Log10'])
        screening=max(changes.values())<=.05 and min(improvement[[0,4]])>=.3 and all(max(abs(r['pressureRelativeError']),abs(r['temperatureRelativeError']))<=.1 for r in comparisons)
        converged=bool(max(frames[-1]['residualVolumeL1Scaled'])<=1e-8)
    result=dict(passed=bool(run['passed'] and complete and screening and converged),model=MODEL,
        planAcceptance=run['passed'],actualSolverDispatchVerified=run['actualSolverDispatchVerified'],progress=progress,
        completeSavedTimes=bool(complete),frames=frames,scalarScales=scale.tolist(),velocityScaleMps=speed,
        scaleDefinition='rho_inf, rho_inf U_inf (each momentum), rho_inf(abs(es_inf)+U_inf^2); L=1m',
        lastWallAreaL1Change=changes,comparisons=comparisons,residualLogImprovement=None if improvement is None else improvement.tolist(),
        numericallyConverged=converged,developmentScreeningPassed=bool(screening),originalControlWallExact=True,
        originalRocketRun=False,physicalQualification=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    if walls:
        with (out/'wall.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.writer(stream);writer.writerow(['faceID','s_m','y_m','z_m','area_m2','p_Pa','Tfluid_K','Ux_mps','Uy_mps','Uz_mps'])
            for i in range(len(ids)):writer.writerow([i+1,s[i],*center[i,1:],area[i],*walls[-1][i]])
        # 模块二：保持量纲/参照和温度身份；无黏控制不展示壁热或固体温度。
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,3,figsize=(12,3.7))
        for i,label in enumerate(['Pressure (Pa)','Adjacent-fluid temperature (K)']):
            axs[i].scatter(s,walls[-1][:,i],s=7,label='Native wall owner')
            for row in comparisons:axs[i].hlines(row['referencePressurePa' if i==0 else 'referenceTemperatureK'],*row['streamwiseRangeM'],color='#b45430',lw=2)
            axs[i].set(xlabel='Streamwise s (m)',ylabel=label);axs[i].grid(alpha=.2)
        for j,label in [(0,'Mass'),(4,'Total energy')]:axs[2].plot([f['timeSeconds']*1000 for f in frames],[f['residualVolumeL1Log10'][j] for f in frames],'-o',label=label)
        axs[2].axhline(-8,color='gray',ls='--');axs[2].set(xlabel='Physical development time (ms)',ylabel='log10 normalized residual');axs[2].legend(frameon=False);axs[2].grid(alpha=.2)
        fig.suptitle('OpenFOAM 13 / NIST integrated Cp / synchronized Euler + constant face states\nM3, 180 K, 1000 Pa | fluid temperature, no wall heat')
        fig.tight_layout();fig.savefig(out/'control.png',dpi=180);plt.close(fig)
    print(json.dumps({k:result[k] for k in ['passed','completeSavedTimes','numericallyConverged','developmentScreeningPassed','lastWallAreaL1Change','comparisons','elapsedSeconds']}))


if __name__=='__main__':main()
