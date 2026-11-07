"""真实NK候选的独立终态审核：保留失败场，不将研究值导出为可用外缘。"""
from pathlib import Path
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import sha,need,read_grid,read_field,compare_source
from prepareHybridRocketEulerNative import local,save
from hybridEdgeSampling import map_original_wall
from eilmerAirThermo import evaluate
from auditHybridNKControl import read_diagnostics,snapshot_steps,normal_terminal

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v197_validation'


def main():
    """模块一：仅对原生已终态任务读取，按保存步、原面与同EOS分别审核。"""
    start=time.monotonic();case=BASE/'rocket01';rp=BASE/'rocket_run01/report.json'
    r=json.loads(rp.read_text());p=r['process'];out=BASE/'rocket_audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_audit.py')
    need(p['inputsStable'] and not p['liveGroupAfterCleanup'],'输入改变或进程未清点')
    for path,h in r['outputSHA256'].items():need(sha(local(path))==h,'原生结果改变')
    text=(rp.parent/'native/stdout.txt').read_text()
    need(re.findall(r'^Number-type: (.+)',text,re.M)==['complex'],'实际不是复数NK')
    rows=(case/'lmrsim/diagnostics/nk-diagnostics').read_text().splitlines()
    count=len([line for line in rows[1:] if line.strip()]);hist=read_diagnostics(case/'lmrsim/diagnostics/nk-diagnostics',count) if count else None
    metadata=case/'lmrsim/snapshots/restart-metadata'
    steps=snapshot_steps(case,hist) if metadata.exists() else [0]
    indices=[int(d.name) for d in sorted((case/'lmrsim/snapshots').glob('[0-9][0-9][0-9][0-9]'))]
    need(indices==list(range(len(steps))),'保存帧缺失或未经元数据绑定')
    normal=False;reason=None
    try:normal_terminal(p,text,120);normal=True
    except ValueError as error:reason=str(error)
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,_=compare_source(grid,case/'mesh.su2')
    need(len(grid['cells'])==64982 and np.all(grid['kinds']==10),'原体网格改变')
    surface=json.loads((case/'surface.json').read_text());v=np.array(surface['vertices']);tri=v[np.array(surface['facesOneBased'])-1]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    need(len(area)==19752 and abs(area.sum()-1.9646292565715553)<2e-13,'原面积改变')
    initial=read_field(case,0,64982);positions=np.column_stack([initial['pos.'+a] for a in 'xyz'])
    free=json.loads((case/'native-input-state.json').read_text())
    h0free=float(evaluate(free['p'],free['T'])['enthalpyJkg']+sum(free['vel'+a]**2 for a in 'xyz')/2)
    loc=map_original_wall(surface,grid,mapping['wall'],positions);owners=np.array(loc['nativeOwnerCellsZeroBased'])
    records=[];walls=[]
    for snapshot,step in enumerate(steps):
        f=read_field(case,snapshot,64982)
        need(np.all(f['vol']>0) and np.all(f['massf-air']==1),'网格体积/气体身份改变')
        need(all(np.array_equal(f[key],initial[key]) for key in ['pos.x','pos.y','pos.z','vol']),'保存场位置/体积改变')
        eligible=(f['T']>=200)&(f['T']<=6000)&(f['p']>0)&(f['rho']>0)
        gas=evaluate(f['p'][eligible],f['T'][eligible]);errors={}
        for name,key in [('rho','densityKgm3'),('e','internalEnergyJkg'),('a','soundSpeedMps')]:
            errors[name]=float(np.max(abs(f[name][eligible]-gas[key])/np.maximum(1,abs(gas[key])))) if eligible.any() else None
        h0=gas['enthalpyJkg']+sum(f['vel.'+a][eligible]**2 for a in 'xyz')/2
        j=int(np.argmin(f['T']));records.append(dict(snapshot=snapshot,reportedStep=step,
            temperatureRangeK=[float(f['T'].min()),float(f['T'].max())],minimumPressurePa=float(f['p'].min()),
            minimumDensityKgm3=float(f['rho'].min()),outsideTemperatureCount=int(np.sum((f['T']<200)|(f['T']>6000))),
            nonpositivePressureCount=int(np.sum(f['p']<=0)),nonpositiveDensityCount=int(np.sum(f['rho']<=0)),
            thermoEligibleCells=int(eligible.sum()),eligibleSubsetEOSError=errors,
            eligibleTotalEnthalpyToFreestreamRange=[float(h0.min()/h0free),float(h0.max()/h0free)] if eligible.any() else None,
            totalEnthalpyDiagnosticScope='cell values only, not global boundary-flux energy conservation proof',
            allCellEOSPassed=bool(eligible.all() and max(errors.values())<=2e-10),
            minTemperatureCellZeroBased=j,minTemperaturePositionM=positions[j].tolist(),
            rejectedOwnerFaceCount=int(np.sum(~eligible[owners])),rejectedOwnerAreaM2=float(area[~eligible[owners]].sum())))
        walls.append(np.column_stack([f['p'][owners],f['T'][owners]]))
    improvement=None
    if count>=40:
        improvement={key:float(np.log10(hist[col][:20]).mean()-np.log10(hist[col][-20:]).mean())
                     for key,col in [('mass','mass-abs'),('energy','total-energy-abs')]}
    relative=None;converged=False;changes=None
    if count:
        relative={name:float(np.max(hist[name][-min(count,20):])) for name in
                  ['mass-rel','x-mom-rel','y-mom-rel','z-mom-rel','total-energy-rel']}
        converged=normal and max(relative.values())<=1e-8
    if len(walls)>1:
        changes={key:float(np.sum(abs(walls[-1][:,i]-walls[-2][:,i])*area)/np.sum(abs(walls[-1][:,i])*area))
                 for i,key in enumerate(['pressure','temperature'])}
    screening=bool(normal and steps==list(range(0,121,20)) and all(x['allCellEOSPassed'] for x in records)
                   and improvement is not None and min(improvement.values())>=.05)
    # 模块二：明确诊断身份，不调用生产外缘交换；异常数值完整保留，q/Tw不补零。
    np.savetxt(out/'wall_diagnostic.csv',np.c_[np.arange(1,19753),owners,area,tri.mean(1),walls[-1],eligible[owners]],
        delimiter=',',comments='',fmt='%.17g',header='face_id,owner_zero_based,area_m2,x_m,y_m,z_m,p_raw_Pa,T_raw_fluid_K,thermo_eligible_only')
    result=dict(passed=True,scope='read-only development audit, not physics or engineering acceptance',
        candidateScreeningPassed=screening,normalTerminal=normal,terminalRejectionReason=reason,
        nativeExitCode=p['nativeExitCode'],stopReason=p['stopReason'],reportedFinalSteps=re.findall(r'^FINAL-STEP: (.+)',text,re.M),
        nativeStopLabels=re.findall(r'^STOP-REASON: (.+)',text,re.M),historyRows=count,savedSteps=steps,frames=records,
        meanLog10ImprovementLast20VsFirst20=improvement,nativeRelativeResidualLast20Max=relative,
        nativeGlobalRelativeResidualFinal=float(hist['global-residual-rel'][-1]) if count else None,
        nativeMassBalanceFinal=float(hist['mass-balance'][-1]) if count else None,
        relaxationFactorRange=[float(hist['relaxation-factor'].min()),float(hist['relaxation-factor'].max())] if count else None,
        finalCFL=float(hist['CFL'][-1]) if count else None,lastPairWallAreaL1=changes,
        nativeRelativeResidualGatePassed=converged,numericalConvergencePassed=False,physicalQualification=False,
        attachedEdgeQualified=False,engineeringHeatReady=False,newSolvedHeatAreaM2=0,
        qIdentity='unknown',solidTwIdentity='unknown',originalFaces=19752,originalAreaM2=float(area.sum()),
        processSeconds=p['elapsedSeconds'],entrySeconds=r['elapsedSeconds'],peakGroupRSSBytes=p['peakObservedGroupRSSBytes'],
        inputsStable=p['inputsStable'],stableInputCount=len(p['inputStable']),elapsedSeconds=time.monotonic()-start,
        sourcesSHA256={str(path):sha(path) for path in [Path(__file__),rp,BASE/'rocket_preparation01/report.json']})
    save(out/'report.json',result)
    print(json.dumps({key:result[key] for key in ['passed','candidateScreeningPassed','normalTerminal','savedSteps','historyRows',
        'meanLog10ImprovementLast20VsFirst20','elapsedSeconds']}))
    return result


if __name__=='__main__':main()
