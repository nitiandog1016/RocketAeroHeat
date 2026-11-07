"""修正后真实Euler异常场只读审核；不绕过热物性拒绝、不输出可用外缘。"""
from pathlib import Path
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import sha,read_grid,read_field,compare_source
from hybridRocketEuler import check_inputs,validate_terminal
from hybridEdgeSampling import map_original_wall
from eilmerAirThermo import evaluate
from prepareHybridRocketEulerNative import local,save

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v196_validation'

def main():
    """模块一：完整保留原始异常值，分别统计接口超域与非物理状态。"""
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_audit.py')
    case=BASE/'rocket01';rp=BASE/'run01/report.json';prep=BASE/'preparation01/report.json'
    r=check_inputs(case,prep);run=json.loads(rp.read_text());p=json.loads((rp.parent/'native/report.json').read_text())
    assert p['inputsStable'] and not p['liveGroupAfterCleanup']
    for name,d in run['outputSHA256'].items():assert sha(local(name))==d
    rejected=False
    try:validate_terminal(case,prep,rp)
    except ValueError:rejected=True
    assert rejected and not run['processPassed']
    text=(rp.parent/'native/stdout.txt').read_text()
    assert 'Explicit update failed after 1 attempts' in text and 'FINAL-STEP: 23' in text
    assert 'RAH-STORAGE cells=64982 ghosts=26324 neq=5 nftl=4 scalarBytes=8 stageLength=456530 payloadBytes=29217920' in text
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,_=compare_source(grid,case/'mesh.su2')
    surface=json.loads((case/'surface.json').read_text());V=np.array(surface['vertices']);F=np.array(surface['facesOneBased'])-1
    tri=V[F];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert len(area)==19752 and abs(area.sum()-1.9646292565715551)<2e-13
    entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)\s*\n\s+dt:\s*(\S+)",
                       (case/'lmrsim/snapshots/snapshot-times-metadata').read_text(),re.M)
    records=[]
    for index,step,t,dt in entries:
        f=read_field(case,int(index),64982)
        location=map_original_wall(surface,grid,mapping['wall'],np.column_stack([f['pos.'+k] for k in 'xyz']))
        owners=np.array(location['nativeOwnerCellsZeroBased'])
        eligible=(f['T']>=200)&(f['T']<=6000)&(f['p']>0)&(f['rho']>0)
        gas=evaluate(f['p'][eligible],f['T'][eligible])
        errors={key:float(np.max(abs(f[key][eligible]-gas[target])/np.maximum(1,abs(gas[target]))))
            for key,target in [('rho','densityKgm3'),('e','internalEnergyJkg'),('a','soundSpeedMps')]}
        worst=int(np.argmin(f['T']));ids=np.flatnonzero(~eligible)
        records.append(dict(snapshot=int(index),reportedStep=int(step),pseudoTime=float(t),savedDt=float(dt),
            temperatureRangeK=[float(f['T'].min()),float(f['T'].max())],pressureRangePa=[float(f['p'].min()),float(f['p'].max())],
            negativeTemperatureCount=int(np.sum(f['T']<=0)),negativePressureCount=int(np.sum(f['p']<=0)),
            negativeDensityCount=int(np.sum(f['rho']<=0)),outside200To6000Count=int(np.sum((f['T']<200)|(f['T']>6000))),
            eligibleThermoCells=int(eligible.sum()),subsetEOSErrors=errors,allCellEOSPassed=bool(eligible.all() and max(errors.values())<2e-10),
            minTemperatureCellZeroBased=worst,minTemperaturePositionM=[float(f['pos.'+k][worst]) for k in 'xyz'],
            originalFacesWithRejectedOwner=int(np.sum(~eligible[owners])),
            originalAreaWithRejectedOwnerM2=float(area[~eligible[owners]].sum())))
    assert [r['reportedStep'] for r in records]==[0,23]
    # 模块二：原日志时间尺度已有限；异常保存帧的dt是失败路径缩减值，不能当成功步长。
    history=np.atleast_2d(np.loadtxt(case/'lmrsim/transient-residuals.dat'))
    assert history.shape==(23,10) and np.isfinite(history).all()
    live=np.array([(int(a),float(b)) for a,b in re.findall(r'Step=\s*(\d+)\s+t=\s*\S+\s+dt=\s*(\S+)',text)])
    assert len(live)==23 and np.all((live[:,1]>0)&(live[:,1]<.001))
    free=json.loads((case/'native-input-state.json').read_text());speed=np.linalg.norm([free[k] for k in ['velx','vely','velz']])
    L=float(np.ptp(V[:,0]));rho=free['rho'];scale=np.array([rho*speed/L]+[rho*speed**2/L]*3+[rho*speed**3/L])
    normalized=history[:,3:8]/scale
    np.savetxt(out/'normalized_residuals.csv',np.c_[history[:,0],normalized],delimiter=',',comments='',
               header='step,mass,x_momentum,y_momentum,z_momentum,energy',fmt='%.17g')
    np.savetxt(out/'rejected_wall_diagnostic.csv',np.c_[np.arange(1,19753),owners,area,tri.mean(1),f['p'][owners],f['T'][owners],eligible[owners]],
        delimiter=',',comments='',fmt='%.17g',header='face_id,owner_zero_based,area_m2,x_m,y_m,z_m,p_raw_Pa,T_raw_K,thermo_input_eligible_only')
    # 模块三：科研诊断图明确全候选拒绝；可见颜色不是生产/物理资格。
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig=plt.figure(figsize=(12,5));ax=fig.add_subplot(121,projection='3d')
    temp=f['T'][owners];norm=matplotlib.colors.Normalize(200,2200)
    colors=plt.cm.cividis(norm(temp));colors[~eligible[owners]]=[.8,.15,.2,1]
    ax.add_collection3d(Poly3DCollection(tri,facecolors=colors,linewidth=0))
    lo,hi=V.min(0),V.max(0);ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]),xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    ax.set_box_aspect(hi-lo);ax.view_init(20,-70);ax.set_title('REJECTED exception field\nRed: owner outside declared thermodynamics')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='cividis'),ax=ax,shrink=.55,pad=.1,label='Raw fluid T (K), not wall T')
    ax=fig.add_subplot(122)
    for i,label in enumerate(['mass','x momentum','y momentum','z momentum','energy']):
        ax.semilogy(history[:,0],normalized[:,i],label=label)
    ax.axhline(1e-8,color='gray',linestyle='--',label='Required 1e-8')
    ax.set(xlabel='Reported step',ylabel='Dimensionless residual',title='Time scale fixed; development failed')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'rejected_development.png',dpi=180);plt.close(fig)
    report=dict(passed=True,scope='read-only failure audit, not candidate acceptance',candidateAccepted=False,
        productionExportRejected=True,numericalConvergencePassed=False,physicalQualification=False,newSolvedHeatAreaM2=0,
        frames=records,historyRows=23,firstReportedDtS=float(live[0,1]),liveReportedDtRangeS=[float(live[:,1].min()),float(live[:,1].max())],
        finalNormalizedResiduals=normalized[-1].tolist(),full50WindowBenefit=None,
        nativeException='Explicit update failed after 1 attempts',exceptionStageAndFirstInvalidQuery='not recorded in fast-build log',
        savedFrameIdentity='exception-path state, not an accepted completed step',
        actualStop='checkpoint guard observed invalid saved field after native exception log; supervisor SIGTERM',
        nativeExitCode=p['nativeExitCode'],stopReason=p['stopReason'],nativeSupervisedSeconds=p['elapsedSeconds'],
        peakGroupRSSBytes=p['peakObservedGroupRSSBytes'],inputCount=len(p['inputSHA256']),inputsStable=True,
        originalFaces=19752,originalAreaM2=float(area.sum()),elapsedSeconds=time.monotonic()-start,
        sourceSHA256={str(q):sha(q) for q in [Path(__file__),rp,prep,rp.parent/'native/report.json']})
    save(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['frames','sourceSHA256']}))

if __name__=='__main__':main()
