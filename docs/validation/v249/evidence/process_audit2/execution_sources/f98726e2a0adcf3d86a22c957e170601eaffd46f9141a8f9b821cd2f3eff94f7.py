"""模块职责：独立验证真实单元星、原生NK初末场/完整阶段及原面结果；资源失败独立列出。"""
from pathlib import Path
import csv,subprocess,sys
import numpy as np
from scipy.linalg import qr,svd
from scipy.io import loadmat
from scipy.optimize import brentq
import prepareOuterStencilV249 as task
from terminalCapacityV249 import supervise
from singleDomainFrame import read_frame
from fasRecoveryFieldsV248 import state
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=task.BASE;DOC=task.DOC;CASE=BASE/'native01';CLOSED=False
VOLUME=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def rows(path):
    with Path(path).open(encoding='utf8') as stream:return list(csv.DictReader(stream))


def audit():
    out=BASE/'audit02';out.mkdir(exist_ok=False);errors={}
    def check(key,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),key
        value=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[key]=value;assert value<tol,(key,value)
    # 模块一：原MAT母网格与SU2文本导出独立对照，包含完整入射单元而非仅形状/数量。
    original=loadmat(str(VOLUME),variable_names=['vertices','wedges','tetrahedra']);xyz=original['vertices']
    wedges=original['wedges'].astype(np.int64)-1;tets=original['tetrahedra'].astype(np.int64)-1
    z={k:v for k,v in np.load(str(BASE/'stencil01/actual_stencils.npz')).items()};report=read(BASE/'stencil01/report.json')
    assert sha(BASE/'stencil01/actual_stencils.npz')==report['artifactSHA256']
    ids=z['globalIDs'];assert np.array_equal(z['xyz'],xyz[ids]);starunion=set()
    for row in report['seeds']:
        seed=row['globalID'];w=np.flatnonzero((wedges==seed).any(axis=1));t=np.flatnonzero((tets==seed).any(axis=1))
        expected=np.r_[w,t+len(wedges)];assert np.array_equal(expected,row['cellIDs']);starunion.update(expected)
        nbset=set(np.unique(tets[t]).tolist());nbset.discard(seed)
        # 棱柱从三角面/四边面闭环导出实际边；不把四边面的对角线当作相邻边。
        for nodes in wedges[w]:
            for face in ((0,1,2),(3,5,4),(0,3,4,1),(1,4,5,2),(2,5,3,0)):
                for i,j in zip(face,face[1:]+face[:1]):
                    if nodes[i]==seed:nbset.add(int(nodes[j]))
                    if nodes[j]==seed:nbset.add(int(nodes[i]))
        nb=np.array(sorted(nbset));assert np.array_equal(nb,row['edgeNeighbors'])
        dx=xyz[nb]-xyz[seed]
        normed=dx/np.linalg.norm(dx,axis=1)[:,None];upper=qr(normed,mode='economic')[1]
        values=svd(upper,compute_uv=False,lapack_driver='gesvd');condition=float(values[0]/values[-1])
        check('originalGridQRCondition_%d'%seed,row['normalizedWLSDesignCondition'],condition)
        e=xyz[tets[t,1:]]-xyz[tets[t,0],None];vol=abs(np.sum(np.cross(e[:,0],e[:,1])*e[:,2],axis=1))/6
        if len(vol):check('independentCellVolume_%d'%seed,row['minimumTetVolumeM3'],min(vol))
        else:assert row['minimumTetVolumeM3'] is None
    assert sorted(starunion)==z['motherCells'][:,0].tolist()
    wc=z['motherCells'][z['motherCells'][:,1]==13];tc=z['motherCells'][z['motherCells'][:,1]==10]
    assert np.array_equal(wc[:,2:8],wedges[wc[:,0]]) and np.array_equal(tc[:,2:6],tets[tc[:,0]-len(wedges)])
    fine_records=rows(ROOT/'outputs/v246_validation/initial01/boundary.csv')
    boundary_ids={int(r['globalID']) for r in fine_records if int(r['level'])==0}
    for row in report['seeds']:assert bool(row['boundaryLabels'])==(row['globalID'] in boundary_ids)
    # 模块二：完整真实原生控制；该步骤的工作RHS不是末场残差。
    process=read(BASE/'process_native/report.json');assert process['passed'] and process['inputsStable'] and not process['liveGroupAfterCleanup']
    capacity=read(BASE/'process_native/terminal_capacity.json');assert not capacity['passed'] and capacity['originalSupervisorPassed']
    done=read(CASE/'native_completion.json');assert done['actualNewtonFactory'] and done['updates']==1 and done['SAEvaluations']==1
    stages=rows(CASE/'stages.csv');flow=[r for r in stages if int(r['system'])==0];sa=[r for r in stages if int(r['system'])==2]
    pattern=['after_sources','after_weak','before_strong','after_strong']
    assert [r['stage'] for r in flow]==pattern*done['flowSpaceEvaluations']
    assert [r['stage'] for r in sa]==pattern+['before_time','after_time']
    assert len(rows(CASE/'rejections.csv'))==0 and len(rows(CASE/'steps.csv'))==1
    history=rows(CASE/'history.csv');h={k.strip().strip('"'):v.strip() for k,v in history[-1].items()}
    linear=[(k,v) for k,v in h.items() if 'LinSol' in k or 'Lin_Res' in k or 'Linear' in k]
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};area=g['faceArea'];tri=g['triangles'];wall=g['wallIDs'];frames=[];air=PhaseAwareViscousAir()
    reference=state(task.PREV/'rocket01/state_000200.bin');reference=reference[np.argsort(reference[:,0])]
    for step in (0,1):
        a=read_frame(CASE/('state_%06d.bin'%step),335726,step==0)['values'];a=a[np.argsort(a[:,0])]
        assert np.array_equal(a[:,1:4],g['points'])
        if step==0:assert np.array_equal(a[:,:11].copy().view(np.uint64),reference[:,:11].copy().view(np.uint64))
        velocity=a[:,6:9]/a[:,5,None];e=a[:,9]/a[:,5]-.5*np.sum(velocity**2,axis=1)-300000;t=inverse_energy(air,e)
        phys=air.evaluate(a[:,5]*air.R*t,t);transport=air.interpolated_transport(a[:,15],a[:,11])
        check('temperature%d'%step,a[:,11],t);check('pressure%d'%step,a[:,15],phys['pressurePa'])
        check('velocity%d'%step,a[:,12:15],velocity);check('conductivity%d'%step,a[:,21],transport['conductivityWmK'])
        check('physicalH%d'%step,a[:,17]-300000,np.asarray(phys['enthalpyJkg'])+.5*np.sum(velocity**2,axis=1))
        nu=a[:,10];chi=nu/(np.asarray(transport['dynamicViscosityPas'])/a[:,5]);mut=np.where(nu<0,0,a[:,5]*nu*chi**3/(chi**3+7.1**3))
        check('updatedSAMut%d'%step,a[:,23],mut)
        selected=np.unique(np.r_[np.linspace(0,len(t)-1,11).astype(int),np.argmin(t),np.argmax(t)])
        independent=[brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-e[i],60,2000,xtol=1e-10) for i in selected]
        check('Brent%d'%step,t[selected],independent)
        q=np.full(len(a),np.nan);q[wall]=a[wall,21]*(a[g['normalNeighbor'],11]-300)/g['thermalDistance']
        faces=np.column_stack((np.arange(1,19753),area,a[tri,15].mean(axis=1),q[tri].mean(axis=1),a[tri,11].mean(axis=1),np.full(19752,300)))
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),faces,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
        frames.append(dict(update=step,Tmin=float(min(t)),Tmax=float(max(t)),pressureRange=[float(min(faces[:,2])),float(max(faces[:,2]))],
          qRange=[float(min(faces[:,3])),float(max(faces[:,3]))],researchPowerIntoWallW=float(area@faces[:,3])))
        if step==1:
            native=np.loadtxt(str(CASE/'restart_flow.csv'),delimiter=',',skiprows=1,usecols=range(10));native=native[np.argsort(native[:,0])]
            assert np.array_equal(native[:,:4],a[:,:4]) and np.array_equal(native[:,4:10].copy().view(np.uint64),a[:,5:11].copy().view(np.uint64))
    result=dict(passed=True,originalMeshIndependentStarPassed=True,patchNodes=len(ids),patchCells=len(starunion),
      leadingPointInterior=not report['seeds'][0]['boundaryLabels'],allSeedBoundaryCounts=sum(bool(r['boundaryLabels']) for r in report['seeds']),
      nativeControlCompleted=True,completeUpdates=1,actualFactory='CNewtonIntegration',flowSpaceEvaluations=done['flowSpaceEvaluations'],
      inferredMatrixFreeProducts=done['flowSpaceEvaluations']-1,inferenceUsesFixedSourceAndStartupDisabled=True,
      nativeLinearHistory=linear,flowStageRows=len(flow),SAStageRows=len(sa),rejections=0,
      initialElevenBitEqual=True,finalNativeSixCoordinateBitEqual=True,independentFineStates=671452,frames=frames,errorMetrics=errors,
      originalCapacityPassed=False,originalCapacityReceiptPreserved=True,auditOutputLimitMiB=1024,
      workingRHSIsFinalResidual=False,nonlinearConvergenceProved=False,fullJacobianAccuracyProved=False,
      physicalQualification=False,newQualifiedAreaM2=0,wholeSurfaceGoalComplete=False)
    save(out/'report.json',result);print({k:result[k] for k in ('passed','leadingPointInterior','allSeedBoundaryCounts','flowSpaceEvaluations','nativeLinearHistory','originalCapacityPassed')})


def main():
    if CLOSED:raise SystemExit('v249 audit closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(task.__file__),ROOT/'scripts/terminalCapacityV249.py',ROOT/'scripts/fasRecoveryFieldsV248.py',
      ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/su2PhaseAirFields.py',
      VOLUME,GEOMETRY,DOC/'PLAN.md',DOC/'NATIVE_CONTROL_PLAN.md',DOC/'AUDIT_RECOVERY_PLAN.md',
      BASE/'stencil01/report.json',BASE/'stencil01/actual_stencils.npz',BASE/'process_native/report.json',BASE/'process_native/terminal_capacity.json',
      task.PREV/'rocket01/state_000200.bin',ROOT/'outputs/v246_validation/initial01/boundary.csv']+[p for p in CASE.iterdir() if p.is_file()]
    return supervise(task.setup_guard(),BASE,'audit2',[sys.executable,str(Path(__file__).resolve()),'task'],90,1024,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v249 audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
