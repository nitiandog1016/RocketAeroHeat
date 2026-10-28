# -*- coding: utf-8 -*-
"""独立重建真实三维墙应力/Fourier、墙本构Jacobian、五量散度及实际Krylov更新。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1] if len(sys.argv)>1 else 'integration01';assert case.replace('_','').isalnum()
previous_case=sys.argv[2] if len(sys.argv)>2 else {'krylov02':'krylov01'}.get(case)
assert previous_case is None or (previous_case.replace('_','').isalnum() and previous_case!=case)
assert case in ('integration01','krylov01') or previous_case is not None, 'New continuation requires its explicit predecessor.'
target=VAL/(case+'_independent_audit.json');assert not target.exists()
file=VAL/case/('wall_velocity_integration_evidence.mat' if case=='integration01' else 'wall_velocity_krylov_evidence.mat')
domain_file=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in [file,domain_file]}
d=loadmat(domain_file,simplify_cells=True,variable_names=['domain'])['domain'];s=loadmat(file,simplify_cells=True);r=s.get('result',s)
Q=r['Q'];info=r['info'];gas=s['gas'];scale=r.get('residualScale',s.get('scale'));cv=gas['R']/(gas['gamma']-1);k=gas['mu']*gas['gamma']*cv/gas['Pr']
f=d['openFaces'];w=d['walls'];owner=w['ownerCell'].astype(int)-1;n=w['normal'];volume=d['cells']['volume'];P=info['primitive']
rho=Q[:,0];u=Q[:,1:4]/rho[:,None];T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
relative=lambda a,b:float(np.max(abs(a-b))/max(1,np.max(abs(b))))
assert relative(np.column_stack((rho,u,T)),P)<3e-15
assert np.all(rho>0) and np.all((T>=gas['temperatureRange'][0])&(T<=gas['temperatureRange'][1]))
dn=np.sum((w['centroid']-d['cells']['centroid'][owner])*n,axis=1);assert np.all(dn>0)
uw=np.asarray(s['boundary']['solidWall']['velocity']);assert np.all(uw==0) and s['boundary']['solidWall']['velocityGradient']=='dirichlet-tangent'
U=u[owner];Tw=info['wallTemperatureInput'];temperature=T[owner]
traction=-gas['mu']/dn[:,None]*(U+np.sum(U*n,axis=1)[:,None]*n/3)
pressure=(gas['gamma']-1)*(Q[owner,4]-.5*rho[owner]*np.sum(U*U,axis=1))
flux=np.column_stack((np.zeros(len(owner)),pressure[:,None]*n-traction,k*(temperature-Tw)/dn))
wall_error=relative(flux,info['wallFlux']);assert wall_error<3e-13
heat_error=relative(k*(temperature-Tw)/dn,info['wallHeatFluxIntoWall']);assert heat_error<3e-13
mechanical=-np.sum(U*traction,axis=1);entropy=w['area']*(mechanical/temperature+(1/temperature-1/Tw)*(-info['wallFlux'][:,4]))
assert np.all(mechanical>=0) and np.all(entropy>=-128*np.finfo(float).eps*max(1,np.max(abs(entropy))))
assert np.all(info['wallMechanicalPowerIntoFluid']==0)  # 静止壁外功，与邻壁代理黏性内能转换分开。
assert info['wallVelocityGradientAudit']['policy']=='dirichlet-tangent'

# 模块二：独立闭式静止墙本构Jacobian与有限差分，只是墙项，不是完整NS Jacobian。
selected=sorted(set([0,14280,16339,16613,int(np.argmin(dn)),int(np.argmax(np.linalg.norm(traction,axis=1)))]));jacobian_error=0.
for j in selected:
    q=Q[owner[j]];normal=n[j];dist=dn[j];vel=q[1:4]/q[0];temp=temperature[j]
    H=np.zeros((3,5));H[:,0]=-vel/q[0];H[:,1:4]=np.eye(3)/q[0]
    dp=(gas['gamma']-1)*np.r_[.5*np.dot(vel,vel),-vel,1.]
    dtemp=np.r_[.5*np.dot(vel,vel)-cv*temp,-vel,1.]/(q[0]*cv)
    D=np.zeros((5,5));D[1:4]=normal[:,None]*dp+gas['mu']/dist*(np.eye(3)+np.outer(normal,normal)/3)@H;D[4]=k/dist*dtemp
    def wall_flux(state):
        density=state[0];velocity=state[1:4]/density;internal=state[4]-.5*density*np.dot(velocity,velocity)
        thermal=internal/(density*cv);p=(gas['gamma']-1)*internal
        force=-gas['mu']/dist*(velocity+np.dot(velocity,normal)*normal/3)
        return np.r_[0.,p*normal-force,k*(thermal-Tw[j])/dist]
    finite=np.zeros((5,5));steps=2e-6*np.r_[q[0],np.repeat(q[0]*max(1,np.linalg.norm(vel)),3),max(abs(q[4]),q[0]*cv*temp)]
    for col in range(5):
        plus=q.copy();minus=q.copy();plus[col]+=steps[col];minus[col]-=steps[col]
        finite[:,col]=(wall_flux(plus)-wall_flux(minus))/(2*steps[col])
    jacobian_error=max(jacobian_error,relative(finite,D))
assert jacobian_error<3e-7

# 模块三：真实全部共享面及原墙散度，重新核对五量边界积分。
m=f['cellMinus'].astype(int)-1;p=f['cellPlus'].astype(int)-1;inside=p>=0;summed=np.zeros_like(Q)
integrated=f['area'][:,None]*info['faceFlux'];np.add.at(summed,m,-integrated);np.add.at(summed,p[inside],integrated[inside]);np.add.at(summed,owner,-w['area'][:,None]*info['wallFlux'])
divergence_error=relative(summed/volume[:,None],r['rate']);assert divergence_error<2e-12
outside=~inside;boundary=(f['area'][outside,None]*info['faceFlux'][outside]).astype(np.longdouble).sum(axis=0)+(w['area'][:,None]*info['wallFlux']).astype(np.longdouble).sum(axis=0)
denom=np.maximum(1,(f['area'][outside,None]*abs(info['faceFlux'][outside])).sum(axis=0)+(w['area'][:,None]*abs(info['wallFlux'])).sum(axis=0))
conservation=float(np.max(abs(summed.astype(np.longdouble).sum(axis=0)+boundary)/denom));assert conservation<1e-9
iteration_audit=None
if 'result' in s:
    integration_file=VAL/'integration01/wall_velocity_integration_evidence.mat'
    if previous_case is not None:
        previous_file=VAL/previous_case/'wall_velocity_krylov_evidence.mat';bindings[str(previous_file.relative_to(ROOT))]=sha(previous_file)
        previous_source=loadmat(previous_file,simplify_cells=True,variable_names=['result','displayHistory','report']);previous_result=previous_source['result']
        base=dict(Q=previous_result['Q'],rate=previous_result['rate']);previous_history=previous_source.get('displayHistory',previous_result['history'])
        prior_total_evaluations=previous_source['report'].get('totalNewEquationResidualEvaluations',previous_result['residualEvaluations'])
        assert set(s['displayHistory'])==set(previous_history)==set(r['history'])
        # 模块四：逐字段合并真实前驱的累计历史，不能只核对上一段局部三次。
        for key,old_values in previous_history.items():
            old_values=np.atleast_1d(old_values);tail=np.atleast_1d(r['history'][key])[1:]
            # MAT按值压缩的小整数可能读成uint8，先拓宽再累计，防止249+9绕回2。
            if key=='iteration':tail=tail.astype(np.int64)+int(old_values[-1])
            elif key=='residualEvaluations':tail=tail.astype(np.int64)+int(prior_total_evaluations)
            np.testing.assert_array_equal(s['displayHistory'][key],np.r_[old_values,tail])
        assert np.array_equal(s['displayHistory']['iteration'],np.arange(s['report']['totalNewEquationAcceptedIterations']+1))
        assert s['displayHistory']['residualEvaluations'][-1]<=s['report']['totalNewEquationResidualEvaluations']
        assert s['report']['priorAcceptedNonlinearIterations']==np.atleast_1d(previous_history['iteration'])[-1]
        assert s['report']['totalNewEquationResidualEvaluations']==prior_total_evaluations+r['residualEvaluations']
    else:
        bindings[str(integration_file.relative_to(ROOT))]=sha(integration_file)
        base=loadmat(integration_file,simplify_cells=True,variable_names=['Q','rate'])
    assert np.array_equal(base['Q'],r['initialQ']) and np.array_equal(base['rate'],r['initialRate']);del base
    assert r['options'].get('trialStateCoordinates','conservative')=='conservative'
    # 允许显式收紧线性目标，但不允许越过原0.2上限；物理/稳态门槛不变。
    assert 0<r['options']['linearRelativeTolerance']<=.2 and r['options']['steadyResidualTolerance']==2e-6
    assert r['options']['steadyFluxTolerance']==5e-10 and r['options']['spatialConservationTolerance']==1e-9
    np.testing.assert_array_equal(gas['temperatureRange'],[80,2000])
    trace=r['trace'];trace=[trace] if isinstance(trace,dict) else list(trace);assert len(trace)==r['acceptedNonlinearIterations']
    previous=r['initialQ'];previous_rate=r['initialRate'];update_error=0.;temperature_range=[]
    for item in trace:
        assert np.array_equal(previous,item['initialQ']) and np.array_equal(previous_rate,item['initialRate'])
        expected=item['initialQ']+item['alpha']*item['linearUpdateU']*scale
        update_error=max(update_error,relative(item['Q'],expected));assert update_error<3e-15
        iq=item['Q'];ir=iq[:,0];iv=iq[:,1:4]/ir[:,None];it=(iq[:,4]-.5*ir*np.sum(iv*iv,axis=1))/(ir*cv)
        assert np.all(ir>0) and np.all((it>=gas['temperatureRange'][0])&(it<=gas['temperatureRange'][1]));temperature_range += [float(it.min()),float(it.max())]
        before=float(np.sqrt(np.mean((item['initialRate']/scale)**2)));after=float(np.sqrt(np.mean((item['rate']/scale)**2)))
        assert after<=(1-r['options']['sufficientDecrease']*item['alpha'])*before
        assert item['linearReport']['fullSystemTargetMet'] and item['linearReport']['fullSystemRelativeDefect']<=r['options']['linearRelativeTolerance']
        previous=item['Q'];previous_rate=item['rate']
    assert np.array_equal(previous,r['Q']) and np.array_equal(previous_rate,r['rate'])
    j=r['firstJvpEvidence'];factor=2 if j['scheme']=='central' else 1
    expected=(j['ratePlus']-j['rateMinusOrBase'])/scale/(factor*j['actualScalarStep']);jvp_error=relative(expected,j['value']);assert jvp_error<3e-15
    maximum=float(np.max(abs(r['rate']/scale)));rms=float(np.sqrt(np.mean((r['rate']/scale)**2)))
    assert abs(maximum-r['maximumResidual'])/maximum<3e-15 and abs(rms-r['rmsMerit'])/rms<3e-15
    if r['steadyConverged']:assert maximum<=2e-6 and r['steadyBoundaryFluxBalance']<=5e-10
    assert not s['mapped']['eligibleForSurfaceCoverage'] and r['physicalTimeAdvanced']==0 and r['positivityClippingCount']==0
    iteration_audit=dict(acceptedIterations=len(trace),maximumResidual=maximum,rmsMerit=rms,updateArithmeticRelativeDifference=update_error,
                         firstFullResidualJvpArithmeticDifference=jvp_error,acceptedTemperatureRangeK=[min(temperature_range),max(temperature_range)] if temperature_range else None,
                         maximumKrylovIterations=int(r['options']['maximumKrylovIterations']),linearRelativeTolerance=float(r['options']['linearRelativeTolerance']),
                         originalLinearToleranceUpperBoundNotWeakened=True,
                         residualEvaluationsAfterLastLocalHistory=int(r['residualEvaluations']-np.atleast_1d(r['history']['residualEvaluations'])[-1]),
                         residualEvaluationsAfterLastAcceptedHistory=int(s['report'].get('totalNewEquationResidualEvaluations',r['residualEvaluations'])-np.atleast_1d(s.get('displayHistory',r['history'])['residualEvaluations'])[-1]),
                         originalSteadyAndTemperatureGuardrailsPreserved=True,
                         steadyConverged=bool(r['steadyConverged']),sourceHistoryAppended=False)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
report=dict(status='passed-independent-actual-static-wall-flux-constitutive-Jacobian-conservation-and-update-audit',case=case,processID=os.getpid(),
            predecessorCase=previous_case,wallCount=len(owner),wallFiveQuantityFluxRelativeDifference=wall_error,wallHeatRelativeDifference=heat_error,
            wallConstitutiveJacobianFDRelativeDifference=jacobian_error,constitutiveJacobianSampleWalls=len(selected),
            mechanicalCoolingWalls=int(np.count_nonzero(mechanical<0)),negativeWallEntropyFaces=int(np.count_nonzero(entropy<0)),
            fiveQuantityDivergenceRelativeDifference=divergence_error,independentSpatialConservationRelative=conservation,iterationAudit=iteration_audit,
            fullIndependentNSSolver=False,wholeOpenBoundaryEntropyBalanceProved=False,physicalSolution=False,newProductionSolvedAreaM2=0,
            inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))
