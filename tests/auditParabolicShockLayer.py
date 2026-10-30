# -*- coding: utf-8 -*-
"""独立Brent/QUADPACK/曲面二次基本形式，审核共同激波层输运及真实体形残差。"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.integrate import quad
from scipy.io import loadmat
from scipy.linalg import expm
from scipy.optimize import brentq

from auditShockSurfaceInlet import thermal, jump
from auditThermalPressureSurfaceFlow import sha
from auditPressureDrivenSurfaceFlow import rows

ROOT=Path(__file__).resolve().parents[1]


def reference(column, upstream, gas):
    """模块一：不用保存的压力系数、熵插值、Gauss质量积或曲率矩阵作为参考。"""
    shape=column['shape'];K=np.asarray(shape['K']);frame=np.asarray(shape['frame']);q=np.asarray(column['coordinate'])
    U=np.linalg.norm(upstream['velocityVector']);rho1=upstream['density'];g=K@q
    A=frame@np.vstack([g,np.eye(2)]);n=frame@np.r_[1.,-g];n/=np.linalg.norm(n)
    expected, shock_velocity,_=jump(upstream,gas,A,1)
    upstream_velocity=np.asarray(upstream['velocityVector']);t=upstream_velocity-(upstream_velocity@n)*n;t/=np.linalg.norm(t)
    b=np.cross(n,t);G=A.T@A;II=K/np.sqrt(1+g@g)
    shape_operator=A@np.linalg.solve(G,II)@np.linalg.solve(G,A.T)
    S=np.column_stack([t,b]).T@shape_operator@np.column_stack([t,b])
    np.testing.assert_allclose(column['curvatureMatrix'],S,rtol=2e-12,atol=2e-11)
    ks,kb,twist=S[0,0],S[1,1],S[0,1];largest=np.linalg.eigvalsh(S).max()
    sine=upstream_velocity@n/U;cosine=np.sqrt(1-sine*sine);vs=shock_velocity@n;us=U*cosine
    mass=rho1*U*np.linalg.norm(g)/np.trace(K)
    # 独立积分实际激波线的上游投影质量通量，不把闭式m再次代入作证明。
    label_tangent=np.array([-q[1],q[0]])
    def projected_mass(tau):
        E=expm(K*tau);point=E@q;label=E@label_tangent
        return rho1*U*np.linalg.det(np.column_stack([K@point,label]))
    psi=quad(projected_mass,-np.inf,0,epsabs=1e-13,epsrel=2e-11)[0]
    label_at_fixed_x=label_tangent-g*np.dot(g,label_tangent)/(g@g)
    independent_mass=psi/np.linalg.norm(A@label_at_fixed_x)
    np.testing.assert_allclose(mass,independent_mass,rtol=2e-10,atol=1e-12)
    P1=us*mass*ks;P2=.5*vs*mass*sine/cosine*(ks+kb)
    C=mass/(rho1*U*cosine)*(ks+kb)
    H0=thermal(upstream['temperature'],upstream['pressure'],gas)[1]+U*U/2
    power=np.trace(K)/np.linalg.eigvalsh(K).min();cache={}

    def state(eta):
        key=float(eta)
        if key in cache:return cache[key]
        source=np.zeros(2) if eta==0 else expm(K*(np.log(eta)/np.trace(K)))@q
        source_A=frame@np.vstack([K@source,np.eye(2)])
        entry=jump(upstream,gas,source_A,1)[0]
        p=expected[0]+P1*(eta-1)-P2*(eta*eta-1);v=vs*(1+C*(eta-1));entropy=entry[4]
        T=brentq(lambda T:thermal(T,p,gas)[2]-entropy,*gas['temperatureRangeK'],xtol=2e-11,rtol=1e-14)
        rho=p/(gas['R']*T);h=thermal(T,p,gas)[1];u2=2*(H0-h)-v*v;assert u2>0
        value=dict(pressure=p,temperature=T,density=rho,streamwiseSpeed=np.sqrt(u2),normalSpeed=v,
                   entropy=entropy,source=source);cache[key]=value;return value

    def capacity(n):
        return quad(lambda z:((1-z*ks)*(1-z*kb)-(z*twist)**2)/
                    np.hypot(1-z*ks,z*twist),0,n,epsabs=1e-13,epsrel=2e-12)[0]

    distances=[];state_error=0.;entropy_error=0.
    for j,eta in enumerate(np.atleast_1d(column['eta'])):
        local=state(float(eta))
        def integrand(tau):
            s=state(tau**power)
            return power*tau**(power-1)/(s['density']*s['streamwiseSpeed'])
        mass_distance=mass*quad(integrand,float(eta)**(1/power),1,epsabs=2e-12,epsrel=3e-10)[0]
        n=0. if eta==1 else brentq(lambda n:capacity(n)-mass_distance,0,(1-1e-8)/largest,xtol=2e-13,rtol=1e-13)
        distances.append(n)
        for name in ['pressure','temperature','density','streamwiseSpeed','normalSpeed','entropy']:
            actual=np.atleast_1d(column['state'][name])[j];expected_value=local[name]
            error=abs(actual-expected_value)/max(1,abs(expected_value));state_error=max(state_error,float(error))
        entropy_error=max(entropy_error,float(abs(column['state']['entropy'][j]-local['entropy'])))
        np.testing.assert_allclose(column['state']['sourceShockCoordinate'][j],local['source'],rtol=2e-12,atol=2e-14)
        stored_n=column['normalDistanceM'][j]
        J=(1-stored_n*ks)*(1-stored_n*kb)-(stored_n*twist)**2
        np.testing.assert_allclose(column['jacobianRatio'][j],J,rtol=2e-12,atol=2e-14)
        cp,h,Svalue=thermal(column['state']['temperature'][j],column['state']['pressure'][j],gas)
        u=column['state']['streamwiseSpeed'][j];v=column['state']['normalSpeed'][j]
        np.testing.assert_allclose(h+(u*u+v*v)/2,H0,rtol=3e-13,atol=2e-7)
    distance_error=float(np.max(abs(np.asarray(distances)-column['normalDistanceM'])))
    assert state_error<2e-9 and distance_error<1e-8
    assert not column['bodyMatched'] and not column['boundaryLayerEdgeAvailable'] and not column['productionQualified']
    return dict(maximumRelativeStateDifference=state_error,maximumDistanceDifferenceM=distance_error,
                maximumEntropyDifferenceJPerKgK=entropy_error,wallDistanceM=float(distances[-1]))


def geometric_check(geometry,result):
    """模块二：原面平面方程与重心解，逐列复核实际体形残差，不重做拟合。"""
    V=geometry['V'];F=geometry['F'].astype(int)-1;triangles=V[F]
    a=triangles[:,1]-triangles[:,0];b=triangles[:,2]-triangles[:,0];N=np.cross(a,b)
    aa=np.einsum('ij,ij->i',a,a);bb=np.einsum('ij,ij->i',b,b);ab=np.einsum('ij,ij->i',a,b)
    determinant=aa*bb-ab*ab;m=result['geometricMapping'];values=[]
    for j,(origin,direction) in enumerate(zip(m['rayPosition'],m['rayDirection'])):
        divisor=N@direction
        distance=np.divide(np.einsum('ij,ij->i',triangles[:,0]-origin,N),divisor,
                           out=np.full(len(F),np.nan),where=abs(divisor)>64*np.finfo(float).eps*np.linalg.norm(N,axis=1))
        p=origin+distance[:,None]*direction;d=p-triangles[:,0]
        da=np.einsum('ij,ij->i',d,a);db=np.einsum('ij,ij->i',d,b)
        u=(da*bb-db*ab)/determinant;w=(db*aa-da*ab)/determinant;tol=m['barycentricTolerance']
        valid=(distance>0)&(u>=-tol)&(w>=-tol)&(u+w<=1+tol)
        assert valid.any();nearest=distance[valid].min();ties=np.flatnonzero(valid&(abs(distance-nearest)<=m['distanceToleranceM']))
        index=ties.min();assert index+1==m['originalFaceID'][j]
        assert abs(distance[index]-m['distanceM'][j])<1e-11
        predicted=result['columns'][j]['normalDistanceM'][-1]
        residual=predicted-distance[index];assert abs(residual-result['residualM'][j])<1e-11;values.append(residual)
    return np.asarray(values)


def audit(evidence,process,target):
    assert not target.exists();s=loadmat(evidence,simplify_cells=True)
    p=json.loads(process.read_text(encoding='utf-8-sig'))
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    checks=[]
    for case in rows(s['cases']):
        c=case['column'];checks.append(reference(c,
            {'pressure':s['fs']['pressure'],'temperature':s['fs']['temperature'],'density':s['fs']['density'],
             'velocityVector':s['fs']['flowVelocityVector']},c['shock']['gas']))
        a=case['sameStreamlineUpstream'];b=case['sameStreamlineDownstream']
        np.testing.assert_allclose(a['state']['sourceShockCoordinate'][1],b['state']['sourceShockCoordinate'][1],rtol=0,atol=1e-14)
        np.testing.assert_allclose(a['state']['entropy'][1],b['state']['entropy'][1],rtol=0,atol=1e-9)
    fit=s['bodyFit'];up={'pressure':s['fs']['pressure'],'temperature':s['fs']['temperature'],
                       'density':s['fs']['density'],'velocityVector':s['fs']['flowVelocityVector']}
    residuals=[]
    for key in ['initial','result']:
        result=fit[key];result['columns']=rows(result['columns'])
        for column in result['columns']:checks.append(reference(column,up,column['shock']['gas']))
        residuals.append(geometric_check(s['geometry'],result))
    assert np.linalg.norm(residuals[1])<np.linalg.norm(residuals[0])
    matched=bool(np.max(abs(residuals[1]))<=fit['bodyMatchToleranceM'])
    assert bool(fit['shapeMatchingConverged'])==matched and not fit['productionQualified']
    source_count=0
    for rel,digest in p['sourceSHA256'].items():
        if rel.startswith('outputs/'):continue
        file=ROOT/rel;assert file.exists() and sha(file)==digest,rel;source_count+=1
    report=dict(status='passed-independent-physical-entropy-QUADPACK-thickness-and-real-body-residual-audit',
                evidenceSHA256=sha(evidence),processSHA256=sha(process),columnsChecked=len(checks),
                maximumRelativeStateDifference=max(x['maximumRelativeStateDifference'] for x in checks),
                maximumDistanceDifferenceM=max(x['maximumDistanceDifferenceM'] for x in checks),
                initialBodyRMSM=float(np.linalg.norm(residuals[0])/3),finalBodyRMSM=float(np.linalg.norm(residuals[1])/3),
                maximumBodyResidualM=float(np.max(abs(residuals[1]))),bodyMatchToleranceM=float(fit['bodyMatchToleranceM']),
                sampledBodyMatchingConverged=matched,sourceSHARecordsVerified=source_count,
                stagnationGlobalMatchingSolved=False,boundaryLayerEdgeAvailable=False,experimentalValidation=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('process',type=Path)
    p.add_argument('--report',type=Path,required=True);a=p.parse_args();audit(a.evidence,a.process,a.report)
