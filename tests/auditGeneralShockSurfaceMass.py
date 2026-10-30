"""独立DOP853变分轨迹/投影流量、闭式四次轨迹、Brent与QUADPACK审核一般激波质量。"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.io import loadmat
from scipy.optimize import brentq

from auditShockSurfaceInlet import thermal,jump
from auditPressureDrivenSurfaceFlow import rows
from auditThermalPressureSurfaceFlow import sha
from auditParabolicShockLayer import geometric_check

ROOT=Path(__file__).resolve().parents[1]


def graph(shape,q):
    """模块一：直接微分图函数，不读取保存的梯度或曲率。"""
    q=np.asarray(q);K=np.asarray(shape['K']);B=float(shape.get('profileB',0.));g=K@q
    s=.5*q@g;D=1+2*B*s;assert D>0;root=np.sqrt(D)
    gradient=g/root;H=K/root-B*np.outer(g,g)/root**3
    if 'deformationWeights' in shape:
        w=np.atleast_1d(shape['deformationWeights']);directions=np.asarray(shape['deformationDirections']).reshape(2,-1)
        L=float(shape['deformationReferenceLengthM']);r=directions.T@q
        gradient+=directions@(w*r**3)/L**3
        H+=directions@np.diag(3*w*r*r/L**3)@directions.T
    return gradient,H


class ReferenceMass:
    """模块二：用完整2×2变分矩阵积分实际投影通量，不使用MATLAB网格或标量log面积。"""
    def __init__(self,shape,q):
        self.shape=shape;self.q=np.asarray(q);H0=graph(shape,np.zeros(2))[1]
        eigen=np.linalg.eigvalsh(H0);assert eigen.min()>0;self.end=25/eigen.min()
        maxrate=max(eigen.max(),np.linalg.eigvalsh(graph(shape,self.q)[1]).max());step=.07/maxrate
        def rhs(t,y):
            g,H=graph(shape,y[:2]);assert np.linalg.eigvalsh(H).min()>0
            M=y[2:].reshape(2,2)
            return np.r_[-g,(-H@M).ravel()]
        self.solution=solve_ivp(rhs,(0,self.end),np.r_[self.q,np.eye(2).ravel()],method='DOP853',
                                rtol=2e-12,atol=2e-14,max_step=step,dense_output=True)
        assert self.solution.success
        g0=graph(shape,self.q)[0];label=np.array([-g0[1],g0[0]])/np.linalg.norm(g0)
        def flux(t):
            y=self.solution.sol(t);g=graph(shape,y[:2])[0];v=y[2:].reshape(2,2)@label
            return np.linalg.det(np.column_stack([g,v]))/np.linalg.norm(g0)
        tail=np.linalg.det(self.solution.y[2:,-1].reshape(2,2))/np.trace(H0)
        self.coefficient=quad(flux,0,self.end,epsabs=2e-13,epsrel=3e-11)[0]+tail
        def stable(t,y):
            H=graph(shape,self.solution.sol(t)[:2])[1]
            return [np.trace(H)*y[0]-1,-1/y[0]]
        self.mass=solve_ivp(stable,(self.end,0),[1/np.trace(H0),0.],method='DOP853',
                            rtol=2e-12,atol=[2e-14,2e-12],max_step=step,dense_output=True)
        assert self.mass.success
        self.log0=self.mass.y[1,-1]
        assert abs(self.coefficient/self.mass.y[0,-1]-1)<2e-9

    def source(self,eta):
        if eta==0:return np.zeros(2)
        if eta==1:return self.q
        target=np.log(eta)
        t=brentq(lambda t:self.mass.sol(t)[1]-self.log0-target,0,self.end,xtol=2e-14,rtol=2e-13)
        return self.solution.sol(t)[:2]


def column_reference(column,up,gas):
    """模块三：独立法向跃迁、入口熵、二次基本形式及厚度积分。"""
    shape=column['shape'];q=column['coordinate'];frame=shape['frame'];g,H=graph(shape,q)
    ref=ReferenceMass(shape,q);U=np.linalg.norm(up['velocityVector']);mass=up['density']*U*ref.coefficient*np.linalg.norm(g)
    np.testing.assert_allclose(column['massPerShockWidthKgPerMS'],mass,rtol=2e-9,atol=1e-12)
    A=frame@np.vstack([g,np.eye(2)]);expected,velocity,n=jump(up,gas,A,1)
    tangent=up['velocityVector']-(up['velocityVector']@n)*n;tangent/=np.linalg.norm(tangent)
    transverse=np.cross(n,tangent);G=A.T@A;II=H/np.sqrt(1+g@g)
    operator=A@np.linalg.solve(G,II)@np.linalg.solve(G,A.T);basis=np.column_stack([tangent,transverse]);S=basis.T@operator@basis
    np.testing.assert_allclose(column['curvatureMatrix'],S,rtol=2e-11,atol=2e-10)
    ks,kb,twist=S[0,0],S[1,1],S[0,1];largest=np.linalg.eigvalsh(S).max()
    sine=up['velocityVector']@n/U;cosine=np.sqrt(1-sine*sine);us=U*cosine;vs=velocity@n
    P1=us*mass*ks;P2=.5*vs*mass*sine/cosine*(ks+kb);C=mass*(ks+kb)/(up['density']*U*cosine)
    H0=thermal(up['temperature'],up['pressure'],gas)[1]+U*U/2
    eigen=np.linalg.eigvalsh(graph(shape,np.zeros(2))[1]);power=eigen.sum()/eigen.min();cache={}
    def state(eta):
        if eta in cache:return cache[eta]
        source=ref.source(eta);gg=graph(shape,source)[0];AA=frame@np.vstack([gg,np.eye(2)])
        entry=jump(up,gas,AA,1)[0];entropy=entry[4]
        p=expected[0]+P1*(eta-1)-P2*(eta*eta-1);v=vs*(1+C*(eta-1))
        T=brentq(lambda T:thermal(T,p,gas)[2]-entropy,*gas['temperatureRangeK'],xtol=2e-11,rtol=1e-14)
        rho=p/(gas['R']*T);h=thermal(T,p,gas)[1];u2=2*(H0-h)-v*v;assert u2>0
        result=dict(pressure=p,temperature=T,density=rho,streamwiseSpeed=np.sqrt(u2),normalSpeed=v,entropy=entropy,source=source)
        cache[eta]=result;return result
    def capacity(n):
        return quad(lambda z:((1-z*ks)*(1-z*kb)-(z*twist)**2)/np.hypot(1-z*ks,z*twist),0,n,epsabs=1e-13,epsrel=2e-12)[0]
    state_error=0.;source_error=0.;distances=[]
    for j,eta in enumerate(np.atleast_1d(column['eta'])):
        value=state(float(eta))
        for name in ('pressure','temperature','density','streamwiseSpeed','normalSpeed','entropy'):
            actual=np.atleast_1d(column['state'][name])[j];state_error=max(state_error,abs(actual-value[name])/max(1,abs(value[name])))
        source_error=max(source_error,float(np.max(abs(column['state']['sourceShockCoordinate'][j]-value['source']))))
        def integrand(t):
            v=state(float(t**power));return power*t**(power-1)/(v['density']*v['streamwiseSpeed'])
        target=mass*quad(integrand,float(eta)**(1/power),1,epsabs=2e-12,epsrel=3e-10)[0]
        n=0. if eta==1 else brentq(lambda n:capacity(n)-target,0,(1-1e-8)/largest,xtol=2e-13,rtol=1e-13)
        distances.append(n)
    distance_error=float(np.max(abs(np.asarray(distances)-column['normalDistanceM'])))
    assert state_error<2e-9 and source_error<3e-10 and distance_error<1e-8
    assert not column['productionQualified'] and not column['boundaryLayerEdgeAvailable']
    return dict(maximumRelativeStateDifference=float(state_error),maximumSourceCoordinateDifferenceM=source_error,
                maximumDistanceDifferenceM=distance_error,massCoefficientM=float(ref.coefficient))


def separable_check(s):
    """模块四：可分离四次图的闭式轨迹和闭式流映射Jacobian独立审核。"""
    shape=s['separableShape'];q=s['separableMap']['coordinate'];K=shape['K'];k=np.asarray(np.diag(K),dtype=float)
    assert K[0,1]==K[1,0]==0 and float(shape['profileB'])==0
    w=shape['deformationWeights'];L=shape['deformationReferenceLengthM'];beta=w*q*q/(k*L**3)
    def trajectory(t):
        decay=np.exp(-k*t);return q*decay/np.sqrt(1+beta*(1-decay*decay))
    def W(t):
        decay=np.exp(-k*t);return np.prod(decay/(1+beta*(1-decay*decay))**1.5)
    a=quad(W,0,np.inf,epsabs=1e-13,epsrel=3e-11)[0]
    np.testing.assert_allclose(s['separableMap']['coefficientM'],a,rtol=2e-9,atol=1e-12)
    sources=[]
    for eta in s['separableEta']:
        if eta==0:sources.append(np.zeros(2));continue
        if eta==1:sources.append(q);continue
        t=brentq(lambda t:quad(W,t,np.inf,epsabs=1e-13,epsrel=3e-11)[0]/a-eta,0,30/k.min(),xtol=2e-14)
        sources.append(trajectory(t))
    err=float(np.max(abs(np.asarray(sources).T-s['separableSources'])));assert err<3e-10
    return dict(coefficientM=float(a),maximumSourceCoordinateDifferenceM=err)


def audit(evidence,process,target):
    assert not target.exists();s=loadmat(evidence,simplify_cells=True);p=json.loads(process.read_text(encoding='utf-8-sig'))
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    checks=[];separable=None
    if 'regression' in s:
        columns=[v['column'] for v in rows(s['regression'])]+[s['general'],s['fine']]
        for col in columns:checks.append(column_reference(col,s['up'],col['shock']['gas']))
        separable=separable_check(s)
        for bank in np.asarray(s['neighbors'],dtype=object):
            for entry in bank:
                if not isinstance(entry,dict):entry={name:getattr(entry,name) for name in entry._fieldnames}
                ref=ReferenceMass(s['shape'],entry['coordinate']);np.testing.assert_allclose(entry['coefficientM'],ref.coefficient,rtol=2e-9,atol=1e-12)
        assert abs(s['divergence'][-1]-1)<2e-6
    else:
        for group in rows(s['groups']):
            group['columns']=rows(group['columns'])
            for col in group['columns']:checks.append(column_reference(col,s['up'],s['gas']))
            geometric_check(s['geometry'],group)
        assert s['report']['newShapeFits']==0 and not s['report']['productionQualified']
    count=0;archived=[]
    revision_file=ROOT/'docs/validation/v99/raw/source_revision_bindings.json'
    versions=json.loads(revision_file.read_text(encoding='utf8'))['versions'] if revision_file.exists() else []
    for relative,digest in p['sourceSHA256'].items():
        if relative.startswith('outputs/'):continue
        file=ROOT/relative
        if not file.exists() or sha(file)!=digest:
            matches=[v for v in versions if v['restorePath']==relative and v['sha256']==digest]
            assert len(matches)==1 and sha(ROOT/matches[0]['path'])==digest,relative
            archived.append(matches[0])
        count+=1
    report=dict(status='passed-independent-general-graph-mass-variational-flux-source-and-layer-audit',
                evidenceSHA256=sha(evidence),processSHA256=sha(process),auditProgramSHA256=sha(Path(__file__)),columnsChecked=len(checks),
                maximumRelativeStateDifference=max(c['maximumRelativeStateDifference'] for c in checks),
                maximumSourceCoordinateDifferenceM=max(c['maximumSourceCoordinateDifferenceM'] for c in checks),
                maximumDistanceDifferenceM=max(c['maximumDistanceDifferenceM'] for c in checks),
                separableClosedFormCheck=separable,sourceSHARecordsVerified=count,sourceVersionsReadFromArchive=archived,productionQualified=False,
                experimentalValidation=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('evidence',type=Path);parser.add_argument('process',type=Path);parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();audit(args.evidence,args.process,args.report)
