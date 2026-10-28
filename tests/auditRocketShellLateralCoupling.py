# -*- coding: utf-8 -*-
"""独立重建全部侧向近似作用，并用公开压力曲线抽查精确Euler方向；不是另一套NS。"""
from pathlib import Path
import ast,hashlib,json,math,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq,bisect

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];assert case.replace('_','').isalnum()
file=VAL/case/'lateral_coupling_evidence.mat';target=VAL/(case+'_independent_lateral_coupling.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
bindings={}
def read(path,names=None):
    bindings[str(path.relative_to(ROOT))]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=names)
s=read(file);source=Path(s['sourceFile']);d=read(source,['result','cut','gas']);r=d['result'];gas=d['gas']
if s['sourceKind']=='unit':cut=d['cut']
else:cut=read(ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat',['domain'])['domain']
Q=r['initialQ'];scale=np.asarray(r['residualScale']).reshape(5);rhs=r['initialRate']/scale;op=r['firstOperator']
log=r['linearLog'];log=log if isinstance(log,dict) else log[0];e=log['steadyDirectionEvidence'];direction=e['directionU'];original=e['value']
assert e['baseAcceptedIteration']==0 and r['linearEquation']=='steady-newton'
n=len(Q);nf=int(op['sourceFaces']);v=cut['cells']['volume'];f=cut['openFaces'];centers=cut['cells']['centroid']
m=f['cellMinus'].astype(int)-1;p=f['cellPlus'].astype(int)-1;side=np.flatnonzero((p>=0)&(f['kind']==3))
B=np.einsum('abn,nb->na',op['diagonal'],direction)
B[nf:]+=np.einsum('abn,nb->na',op['lower'][:,:,nf:],direction[:-nf])
B[:-nf]+=np.einsum('abn,nb->na',op['upper'][:,:,:-nf],direction[nf:])
line=direction/op['pseudoTimeStepSeconds']-B
cv=gas['R']/(gas['gamma']-1);conductivity=gas['mu']*gas['gamma']*cv/gas['Pr']
def primitive_direction(q,dq):
    rho=q[:,0];u=q[:,1:4]/rho[:,None];temp=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
    du=(dq[:,1:4]-u*dq[:,0,None])/rho[:,None]
    eint=dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+.5*np.sum(u*u,axis=1)*dq[:,0]
    dp=(gas['gamma']-1)*eint;dt=(eint-cv*temp*dq[:,0])/(rho*cv)
    return rho,u,temp,du,dp,dt
def rusanov(ql,qr,dl,dr,normal,distance):
    rl,ul,tl,dul,dpl,dtl=primitive_direction(ql,dl);rr,ur,tr,dur,dpr,dtr=primitive_direction(qr,dr)
    def euler(q,dq,u,t,du,dp):
        un=np.sum(u*normal,axis=1);dun=np.sum(du*normal,axis=1);pressure=q[:,0]*gas['R']*t
        return np.column_stack((dq[:,0]*un+q[:,0]*dun,
            dq[:,1:4]*un[:,None]+q[:,1:4]*dun[:,None]+dp[:,None]*normal,
            (dq[:,4]+dp)*un+(q[:,4]+pressure)*dun))
    wave=np.maximum(abs(np.sum(ul*normal,axis=1))+np.sqrt(gas['gamma']*gas['R']*tl),
                    abs(np.sum(ur*normal,axis=1))+np.sqrt(gas['gamma']*gas['R']*tr))
    out=.5*(euler(ql,dl,ul,tl,dul,dpl)+euler(qr,dr,ur,tr,dur,dpr))-.5*wave[:,None]*(dr-dl)
    velocity=dul-dur;traction=gas['mu']/distance[:,None]*(velocity+normal*np.sum(velocity*normal,axis=1)[:,None]/3)
    out[:,1:4]+=traction
    out[:,4]+=np.sum(.5*(ul+ur)*traction,axis=1)+conductivity/distance*(dtl-dtr)
    return out
off=np.zeros_like(Q);correction=np.zeros_like(Q)
def add_pair(out,mi,pi,ids,left,right):
    np.add.at(out,mi,-f['area'][ids,None]/v[mi,None]*left)
    np.add.at(out,pi,f['area'][ids,None]/v[pi,None]*right)
for first in range(0,len(side),8192):
    pick=np.arange(first,min(first+8192,len(side)));ids=side[pick];mi=m[ids];pi=p[ids];normal=f['normal'][ids]
    distance=np.sum((centers[pi]-centers[mi])*normal,axis=1)
    ql,qr=Q[mi],Q[pi];dl=direction[mi]*scale;dr=direction[pi]*scale;zero=np.zeros_like(dl)
    left=rusanov(ql,qr,dl,zero,normal,distance)/scale
    right=rusanov(ql,qr,zero,dr,normal,distance)/scale
    add_pair(off,mi,pi,ids,right,left)
    change=s['exactSideFluxDerivative'][pick]-left-right
    add_pair(correction,mi,pi,ids,change,change)
graph=line+off;exact_graph=graph+correction
def relative(actual,expected):return float(np.linalg.norm(actual-expected)/max(1,np.linalg.norm(expected)))
checks={key:relative(actual,s[key]) for key,actual in [('lineAction',line),('sideOffDiagonal',off),
    ('exactSideCorrection',correction),('rusanovGraphAction',graph),('exactSideGraphAction',exact_graph)]}
assert max(checks.values())<2e-11,checks
def metrics(value):
    delta=value-original
    return dict(relativeActionError=float(np.linalg.norm(delta)/np.linalg.norm(original)),
        equationRelativeActionErrors=(np.linalg.norm(delta,axis=0)/np.linalg.norm(original,axis=0)).tolist(),
        steadyLinearDefect=float(np.linalg.norm(-value-rhs)/np.linalg.norm(rhs)),
        cosineToOriginalAction=float(np.sum(value*original)/(np.linalg.norm(value)*np.linalg.norm(original))))
comparison={key:metrics(value) for key,value in [('line',line),('rusanovGraph',graph),('exactSideGraph',exact_graph)]}
for key,values in comparison.items():
    for field,value in values.items():np.testing.assert_allclose(value,s['comparison'][key][field],rtol=3e-10,atol=2e-12)

# 公开Clawpack压力曲线和独立Brent求根的抽样方向差分；不从MATLAB导数抄写公式。
public=ROOT/'outputs/v82_validation/clawpack_exact_euler.py';probe=ROOT/'outputs/v82_validation/probe_exact_riemann.py'
for path in (public,probe):bindings[str(path.relative_to(ROOT))]=sha(path)
assert sha(public)=='2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'
assert sha(probe)=='402ffeee97aed006c6784fa6b46d3276a3ab428897d77d07921308d9708c067a'
allowed={'pospart','sound_speed','beta','integral_curve_1','integral_curve_3','hugoniot_locus_1','hugoniot_locus_3'}
nodes=[node for node in ast.parse(public.read_text(encoding='utf8')).body if isinstance(node,ast.FunctionDef) and node.name in allowed]
reference={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(public),'exec'),reference)
nodes=[node for node in ast.parse(probe.read_text(encoding='utf8')).body if isinstance(node,ast.FunctionDef) and node.name=='exact']
namespace=dict(np=np,math=math,brentq=brentq,bisect=bisect,reference=reference)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(probe),'exec'),namespace)
def exact_euler(ql,qr,normal):
    pairs=[]
    for q in (ql,qr):
        rho=q[0];u=q[1:4]/rho;p=(gas['gamma']-1)*(q[4]-.5*rho*np.dot(u,u));pairs.append((rho,u,p))
    rl,ul,pl=pairs[0];rr,ur,pr=pairs[1];unl=float(ul@normal);unr=float(ur@normal)
    answer=namespace['exact']((rl,unl,pl),(rr,unr,pr),gas['gamma']);rho,un,p=answer['sample']
    u=(ul-unl*normal if answer['side']=='left' else ur-unr*normal)+un*normal
    energy=p/(gas['gamma']-1)+.5*rho*np.dot(u,u)
    return np.r_[rho*un,rho*u*un+p*normal,(energy+p)*un],answer['sampleRegion']
worst=0;sample_checks=[];branch_reductions=0
for pick in np.atleast_1d(s['sampleSideIndices']).astype(int)-1:
    face=side[pick];mi,pi=m[face],p[face];ql,qr=Q[mi],Q[pi];dl,dr=direction[mi]*scale,direction[pi]*scale;normal=f['normal'][face]
    distance=float((centers[pi]-centers[mi])@normal)
    base,region=exact_euler(ql,qr,normal)
    rate=0.
    for q,dq in ((ql,dl),(qr,dr)):
        rho=q[0];u=q[1:4]/rho;t=(q[4]-.5*rho*(u@u))/(rho*cv);speed=max(np.linalg.norm(u),math.sqrt(gas['gamma']*gas['R']*t))
        basis=np.r_[rho,np.repeat(rho*speed,3),rho*speed*speed];rate=max(rate,float(np.max(abs(dq)/basis)))
    assert rate>0
    h=2e-6/rate
    for attempt in range(7):
        plus,rp=exact_euler(ql+h*dl,qr+h*dr,normal);minus,rm=exact_euler(ql-h*dl,qr-h*dr,normal)
        if rp==rm==region:break
        h/=10;branch_reductions+=1
    assert rp==rm==region,'Sample crosses a nonsmooth branch; do not silently call it checked.'
    # 冻结黏性功速度；Newton应力和Fourier方向由独立本构表达恢复。
    zero=np.zeros((1,5));normal2=normal[None,:];distance2=np.array([distance]);ql2=ql[None,:];qr2=qr[None,:]
    _,ul,_,dul,_,dtl=primitive_direction(ql2,dl[None,:]);_,ur,_,dur,_,dtr=primitive_direction(qr2,dr[None,:])
    du=dul[0]-dur[0];traction=gas['mu']/distance*(du+normal*(du@normal)/3)
    actual=(plus-minus)/(2*h);actual[1:4]+=traction
    actual[4]+=.5*(ul[0]+ur[0])@traction+conductivity/distance*(dtl[0]-dtr[0]);actual/=scale
    expected=s['exactSideFluxDerivative'][pick];error=float(np.linalg.norm(actual-expected)/max(1,np.linalg.norm(expected)))
    worst=max(worst,error);sample_checks.append(dict(sideIndex=int(pick+1),faceID=int(face+1),relativeDifference=error))
assert worst<2e-5,worst
conservation=float(np.max(abs(np.sum(v[:,None]*correction,axis=0))/np.maximum(1,np.sum(abs(v[:,None]*correction),axis=0))))
assert conservation<1e-11 and s['report']['newOriginalNSEvaluations']==0
assert all(sha(ROOT/name)==digest for name,digest in bindings.items())
out=dict(status='passed-independent-all-lateral-actions-and-public-Riemann-direction-samples',case=case,processID=os.getpid(),
    cellCount=n,sideFaceCount=len(side),allActionRelativeDifferences=checks,comparison=comparison,
    exactDirectionSamples=sample_checks,maximumPublicExactDirectionRelativeDifference=worst,branchStepReductions=branch_reductions,
    sideReplacementVolumeConservationRelative=conservation,newOriginalNSEvaluations=0,newAcceptedIterations=0,
    fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'))
print(json.dumps({key:out[key] for key in ('status','case','comparison','maximumPublicExactDirectionRelativeDifference')}))
