# -*- coding: utf-8 -*-
"""独立组装WLS方向、Newton/Fourier/功链式项、共享散度和边界样本差分。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,kind=sys.argv[1:3]
assert case.replace('_','').isalnum() and kind in ('unit','real');target=VAL/(case+'_independent_wls_direction.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/('wls_direction_unit_evidence.mat' if kind=='unit' else 'wls_direction_real_evidence.mat');bindings={str(file.relative_to(ROOT)):sha(file)}
s=loadmat(file,simplify_cells=True);Q=s['Q'];gas=s['gas'];cut=s['cut'];base=s['baseInfo'];state=s['state'];scale=s['scale'];g=gas['gamma'];gm=g-1;cv=gas['R']/gm;kappa=gas['mu']*g*cv/gas['Pr']
c=cut['cells']['centroid'];v=cut['cells']['volume'];f=cut['openFaces'];w=cut['walls'];n=len(v);m=f['cellMinus'].astype(int)-1;p=f['cellPlus'].astype(int)-1;owner=w['ownerCell'].astype(int)-1
inside=np.flatnonzero((f['area']>0)&(m>=0)&(p>=0));outside=np.flatnonzero((f['area']>0)&((m>=0)^(p>=0)));active=np.flatnonzero(f['area']>0)
def rel(a,b):return float(np.max(abs(a-b))/max(1,np.max(abs(b))))
def rows(x):return [x] if isinstance(x,dict) else list(x)
def resolve(P,b,normal):
    rho=P[0];u=P[1:4];T=P[4];a=np.sqrt(g*gas['R']*T);un=float(u@normal);r=b['reservoirPrimitive'];ar=np.sqrt(g*gas['R']*r[4]);ur=float(r[1:4]@normal)
    if un>a:return P.copy(),'supersonic-outflow'
    if un<-a:return np.array(r),'supersonic-inflow'
    K=gas['R']*T/rho**gm;Kr=gas['R']*r[4]/r[0]**gm;Jp=un+2*a/gm;Jm=ur-2*ar/gm
    def acoustic():
        ub=(Jp+Jm)/2;ab=gm*(Jp-Jm)/4;selectedK,tangent=(Kr,r[1:4]-ur*normal) if ub<0 else (K,u-un*normal)
        return np.r_[(ab*ab/(g*selectedK))**(1/gm),tangent+ub*normal,ab*ab/(g*gas['R'])],('reservoir-backflow' if ub<0 else 'subsonic-outflow')
    if b['mode']=='riemannFarfield':return acoustic()
    pressure=float(b['pressure']);rb=(pressure/K)**(1/g);Tb=pressure/(rb*gas['R']);ab=np.sqrt(g*gas['R']*Tb);ub=Jp-2*ab/gm
    trace=np.r_[rb,u-un*normal+ub*normal,Tb]
    if un>=0 and ub>=0:return trace,'subsonic-outflow'
    astate,branch=acoustic()
    return (astate,branch) if branch=='reservoir-backflow' else (trace,'subsonic-outflow')
def force(G,normal):
    V=G[:,:3,:];stress=gas['mu']*(V+V.swapaxes(1,2)-(2/3)*np.trace(V,axis1=1,axis2=2)[:,None,None]*np.eye(3))
    return np.einsum('nij,nj->ni',stress,normal)
def divergence(face,wall):
    out=np.zeros((n,face.shape[1]));left=(m>=0)&(f['area']>0);right=(p>=0)&(f['area']>0);valid=owner>=0
    np.add.at(out,m[left],f['area'][left,None]/v[m[left],None]*face[left]);np.add.at(out,p[right],-f['area'][right,None]/v[p[right],None]*face[right])
    np.add.at(out,owner[valid],w['area'][valid,None]/v[owner[valid],None]*wall[valid]);return out
checks=[];one_sided=0
for ordinal,item in enumerate(rows(s['allEvidence'])):
    dq=item['directionU']*scale;rho=Q[:,0];u=Q[:,1:4]/rho[:,None];T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
    du=(dq[:,1:4]-u*dq[:,0,None])/rho[:,None];dp=gm*(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+.5*np.sum(u*u,axis=1)*dq[:,0])
    dT=dp/(rho*gas['R'])-T*dq[:,0]/rho;dphi=np.column_stack((du,dT));dGiven=item['boundaryDirection']['primitiveDerivative']
    # 直接几何正规方程，独立于MATLAB稀疏C和边界提升缓存。
    A=np.zeros((n,3,3));B=np.zeros((n,3,4))
    for first in range(0,len(inside),8192):
        ids=inside[first:first+8192];mi,pi=m[ids],p[ids];delta=c[pi]-c[mi];weight=f['area'][ids]/np.sum(delta*delta,axis=1)
        aa=weight[:,None,None]*delta[:,:,None]*delta[:,None,:];bb=weight[:,None,None]*delta[:,:,None]*(dphi[pi]-dphi[mi])[:,None,:]
        for cells in (mi,pi):np.add.at(A,cells,aa);np.add.at(B,cells,bb)
    incoming=outside[state['incoming'][outside].astype(bool)];cells=np.maximum(m[incoming],p[incoming]);delta=f['centroid'][incoming]-c[cells];weight=f['area'][incoming]/np.sum(delta*delta,axis=1)
    np.add.at(A,cells,weight[:,None,None]*delta[:,:,None]*delta[:,None,:]);np.add.at(B,cells,weight[:,None,None]*delta[:,:,None]*(dGiven[incoming,1:5]-dphi[cells])[:,None,:])
    valid=owner>=0;cells=owner[valid];delta=w['centroid'][valid]-c[cells];weight=w['area'][valid]/np.sum(delta*delta,axis=1)
    np.add.at(A,cells,weight[:,None,None]*delta[:,:,None]*delta[:,None,:]);np.add.at(B,cells,-weight[:,None,None]*delta[:,:,None]*dphi[cells,None,:])
    gradient=np.linalg.solve(A,B).swapaxes(1,2);vd=item['viscousDetail'];saved_gradient=np.asarray(vd['gradientDerivative']).transpose(2,0,1)
    gradient_error=rel(gradient,saved_gradient);assert gradient_error<2e-9,gradient_error
    viscous_flux=np.zeros((len(m),5));momentum_flux=np.zeros((len(m),3));heat=np.zeros(len(m));wg=heat.copy();wv=heat.copy()
    for first in range(0,len(active),8192):
        ids=active[first:first+8192];mi,pi=m[ids],p[ids];mask=(mi>=0)&(pi>=0);G=np.zeros((len(ids),4,3));velocity=np.zeros((len(ids),3));dvelocity=velocity.copy()
        ii=np.flatnonzero(mask);a,b=mi[ii],pi[ii];delta=c[b]-c[a];G[ii]=.5*(gradient[a]+gradient[b]);diff=dphi[b]-dphi[a]
        G[ii]+=(diff-np.einsum('nij,nj->ni',G[ii],delta))[:,:,None]*delta[:,None,:]/np.sum(delta*delta,axis=1)[:,None,None]
        velocity[ii]=.5*(u[a]+u[b]);dvelocity[ii]=.5*(du[a]+du[b]);oo=np.flatnonzero(~mask);ff=ids[oo];cells=np.maximum(mi[oo],pi[oo]);G[oo]=gradient[cells]
        inflow=state['incoming'][ff].astype(bool);ii=oo[inflow];faces=ff[inflow];selected=cells[inflow];delta=f['centroid'][faces]-c[selected]
        G[ii]+=(dGiven[faces,1:5]-dphi[selected]-np.einsum('nij,nj->ni',G[ii],delta))[:,:,None]*delta[:,None,:]/np.sum(delta*delta,axis=1)[:,None,None]
        ii=oo[~inflow];normal=f['normal'][ff[~inflow]];G[ii]-=np.einsum('nij,nj->ni',G[ii],normal)[:,:,None]*normal[:,None,:]
        velocity[oo]=state['primitive'][ff,1:4];dvelocity[oo]=dGiven[ff,1:4];normal=f['normal'][ids];traction=force(G,normal)
        momentum_flux[ids]=traction;heat[ids]=kappa*np.sum(G[:,3,:]*normal,axis=1);wg[ids]=np.sum(velocity*traction,axis=1);wv[ids]=np.sum(dvelocity*base['viscousFaceFlux'][ids,1:4],axis=1)
        viscous_flux[ids,1:4]=traction;viscous_flux[ids,4]=heat[ids]+wg[ids]+wv[ids]
    wall_flux=np.zeros((len(owner),5));cells=owner[valid];normal=w['normal'][valid];dn=np.sum((w['centroid'][valid]-c[cells])*normal,axis=1)
    G=-dphi[cells,:,None]*normal[:,None,:]/dn[:,None,None];wall_flux[valid,1:4]=force(G,normal);wall_flux[valid,4]=kappa*np.sum(G[:,3,:]*normal,axis=1)
    flux_error=rel(viscous_flux,vd['faceViscousFluxDerivative']);wall_error=rel(wall_flux,vd['wallViscousFluxDerivative']);viscous=divergence(viscous_flux,wall_flux)/scale
    assert max(flux_error,wall_error,rel(viscous,item['viscous']))<2e-9
    work_velocity=np.zeros((n,5));work_velocity[:,4]=divergence(wv[:,None],np.zeros((len(owner),1)))[:,0]/scale[4]
    np.testing.assert_allclose(work_velocity,vd['workVelocityU'],rtol=3e-10,atol=3e-12)
    ed=item['eulerDetail'];euler=-divergence(ed['faceEulerFluxDerivative'],ed['wallEulerFluxDerivative'])/scale
    assert rel(euler,item['euler'])<2e-10;analytic=euler+viscous;fd_checks=[]
    for refinement in rows(item['refinements']):
        h=float(refinement['h']);full=(refinement['ratePlus']-refinement['rateMinus'])/(2*h)/scale
        fv=(refinement['faceViscousPlus']-refinement['faceViscousMinus'])/(2*h)
        fc=(refinement['faceEulerPlus']-refinement['faceEulerMinus'])/(2*h)
        fd_checks.append(dict(fullNSRelativeDifference=rel(analytic,full),viscousFluxRelativeDifference=rel(viscous_flux,fv),eulerFluxRelativeDifference=rel(ed['faceEulerFluxDerivative'],fc)))
        assert max(fd_checks[-1].values())<(2e-6 if kind=='unit' else 1e-3)
    # 用直接幂形式特征关系作开放样本差分，和原实现exp/log算术分开。
    sample=np.unique(np.r_[np.linspace(0,len(outside)-1,min(65,len(outside))).astype(int)]);boundary_errors=[]
    for index in sample:
        face=int(outside[index]);cell=int(max(m[face],p[face]));normal=f['normal'][face]*(1 if m[face]>=0 else -1);side=1 if m[face]>=0 else 0
        name=[['xmin','xmax'],['ymin','ymax'],['zmin','zmax']][int(f['axis'][face])-1][side];b=s['boundary'][name]
        if b['type']!='characteristic':continue
        P=base['primitive'][cell];dP=np.r_[dq[cell,0],du[cell],dT[cell]];ref,branch=resolve(P,b,normal)
        basis=np.r_[P[0],np.repeat(max(np.linalg.norm(P[1:4]),np.sqrt(g*gas['R']*P[4])),3),P[4]];rate=np.max(abs(dP)/basis)
        if rate==0:actual=np.zeros(5)
        else:
            h=2e-6/rate;plus,bp=resolve(P+h*dP,b,normal);minus,bm=resolve(P-h*dP,b,normal)
            if bp==bm==branch:actual=(plus-minus)/(2*h)
            elif bp==branch:actual=(plus-ref)/h;one_sided+=1
            elif bm==branch:actual=(ref-minus)/h;one_sided+=1
            else:raise AssertionError('Independent boundary samples cross both selected branches.')
        boundary_errors.append(rel(actual,dGiven[face]));assert boundary_errors[-1]<2e-4
    # 共享面导数符号/五量边界账；不是物理稳态证明。
    full_flux=ed['faceEulerFluxDerivative']-viscous_flux;full_wall=ed['wallEulerFluxDerivative']-wall_flux
    signs=np.where(m[outside]>=0,1,-1);net=np.sum(full_flux[outside].astype(np.longdouble)*(signs*f['area'][outside])[:,None],axis=0)+np.sum(full_wall.astype(np.longdouble)*w['area'][:,None],axis=0)
    volume=np.sum(analytic.astype(np.longdouble)*scale*v[:,None],axis=0);denom=np.maximum(1,np.sum(abs(full_flux[outside])*f['area'][outside,None],axis=0)+np.sum(abs(full_wall)*w['area'][:,None],axis=0))
    balance=float(np.max(abs(volume+net)/denom));assert balance<1e-10
    checks.append(dict(direction=ordinal+1,gradientRelativeDifference=gradient_error,viscousFluxRelativeDifference=flux_error,wallViscousFluxRelativeDifference=wall_error,
        boundarySampleCount=len(boundary_errors),maximumBoundarySampleRelativeDifference=max(boundary_errors,default=0),independentDerivativeBalance=balance,originalNSFiniteDifferences=fd_checks,
        analyticL2RelativeDifference=float(np.linalg.norm(analytic-full)/np.linalg.norm(full))))
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-WLS-constitutive-chain-boundary-samples-and-complete-direction-audit',case=case,kind=kind,processID=os.getpid(),directions=checks,
    oneSidedSelectedBoundaryChecks=one_sided,fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
