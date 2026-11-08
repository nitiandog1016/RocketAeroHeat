"""模块职责：独立解析SA-Edwards源项链式导数与六式散度组装；不是生产Jacobian。"""
from pathlib import Path
import json,time
import numpy as np
from readEilmerNative import read_grid,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
from auditFrozenMappedGeometryV213 import FACES
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v214_validation'

def sa_source(rho,drho,mu,dmu,nh,dnh,A,dA,b,db,d):
    # 模块一：按Edwards/Chandra变体的张量表达式和显式链式法则，独立计算值/切向。
    k=.41;cb1=.1355;cb2=.622;sigma=2/3;cw1=cb1/k/k+(1+cb2)/sigma
    nu=mu/rho;dnu=dmu/rho-mu*drho/rho**2
    chi=nh/nu;dc=dnh/nu-nh*dnu/nu**2
    chi3=chi**3;c=7.1**3;fv1=chi3/(chi3+c);dfv1=3*chi**2*c/(chi3+c)**2*dc
    tr=np.trace(A,axis1=1,axis2=2);dtr=np.trace(dA,axis1=1,axis2=2)
    S=np.sum((A+A.transpose(0,2,1))*A,axis=(1,2))-2/3*tr**2
    dS=np.sum((dA+dA.transpose(0,2,1))*A+(A+A.transpose(0,2,1))*dA,axis=(1,2))-4/3*tr*dtr
    assert np.all(S>0)
    sh=np.sqrt(S)*(nu+fv1*nh)
    dsh=.5*dS/np.sqrt(S)*(nu+fv1*nh)+np.sqrt(S)*(dnu+dfv1*nh+fv1*dnh)
    ell=nh**2/(k*k*d*d);dell=2*nh*dnh/(k*k*d*d)
    capped=sh*10<=ell
    normal=np.where(capped,10,ell/sh);dn=np.where(capped,0,dell/sh-ell*dsh/sh**2)
    r=np.tanh(normal)/np.tanh(1);dr=(1-np.tanh(normal)**2)*dn/np.tanh(1)
    g=.7*r+.3*r**6;dg=(.7+1.8*r**5)*dr
    factor=(65/(g**6+64))**(1/6);fw=g*factor;dfw=factor*dg*64/(g**6+64)
    prod=rho*cb1*sh;dprod=cb1*(drho*sh+rho*dsh)
    destruction=rho*cw1*fw*nh**2/d**2
    ddes=cw1/d**2*(drho*fw*nh**2+rho*dfw*nh**2+2*rho*fw*nh*dnh)
    norm=np.sum(b*b,axis=1);diss=cb2/sigma*rho*norm
    ddis=cb2/sigma*(drho*norm+2*rho*np.sum(b*db,axis=1))
    return prod-destruction+diss,dprod-ddes+ddis,capped,S

def main():
    start=time.monotonic();case=BASE/'control01';out=BASE/'assembly01';out.mkdir(exist_ok=False)
    source=ROOT/'outputs/v211_validation/source01/source/lmr/turbulence.d'
    build=json.loads((BASE/'build01/report.json').read_text());assert sha(source)==build['inputSHA256']['/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/turbulence.d']
    grid=read_grid(case/'lmrsim/grid/grid-0000');a=np.genfromtxt(str(case/'nk_samples.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    state=np.genfromtxt(str(case/'directions.csv'),delimiter=',',names=True)
    geometry=np.genfromtxt(str(ROOT/'outputs/v213_validation/control03/geometry.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    c=geometry[geometry['kind']=='cell'];f=geometry[geometry['kind']=='face'];vol=c['volumeOrArea'];d=c['wallDistance']
    lookup={tuple(sorted(row[row>=0])):i for i,row in enumerate(grid['faces'])}
    ids=np.array([[lookup[tuple(sorted(row[list(local)]))] for local in FACES] for row in grid['cells'][:,:6]])
    xyz=lambda x:np.stack([x[k] for k in ['x','y','z']],axis=1)
    normals=np.stack([f[k] for k in ['nx','ny','nz']],axis=1)
    signs=np.sign(np.einsum('ijk,ijk->ij',xyz(f)[ids]-xyz(c)[:,None,:],normals[ids]));assert np.all(np.abs(signs)==1)
    for boundary in grid['boundaries'].values():
        for fid,s,owner in zip(boundary['ids'],boundary['signs'],boundary['owners']):assert signs[owner][ids[owner]==fid][0]==s
    gas=PhaseAwareViscousAir();records=[]
    for regime in range(3):
        r=a[a['regime']==regime];q=state[state['regime']==regime]
        get=lambda kind,col,width:r[r['kind']==kind][col].reshape(-1,width)
        A=get('velocityGradient','frozen',9).reshape(-1,3,3);dA=get('velocityGradient','frozenDerivative',9).reshape(-1,3,3)
        b=get('saGradient','frozen',3);db=get('saGradient','frozenDerivative',3)
        mu=gas.interpolated_transport(q['p'],q['T'])['dynamicViscosityPas']
        ii=np.searchsorted(gas.T,q['T'],side='right')-1;dmu=np.diff(gas.mu)[ii]/np.diff(gas.T)[ii]*.2
        source0,sourceJ,capped,S=sa_source(q['rho'],q['drho'],mu,dmu,q['nuhat'],np.full(416,1e-7),A,dA,b,db,d)
        flux0=get('integratedFlux','frozen',6);fluxJ=get('integratedFlux','frozenDerivative',6)
        residual0=get('residual','frozen',6);residualJ=get('residual','frozenDerivative',6)
        divergence=lambda flux:-np.sum(flux[ids]*signs[:,:,None],axis=1)/vol[:,None]
        expected0=divergence(flux0);expectedJ=divergence(fluxJ)
        expected0[:,5]+=source0;expectedJ[:,5]+=sourceJ
        err0=np.abs(expected0-residual0)/np.maximum(1,np.maximum(np.abs(expected0),np.abs(residual0)))
        errJ=np.abs(expectedJ-residualJ)/np.maximum(1,np.maximum(np.abs(expectedJ),np.abs(residualJ)))
        # 模块二：保存差分缺口经过原有向面/体积算子后的量值，不将放大直接归因舍入。
        fdflux=(get('integratedFlux','frozenPlus',6)-get('integratedFlux','frozenMinus',6))/2e-4
        fdres=(get('residual','frozenPlus',6)-get('residual','frozenMinus',6))/2e-4
        propagated=divergence(fdflux-fluxJ)[:,:5];actual=(fdres-residualJ)[:,:5]
        discrepancy=np.abs(propagated-actual)
        records.append(dict(regime=regime,valueMax=float(err0.max()),assemblyDerivativeMax=float(errJ.max()),
            derivativeByEquation=errJ.max(axis=0).tolist(),cappedCells=int(capped.sum()),minStrainInvariant=float(S.min()),
            firstFiveFDGapMax=float(np.abs(actual).max()),firstFiveFluxPropagationMax=float(np.abs(propagated).max()),
            propagationRemainderMax=float(discrepancy.max()),
            derivativeGatePassed=bool(errJ.max()<=1e-5),valueGatePassed=bool(err0.max()<=1e-5)))
    result=dict(auditCompleted=True,records=records,originalFiniteDifferenceAcceptance=False,
        independentScope='analytic SA source chain rule and FV assembly; flux/gradient tangents remain native inputs',
        sourceSHA256=sha(source),inputSHA256=sha(case/'nk_samples.csv'),seconds=time.monotonic()-start,
        nativeRuns=0,newQualifiedAreaM2=0)
    (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
