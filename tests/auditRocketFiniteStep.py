"""独立公开Riemann参考、保存逐步更新和凸选择KKT核对；不替代整箭验证。"""
from pathlib import Path
import json, sys, hashlib
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from riemannReservoirReference import exact, PUBLIC_BINDINGS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v84_validation'
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def prim(q,g):
    rho=q[:,0];v=q[:,1:4]/rho[:,None]
    t=(q[:,4]/rho-.5*np.sum(v*v,axis=1))/(g['R']/(g['gamma']-1))
    return np.column_stack((rho,v,t))

def main():
    p=loadmat(OUT/'diagnosis01/finite_step_patch.mat',simplify_cells=True)['patch']
    m=p['fluxMinus'].astype(int)-1;n=p['fluxPlus'].astype(int)-1
    maxerr=0.;checks=0
    original=loadmat(ROOT/'outputs/v83_validation/krylov24/wall_velocity_krylov_evidence.mat',simplify_cells=True,variable_names=['result'])['result']
    # 与独立公开压力曲线+Brent根比较两个真实端点的全部276面。
    for name in ['Q0','Q1']:
        v=prim(p[name],p['gas'])
        expected=original['initialInfo' if name=='Q0' else 'info']['convectiveFaceFlux'][p['fluxGlobalFaceID'].astype(int)-1]
        for face,(a,b,normal) in enumerate(zip(m,n,p['fluxNormal'])):
            state,flux,region=exact(v[a],v[b],normal,p['gas'])
            assert np.all(np.isfinite(flux)) and state[0]>0 and 80<=state[4]<=2000
            error=np.linalg.norm((flux-expected[face])/p['scale'])/max(np.linalg.norm(expected[face]/p['scale']),1.)
            assert error<1e-10
            maxerr=max(maxerr,error)
            checks+=1
    case=sys.argv[1] if len(sys.argv)>1 else 'patch02';assert case in ['patch02','patch03']
    s=loadmat(OUT/case/'patch_comparison.mat',simplify_cells=True)
    stats={}
    for mode in ['baseline','guarded']:
        run=s[mode];h=run['summary']['history'];trials=run['trials'];dirs=run['directions']
        assert len(h)==7 and len(dirs)==6
        for i,d in enumerate(dirs):
            q=d['Q'];v=prim(q,p['gas']);assert np.all(v[:,0]>0) and np.min(v[:,4])>=80 and np.max(v[:,4])<=2000
            accepted=[t for t in trials if t['iteration']==i+1 and t['accepted']]
            assert len(accepted)==1
            t=accepted[0];expect=q+t['alpha']*d['direction']*p['scale']
            actual=dirs[i+1]['Q'] if i+1<len(dirs) else run['Q']
            assert np.allclose(expect,actual,rtol=2e-15,atol=0)
            assert d['nsDefect']<=.1 and h[i+1]['rms']<=(1-1e-4*t['alpha'])*h[i]['rms']
            assert h[i+1]['spatialConservationRelative']<=1e-9
        stats[mode]={'initialRMS':h[0]['rms'],'finalRMS':h[-1]['rms'],'initialMaximum':h[0]['maximum'],'finalMaximum':h[-1]['maximum'],
            'seconds':run['summary']['elapsedSeconds'],'NSCalls':run['summary']['residualEvaluations'],
            'allSixMaximumNonIncreasing':all(h[i+1]['maximum']<=h[i]['maximum'] for i in range(6)),
            'lastFourRMSReductionFraction':1-h[-1]['rms']/h[2]['rms'],'steadyConverged':False}
    u=loadmat(OUT/'unit01/minimum_state_unit.mat',simplify_cells=True)
    A=u['A'];b=u['b'];z=u['z'];r=u['r'];U,sv,Vt=np.linalg.svd(A,full_matrices=False);c=U.T@b
    perp=np.linalg.norm(b-U@c)**2;den=np.linalg.norm(b)
    fun=lambda l: np.sqrt(perp+np.sum((l/(sv*sv+l)*c)**2))/den-.0999
    lam=brentq(fun,0,max(sv)**2,xtol=1e-12)
    reference=Vt.T@((sv/(sv*sv+lam))*c);err=np.linalg.norm(reference-z)/np.linalg.norm(z);assert err<1e-10
    report={'status':'passed-independent-saved-six-step-state-gates-and-SVD-reference','publicRiemannFacesChecked':checks,
        'publicBindings':PUBLIC_BINDINGS,'publicRiemannMaximumRelativeFluxDifference':maxerr,'KKTIndependentSolutionRelativeError':err,'patchComparison':stats,
        'physicalSolution':False,'newProductionSolvedAreaM2':0}
    report['case']=case
    target='independent_audit.json' if case=='patch02' else 'independent_audit_patch03.json'
    (OUT/target).write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
