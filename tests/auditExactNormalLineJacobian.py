# -*- coding: utf-8 -*-
"""用独立标量Brent/公开压力曲线复算原面Euler导数，不运行另一套NS求解器。"""
from pathlib import Path
import ast,hashlib,json,math,os,sys
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq,bisect

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';REF=ROOT/'outputs/v82_validation'
case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_jacobian.json');assert not target.exists()
file=VAL/case/'exact_branch_jacobian_evidence.mat';source=REF/'clawpack_exact_euler.py';probe=REF/'probe_exact_riemann.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in [file,source,probe]}
assert sha(source)=='2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'
assert sha(probe)=='402ffeee97aed006c6784fa6b46d3276a3ab428897d77d07921308d9708c067a'
allowed={'pospart','sound_speed','beta','integral_curve_1','integral_curve_3','hugoniot_locus_1','hugoniot_locus_3'}
nodes=[n for n in ast.parse(source.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name in allowed]
assert {n.name for n in nodes}==allowed
reference={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),reference)
nodes=[n for n in ast.parse(probe.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name=='exact']
assert len(nodes)==1
namespace=dict(np=np,math=math,brentq=brentq,bisect=bisect,reference=reference)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(probe),'exec'),namespace)
s=loadmat(file,simplify_cells=True);gas=s['gas'];gamma=gas['gamma'];R=gas['R'];cv=R/(gamma-1)
regions=['left-state','left-star','left-rarefaction-fan','right-state','right-star','right-rarefaction-fan']
def flux(ql,qr,n):
    primitives=[]
    for q in [ql,qr]:
        rho=q[0];u=q[1:4]/rho;p=(gamma-1)*(q[4]-.5*rho*np.dot(u,u));primitives.append((rho,u,p))
    rl,vl,pl=primitives[0];rr,vr,pr=primitives[1];ul=float(vl@n);ur=float(vr@n)
    result=namespace['exact']((rl,ul,pl),(rr,ur,pr),gamma);rho,un,p=result['sample']
    tangent=vl-ul*n if result['side']=='left' else vr-ur*n;u=tangent+un*n;energy=p/(gamma-1)+.5*rho*np.dot(u,u)
    return np.r_[rho*un,rho*u*un+p*n,(energy+p)*un],regions.index(result['sampleRegion'])+1
maximum=0.;base_error=0.;tested=0;excluded=0;case_errors=[]
for j,(ql,qr,n) in enumerate(zip(s['QL'],s['QR'],s['normal'])):
    base,region=flux(ql,qr,n);assert region==s['audit']['sampleRegion'][j]
    base_error=max(base_error,float(np.max(abs(base-s['audit']['baseFlux'][j]))/max(1,np.max(abs(base)))))
    worst=0.
    for side in range(2):
        q=[ql,qr][side];rho=q[0];u=q[1:4]/rho;T=(q[4]-.5*rho*np.dot(u,u))/(rho*cv)
        speed=max(np.linalg.norm(u),math.sqrt(gamma*R*T));basis=np.r_[rho,np.repeat(rho*speed,3),rho*speed*speed]
        for b in range(5):
            h=1e-6*basis[b];plus=q.copy();minus=q.copy();plus[b]+=h;minus[b]-=h
            fp,rp=flux(plus,qr,n) if side==0 else flux(ql,plus,n)
            fm,rm=flux(minus,qr,n) if side==0 else flux(ql,minus,n)
            if rp!=region or rm!=region:excluded+=1;continue
            actual=(fp-fm)/(2*h);expected=s['left' if side==0 else 'right'][:,b,j]
            error=float(np.max(abs(actual-expected))/max(1,np.max(abs(expected))));worst=max(worst,error);tested+=1
    case_errors.append(worst);maximum=max(maximum,worst)
assert base_error<3e-12 and maximum<5e-7 and tested==150 and excluded==0
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-Brent-public-pressure-curve-selected-branch-Jacobian-FD',case=case,processID=os.getpid(),
         baseFluxRelativeDifference=base_error,maximumIndependentFDRelativeDifference=maximum,caseRelativeDifferences=case_errors,
         finiteDifferenceColumns=tested,crossBranchColumnsExcluded=excluded,globalSmoothJacobianClaimed=False,
         fullIndependentNSSolver=False,physicalDiscretizationReplaced=False,newProductionSolvedAreaM2=0,
         inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
