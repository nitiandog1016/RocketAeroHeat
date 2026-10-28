# -*- coding: utf-8 -*-
"""对真实收缩段受限正状态作精确理想气体Riemann局部比较，不导入外部CFD或替换原残差。"""
from pathlib import Path
import ast,hashlib,json,math,os
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq,bisect
ROOT=Path(__file__).resolve().parents[2];VAL=Path(__file__).resolve().parent
target=VAL/'actual_exact_riemann_comparison.json';assert not target.exists()
source=VAL/'clawpack_exact_euler.py';assert hashlib.sha256(source.read_bytes()).hexdigest()=='2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'
# 只提取已逐行检查的公开纯数学压力曲线；不导入上游绘图/交互/文件功能。
allowed={'pospart','sound_speed','beta','integral_curve_1','integral_curve_3','hugoniot_locus_1','hugoniot_locus_3'}
tree=ast.parse(source.read_text(encoding='utf8'));nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in allowed]
assert {n.name for n in nodes}==allowed
for definition in nodes:
    for node in ast.walk(definition):
        assert not isinstance(node,(ast.Import,ast.ImportFrom,ast.Global,ast.Nonlocal,ast.With,ast.Try))
        if isinstance(node,ast.Call):
            assert (isinstance(node.func,ast.Name) and node.func.id in allowed) or (isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id=='np' and node.func.attr in ['maximum','sqrt'])
reference={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),reference)

def exact(left,right,gamma):
    rl,ul,pl=left;rr,ur,pr=right;al=math.sqrt(gamma*pl/rl);ar=math.sqrt(gamma*pr/rr);gm=gamma-1
    margin=2*(al+ar)/gm-(ur-ul)
    if margin<=0:raise ValueError('vacuum Riemann problem: not assigned a fabricated positive pressure')
    exponent=gm/(2*gamma)
    def wave(p,rho,pk,a):
        if p>pk:return (p-pk)*math.sqrt(2/((gamma+1)*rho)/(p+gm/(gamma+1)*pk))
        return 2*a/gm*((p/pk)**exponent-1)
    equation=lambda p:wave(p,rl,pl,al)+wave(p,rr,pr,ar)+ur-ul
    hi=max(pl,pr)
    while equation(hi)<0:hi*=2
    ps=brentq(equation,0,hi,xtol=1e-14,rtol=8*np.finfo(float).eps)
    us=.5*(ul+ur+wave(ps,rr,pr,ar)-wave(ps,rl,pl,al))
    # 第三方原压力曲线采用独立二分交点。参考的pospart在此正压括区没有触发。
    phiL=lambda p:reference['hugoniot_locus_1' if p>=pl else 'integral_curve_1'](p,rl,ul,pl,gamma)
    phiR=lambda p:reference['hugoniot_locus_3' if p>=pr else 'integral_curve_3'](p,rr,ur,pr,gamma)
    lower=min(pl,pr)*1e-12;assert lower>1e-15
    ref=bisect(lambda p:float(phiL(p)-phiR(p)),lower,hi,xtol=1e-13,rtol=8*np.finfo(float).eps,maxiter=200)
    assert abs(ps-ref)/max(1,abs(ps))<2e-12 and abs(us-float(phiL(ref)))/max(1,abs(us))<2e-12
    def density(rho,pk):
        if ps>pk:
            ratio=ps/pk;return rho*(ratio+gm/(gamma+1))/(gm/(gamma+1)*ratio+1)
        return rho*(ps/pk)**(1/gamma)
    if us>=0:
        side='left';rs=density(rl,pl)
        if ps>pl:
            speed=ul-al*math.sqrt((gamma+1)/(2*gamma)*ps/pl+gm/(2*gamma));sample=left if speed>=0 else (rs,us,ps);region='left-state' if speed>=0 else 'left-star'
        else:
            head=ul-al;tail=us-al*(ps/pl)**exponent
            if head>=0:sample=left;region='left-state'
            elif tail<=0:sample=(rs,us,ps);region='left-star'
            else:
                sound=2/(gamma+1)*(al+.5*gm*ul);sample=(rl*(sound/al)**(2/gm),sound,pl*(sound/al)**(2*gamma/gm));region='left-rarefaction-fan'
    else:
        side='right';rs=density(rr,pr)
        if ps>pr:
            speed=ur+ar*math.sqrt((gamma+1)/(2*gamma)*ps/pr+gm/(2*gamma));sample=right if speed<=0 else (rs,us,ps);region='right-state' if speed<=0 else 'right-star'
        else:
            head=ur+ar;tail=us+ar*(ps/pr)**exponent
            if head<=0:sample=right;region='right-state'
            elif tail>=0:sample=(rs,us,ps);region='right-star'
            else:
                sound=2/(gamma+1)*(ar-.5*gm*ur);sample=(rr*(sound/ar)**(2/gm),-sound,pr*(sound/ar)**(2*gamma/gm));region='right-rarefaction-fan'
    assert all(math.isfinite(x) for x in sample) and sample[0]>0 and sample[2]>0
    return dict(pressureStar=ps,velocityStar=us,sample=list(sample),sampleRegion=region,side=side,vacuumVelocityMargin=margin,
                originalClawpackCurvePressureReference=ref,pressureRootResidual=equation(ps))

gamma=1.4;sod=exact((1,0,1),(.125,0,.1),gamma);assert abs(sod['pressureStar']-.30313017805064685)<2e-13 and abs(sod['velocityStar']-.9274526200489499)<2e-13
vacuumRejected=False
try:exact((1,-5,.4),(1,5,.4),gamma)
except ValueError:vacuumRejected=True
assert vacuumRejected
domainPath=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat';statePath=ROOT/'outputs/v81_validation/arithmetic02/steady_krylov_evidence.mat'
d=loadmat(domainPath,simplify_cells=True)['domain'];s=loadmat(statePath,simplify_cells=True);f=d['openFaces'];face=70995-1;m=int(f['cellMinus'][face])-1;p=int(f['cellPlus'][face])-1;n=f['normal'][face];Q=s['result']['Q'];P=s['result']['info']['primitive'];R=s['gas']['R']
def evaluate_normal(normal,ql,qr):
    def prim(q):
        rho=q[0];u=q[1:4]/rho;pressure=(gamma-1)*(q[4]-.5*rho*np.dot(u,u));return rho,u,pressure
    rl,vl,pl=prim(ql);rr,vr,pr=prim(qr);ul=float(vl@normal);ur=float(vr@normal);answer=exact((rl,ul,pl),(rr,ur,pr),gamma)
    rho,un,pressure=answer['sample'];tangent=(vl-ul*normal) if answer['side']=='left' else (vr-ur*normal);u=tangent+un*normal
    energy=pressure/(gamma-1)+.5*rho*np.dot(u,u);flux=np.r_[rho*un,rho*u*un+pressure*normal,(energy+pressure)*un]
    answer['sampleTemperatureK']=pressure/(rho*R);answer['flux']=flux.tolist();return answer
actual=evaluate_normal(n,Q[m],Q[p]);assert 80<=actual['sampleTemperatureK']<=2000
reverse=evaluate_normal(-n,Q[p],Q[m]);assert np.max(abs(np.asarray(reverse['flux'])+actual['flux']))/max(1,np.max(np.abs(actual['flux'])))<3e-12
angle=.47;rotation=np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1.]])
def rotate(q):
    r=q.copy();r[1:4]=rotation@r[1:4];return r
rot=evaluate_normal(rotation@n,rotate(Q[m]),rotate(Q[p]));expected=np.array(actual['flux']);expected[1:4]=rotation@expected[1:4]
assert np.max(abs(np.array(rot['flux'])-expected))/max(1,np.max(abs(expected)))<3e-12
margin=json.loads((ROOT/'outputs/v81_validation/hllc_margin_diagnostic.json').read_text(encoding='utf8'))['cases'][-1]
report=dict(status='local-actual-contraction-exact-Riemann-math-comparison-not-integrated-flow-solver',processID=os.getpid(),
    sourceReferences=[dict(url='https://www.clawpack.org/riemann_book/html/Euler.html',scope='ideal-gas-Euler-pressure-curves-and-Riemann-waves'),
                      dict(path=source.name,sha256=hashlib.sha256(source.read_bytes()).hexdigest(),onlyPureMathCurvesExecuted=True,referenceFloorActiveAtRoot=False)],
    classicSodNondimensionalAnalyticCheck=sod,trueVacuumExplicitlyRejected=vacuumRejected,
    actualFaceID=face+1,ownerCV=m+1,neighbourCV=p+1,centroidM=f['centroid'][face].tolist(),leftPrimitiveRhoUvwT=P[m].tolist(),rightPrimitiveRhoUvwT=P[p].tolist(),
    actualExactRiemann=actual,davisHllcApproximateStarPressurePa=margin['minimumStarPressurePa'],
    exactStarPressureIsWallPressure=False,fullThreeDimensionalNSOrHeatSolved=False,originalHllcGuardChanged=False,
    reverseNormalAndRotationChecked=True,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [domainPath,statePath]})
target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');print(json.dumps(report))
