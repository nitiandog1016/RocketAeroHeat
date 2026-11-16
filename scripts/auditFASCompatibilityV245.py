"""模块职责：独立复核两层状态、FAS代数与热壁几何；运行成功与物理适用分别判定。"""
from pathlib import Path
import csv,struct,sys
import numpy as np
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import validate
import runFASCompatibilityV245 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CLOSED=True


def rows(p):
    with Path(p).open(encoding='utf8') as f:return list(csv.DictReader(f))


def frame(p):
    b=Path(p).read_bytes()
    if len(b)<24:raise ValueError('HEADER')
    magic,n,k=struct.unpack('<8sQQ',b[:24])
    if magic!=b'RAH245LE' or k!=45 or n not in (4,24) or len(b)!=24+n*k*8:raise ValueError('FRAME_IDENTITY_OR_LENGTH')
    a=np.frombuffer(b[24:],dtype='<f8').reshape(n,k)
    if not np.array_equal(a[:,0],np.arange(n)) or not np.isfinite(a[:,1:9]).all() or np.any(a[:,1:3]<=0):raise ValueError('CURRENT_STATE')
    return a


def audit():
    """模块一：独立体积加权与物理能量线性变换，不复用C++层间实现。"""
    out=BASE/'audit01';out.mkdir(exist_ok=False);case=BASE/'observe01';errors={};checks=0
    def check(name,a,b,tol=2e-9):
        nonlocal checks
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[name]=max(errors.get(name,0),e)
        assert e<tol,(name,e);checks+=int(a.size)
    def physical(u):
        u=np.array(u,copy=True);u[:,4]-=300000*u[:,0];return u
    links=np.loadtxt(str(case/'children.csv'),delimiter=',',skiprows=1);co=links[:,0].astype(int);fi=links[:,1].astype(int)
    assert sorted(fi)==list(range(24));weight=links[:,2]/links[:,3]
    P=np.zeros((4,24));P[co,fi]=weight;check('volumePartition',P.sum(axis=1),np.ones(4))
    W=[np.loadtxt(str(case/'walls.csv'),delimiter=',',skiprows=1)]
    walls=[W[0][W[0][:,0]==l] for l in range(2)];wi=[w[:,1].astype(int) for w in walls]
    observed=rows(case/'frames.csv');off=rows(BASE/'off01/frames.csv')
    assert len(observed)==232 and len(off)==12
    lookup={(int(r['completed']),r['label'],int(r['level'])):r for r in observed if r['label'] in ('initial','end')}
    for r in off:
        other=lookup[(int(r['completed']),r['label'],int(r['level']))]
        assert sha(BASE/'off01'/r['file'])==sha(case/other['file'])
    for n in ('history.csv','restart_flow.csv','children.csv','walls.csv'):
        assert sha(case/n)==sha(BASE/'off01'/n),n
    for name in ('off01','observe01'):
        completion=read(BASE/name/'native_completion.json')
        assert completion['normalReturn'] and completion['updates']==5 and completion['coarsePoints']==4
    data={int(r['seq']):frame(case/r['file']) for r in observed}
    def get(step,label,level):
        matches=[data[int(r['seq'])] for r in observed if int(r['inner'])==step and r['label']==label and int(r['level'])==level]
        assert len(matches)==1,(step,label,level);return matches[0]
    for step in range(5):
        f=get(step,'restrict_in',0);c=get(step,'restrict_out',1);expected=P@f[:,2:7];expected[wi[1],1:4]=0
        check('restrictedFive',c[:,2:7],expected);check('restrictionPhysicalEnergy',physical(c[:,2:7])[:,[0,4]],(P@physical(f[:,2:7]))[:,[0,4]])
        f=get(step,'forcing_in',0);c=get(step,'forcing_in',1);a=get(step,'forcing_out',1)
        expected=np.zeros((4,5));np.add.at(expected,co,.5*f[fi,34:39]);expected[wi[1],1:4]=0;expected-=c[:,34:39]
        check('forcingFive',a[:,39:44],expected);check('forcingPhysicalEnergy',physical(a[:,39:44]),physical(expected))
        f=get(step,'correction_in',0);c=get(step,'correction_in',1);a=get(step,'correction_out',1);b=get(step,'correction_out',0)
        expected=c[:,2:7]-P@f[:,2:7];expected[wi[1],1:4]=0
        check('coarseCorrection',a[:,29:34],expected)
        injected=np.empty((24,5));injected[fi]=expected[co];check('injectedCorrection',b[:,34:39],injected)
        f=get(step,'prolong_in',0);a=get(step,'prolong_out',0)
        check('prolongedFive',a[:,2:7],f[:,2:7]+.5*f[:,34:39]);check('prolongPhysicalEnergy',physical(a[:,2:7]),physical(f[:,2:7])+.5*physical(f[:,34:39]))
        for l in range(2):
            a=get(step,'end',l);check('wallNoSlip',a[wi[l],3:6],np.zeros((len(wi[l]),3)))
        fine=get(step,'end',0);coarse=get(step,'end',1)
        check('SARestriction',coarse[:,7],P@fine[:,7]);mut=P@fine[:,8];mut[wi[1]]=0;check('eddyViscosityRestriction',coarse[:,8],mut)
        check('fineWallSA',fine[wi[0],7],np.zeros(len(wi[0])))
    # 模块二：Cv分段积分与Brent独立反解；验证所有保存守恒状态，不拿缓存温度作输入。
    gas=read(ROOT/'data/thermo/phase_aware_dry_air_v1.json');T=np.array(gas['temperatureK']);cp=np.array(gas['cpJkgK']);R=gas['gasConstantJkgK'];cv=cp-R
    knots=np.r_[0,np.cumsum(np.diff(T)*(cv[:-1]+cv[1:])/2)];ref=int(np.flatnonzero(T==298.15)[0]);knots-=knots[ref]+R*298.15
    def energy(temp):
        i=min(len(T)-2,int(np.searchsorted(T,temp,side='right')-1));dt=temp-T[i]
        return knots[i]+cv[i]*dt+.5*(cv[i+1]-cv[i])/(T[i+1]-T[i])*dt*dt
    cache={};temps={};statecount=0;air=PhaseAwareViscousAir()
    for r in observed:
        a=data[int(r['seq'])];u=a[:,3:6]/a[:,2,None];ep=a[:,6]/a[:,2]-.5*np.sum(u*u,axis=1)-300000
        ts=[]
        for v in ep:
            if float(v) not in cache:cache[float(v)]=brentq(lambda t:energy(t)-v,T[0],T[-1],xtol=1e-10)
            ts.append(cache[float(v)])
        ts=np.array(ts);validate(a[:,2]*R*ts,ts);temps[int(r['seq'])]=ts;statecount+=len(a)
        if r['label']=='coarse_preprocessed' and r['level']=='1':
            check('coarsePrimitiveT',a[:,9],ts);check('coarsePrimitiveP',a[:,13],a[:,2]*R*ts);check('coarseVelocity',a[:,10:13],u)
            cpv=np.interp(ts,T,cp);cvv=cpv-R;tr=air.interpolated_transport(a[:,13],ts)
            check('coarseCp',a[:,20],cpv);check('coarseSound',a[:,16],np.sqrt(cpv/cvv*R*ts))
            check('coarsePhysicalH',a[:,15]-300000,ep+R*ts+.5*np.sum(u*u,axis=1))
            j=np.minimum(len(T)-2,np.searchsorted(T,ts,side='right')-1)
            secondary=np.column_stack((R*ts,a[:,2]*R/cvv,np.zeros(len(a)),1/cvv,np.zeros(len(a)),
                np.diff(air.mu)[j]/np.diff(T)[j],np.zeros(len(a)),np.diff(air.kappa)[j]/np.diff(T)[j]))
            check('coarseSecondary',a[:,21:29],secondary);check('coarseTransport',a[:,[17,19]],np.column_stack((tr['dynamicViscosityPas'],tr['conductivityWmK'])))
    # 原生空间阶段的热壁增量：采用相同阶段记录，只复算Fourier通量，不改变运行缓存。
    heat=[]
    for r in observed:
        if r['label']!='flow_before_strong':continue
        k=int(r['seq']);nxt=observed[k+1];assert nxt['label']=='flow_after_strong' and nxt['level']==r['level']
        a=data[k];b=data[k+1];l=int(r['level']);w=walls[l];ids=wi[l];nb=w[:,2].astype(int);normal=w[:,3:6]
        area=np.linalg.norm(normal,axis=1);dist=abs(np.sum(normal*(w[:,9:12]-w[:,6:9]),axis=1))/area
        q=(a[ids,19]+a[ids,20]*a[ids,18]/.9)*(a[nb,9]-300)/dist*area
        check('nativeWallHeatIncrement',physical(b[:,34:39])[:,4]-physical(a[:,34:39])[:,4],np.bincount(ids,weights=q,minlength=len(a)))
        check('strongWallMomentum',b[ids,35:38],np.zeros((len(ids),3)))
        heat.append(dict(inner=int(r['inner']),level=l,heatIntoWallW=float(sum(q))))
    # 模块三：相反方向整壳壁被聚成一点，不将近零向量法向当有效热壁。
    fine_normals={int(row[1]):row[3:6] for row in walls[0]};cancellation=[]
    for row in walls[1]:
        c=int(row[1]);child=[int(i) for i in fi[co==c] if int(i) in fine_normals]
        vectors=np.array([fine_normals[i] for i in child]);norms=np.linalg.norm(vectors,axis=1)
        check('wallNormalAggregation',row[3:6],vectors.sum(axis=0))
        ratio=float(np.linalg.norm(row[3:6])/sum(norms));unresolved=ratio<=64*np.finfo(float).eps
        cancellation.append(dict(coarsePoint=c,childWallPoints=child,sumFineNormalMagnitudes=float(sum(norms)),
            coarseNormalMagnitude=float(np.linalg.norm(row[3:6])),ratio=ratio,numericallyUnresolved=bool(unresolved)))
    geometry_ok=not any(r['numericallyUnresolved'] for r in cancellation)
    # 真实文件读取边界：错误魔数、截断、尺寸和非有限当前状态均拒绝。
    original=(case/observed[0]['file']).read_bytes();bad=[b'BAD245LE'+original[8:],original[:-1],original[:16]+struct.pack('<Q',44)+original[24:]]
    changed=bytearray(original);changed[24+2*8:24+3*8]=struct.pack('<d',float('nan'));bad.append(bytes(changed))
    for j,b in enumerate(bad):
        p=out/('bad%d.bin'%j);p.write_bytes(b)
        try:frame(p)
        except ValueError:continue
        raise AssertionError('MALFORMED_FRAME_ACCEPTED')
    report=dict(passed=True,softwareAlgebraPassed=True,integratedThermalCompatibilityPassed=geometry_ok,
      productionAdmitted=False,decision='reject-this-small-shell-wall-agglomeration',finePoints=24,coarsePoints=4,
      updatesEach=5,observedFrames=len(observed),independentStateDomainChecks=statecount,uniqueEnergyInversions=len(cache),
      independentComparedScalars=checks,maxScaledErrors=errors,observationNeutralBitEqual=True,
      sameNativeHistoryAndRestartBitEqual=True,wallCancellation=cancellation,wallHeatStages=heat,
      malformedFramesRejected=len(bad),physicalQualification=False,newQualifiedAreaM2=0,
      unknownFieldsZeroFilled=False,initialRHSUsedForScience=False,originalRocketRun=False,
      requestedPostSmooth=[1,0],effectivePostSmooth=[0,0],nativeFinePostSmoothReset=True)
    save(out/'report.json',report);print({k:report[k] for k in ('passed','integratedThermalCompatibilityPassed','independentStateDomainChecks','independentComparedScalars')})


def main():
    if CLOSED:raise SystemExit('v245 FAS audit closed; follow NEXT_DECISION')
    g=run.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    locks += [p for name in ('observe01','off01') for p in (BASE/name).iterdir() if p.is_file()]
    return g.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,256,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v245 FAS audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
