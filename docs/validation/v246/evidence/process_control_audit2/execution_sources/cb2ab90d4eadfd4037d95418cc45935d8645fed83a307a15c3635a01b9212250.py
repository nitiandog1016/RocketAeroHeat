"""模块职责：独立审核真实两层几何、全部原壁与相态初场；几何通过不授CFD/热资格。"""
from pathlib import Path
import csv,struct,sys
import numpy as np
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import validate
import runFASRocketPreflightV246 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def frame(path):
    """读取只接受固定格式，保留原生局部顺序与全球号。"""
    with Path(path).open('rb') as f:
        h=f.read(24)
        if len(h)!=24:raise ValueError('FRAME_HEADER')
        magic,n,k=struct.unpack('<8sQQ',h)
        if magic!=b'RAH246LE' or k!=33 or not 500<=n<=335726:raise ValueError('FRAME_IDENTITY')
        a=np.fromfile(f,dtype='<f8')
    if a.size!=n*k:raise ValueError('FRAME_SIZE')
    a=a.reshape(n,k)
    if not np.isfinite(a).all() or np.any(a[:,4:6]<=0):raise ValueError('FRAME_STATE')
    return a


def mesh(path):
    """流式独立原文本读取，避免复用原生重编号或粗层数组。"""
    with Path(path).open(encoding='utf8') as f:
        assert next(f).strip()=='NDIME= 3';cells=int(next(f).split('=')[1]);assert cells==705490
        for _ in range(cells):next(f)
        points=int(next(f).split('=')[1]);assert points==335726
        xyz=np.empty((points,3))
        for j in range(points):
            row=next(f).split();assert int(row[3])==j;xyz[j]=list(map(float,row[:3]))
        markers={};n=int(next(f).split('=')[1])
        for _ in range(n):
            name=next(f).split('=')[1].strip();count=int(next(f).split('=')[1]);faces=[]
            for j in range(count):
                row=list(map(int,next(f).split()));assert row[0] in (5,9);faces.append(row[1:])
            markers[name]=np.asarray(faces,dtype=np.int64)
        assert not f.read().strip()
    return xyz,markers


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);case=run.CASE;errors={};checks=0
    def check(name,a,b,tol=2e-9):
        nonlocal checks
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[name]=max(errors.get(name,0),e);assert e<tol,(name,e);checks+=a.size
    complete=read(case/'native_completion.json');assert complete['normalReturn'] and complete['flowUpdates']==complete['SAUpdates']==complete['rejections']==0
    a=frame(case/'fine_constructed.bin');b=frame(case/'coarse_constructed.bin');f=frame(case/'fine_initialized.bin');c=frame(case/'coarse_restricted.bin')
    N=len(f);M=len(c);assert (N,M)==(complete['finePoints'],complete['coarsePoints']) and M>=500 and M<N
    assert np.array_equal(a[:,0:5],f[:,0:5]) and np.array_equal(b[:,0:5],c[:,0:5])
    ids=f[:,0].astype(np.int64);assert np.array_equal(np.sort(ids),np.arange(N));xyz,markers=mesh(case/'mesh.su2')
    assert np.array_equal(f[:,1:4],xyz[ids]);assert sha(case/'mesh.su2')==sha(run.MESH)
    fw=np.loadtxt(str(case/'fine_walls.csv'),delimiter=',',skiprows=1,dtype=np.int64);assert len(fw)==19752 and np.all(fw[:,2]==3) and np.all(fw[:,-1]==-1)
    def keys(v):return sorted(map(tuple,np.sort(v,axis=1)))
    assert keys(fw[:,3:6])==keys(markers['wall'])
    tri=xyz[markers['wall']];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    check('originalWallAreaM2',areas.sum(),1.9646292565715553,2e-12)
    with (case/'children.bin').open('rb') as stream:
        magic,n,k=struct.unpack('<8sQQ',stream.read(24));assert magic==b'RAH246ID' and n==N and k==2
        pairs=np.fromfile(stream,dtype='<u8').reshape(N,2).astype(np.int64)
    co,fi=pairs.T;assert np.array_equal(np.sort(fi),np.arange(N)) and co.min()==0 and co.max()==M-1
    parents=np.empty(N,dtype=np.int64);parents[fi]=co;count=np.bincount(co,minlength=M);assert count.min()>0
    volumes=np.bincount(co,weights=f[fi,4],minlength=M);check('coarseVolume',c[:,4],volumes)
    expected=np.zeros((M,5));weights=f[fi,4]/c[co,4];np.add.at(expected,co,f[fi,5:10]*weights[:,None]);check('nativePublicRestriction',c[:,5:10],expected)
    physical=f[:,9]-300000*f[:,5];expected_energy=np.bincount(co,weights=physical[fi]*weights,minlength=M)
    check('restrictionPhysicalEnergy',c[:,9]-300000*c[:,5],expected_energy)
    assert np.array_equal(a[:,5:12],f[:,5:12]);assert np.array_equal(b[:,10:12],c[:,10:12])
    # 模块二：每种真实边界独立求子法向和，检查方向与邻点，单列近壁责任。
    boundary=np.loadtxt(str(case/'boundary.csv'),delimiter=',',skiprows=1);summary=[];reject=[];wall_parent_stats=[]
    eps=np.finfo(float).eps
    for marker in np.unique(boundary[:,1]).astype(int):
        fine=boundary[(boundary[:,0]==0)&(boundary[:,1]==marker)];coarse=boundary[(boundary[:,0]==1)&(boundary[:,1]==marker)]
        find=fine[:,3].astype(int);cind=coarse[:,3].astype(int);p=parents[find]
        assert len(set(find))==len(find) and len(set(cind))==len(cind) and set(p)==set(cind)
        assert len(np.unique(fine[:,2]))==1 and set(fine[:,2])==set(coarse[:,2]) and set(fine[:,15])==set(coarse[:,15])
        summed=np.zeros((M,3));np.add.at(summed,p,fine[:,6:9]);magnitudes=np.bincount(p,weights=np.linalg.norm(fine[:,6:9],axis=1),minlength=M)
        check('boundaryVectorAggregation',coarse[:,6:9],summed[cind])
        area=np.linalg.norm(coarse[:,6:9],axis=1);ratio=area/magnitudes[cind];bad=area<=64*eps*magnitudes[cind]
        for level,table,state in ((0,fine,f),(1,coarse,c)):
            point=table[:,3].astype(int);nb=table[:,5].astype(int);assert np.all((nb>=0)&(nb<len(state)))
            assert np.array_equal(table[:,9:12],state[point,1:4]) and np.array_equal(table[:,12:15],state[nb,1:4])
            normal=table[:,6:9];A=np.linalg.norm(normal,axis=1);delta=table[:,12:15]-table[:,9:12];length=np.linalg.norm(delta,axis=1)
            distance=abs(np.sum(normal*delta,axis=1))/A;bad_distance=(~np.isfinite(distance))|(distance<=64*eps*length)|(length==0)
            # 只有热壁使用法向邻点；外场邻点不是本次热壁拒绝量。
            if fine[0,15] and np.any(bad_distance):reject.append(dict(marker=marker,level=level,reason='unresolved-wall-normal-neighbor',points=point[bad_distance].tolist()))
            if fine[0,15] and level==1:
                for j,i in enumerate(cind):wall_parent_stats.append([int(i),int(count[i]),float(ratio[j]),float(distance[j]),float(A[j])])
        if np.any(bad):reject.append(dict(marker=marker,level=1,reason='unresolved-boundary-normal',points=cind[bad].tolist()))
        summary.append(dict(marker=int(marker),kind=int(fine[0,2]),viscous=bool(fine[0,15]),finePoints=len(fine),coarsePoints=len(coarse),
            minimumNormalRatio=float(min(ratio)),normalRatioQuantiles=np.quantile(ratio,[0,.01,.5,.99,1]).tolist(),unresolvedNormals=int(sum(bad))))
    np.savetxt(str(out/'coarse_wall_geometry.csv'),np.asarray(wall_parent_stats),delimiter=',',header='coarsePoint,children,normalRatio,normalDistanceM,normalAreaM2',comments='',fmt='%.17g')
    # 模块三：独立Cv积分/Brent，重复能量只反解一次；每个保存状态都检查联合域。
    gas=read(ROOT/'data/thermo/phase_aware_dry_air_v1.json');T=np.array(gas['temperatureK']);cp=np.array(gas['cpJkgK']);R=gas['gasConstantJkgK'];cv=cp-R
    knots=np.r_[0,np.cumsum(np.diff(T)*(cv[:-1]+cv[1:])/2)];ref=int(np.flatnonzero(T==298.15)[0]);knots-=knots[ref]+R*298.15
    def energy(temp):
        i=min(len(T)-2,int(np.searchsorted(T,temp,side='right')-1));dt=temp-T[i]
        return knots[i]+cv[i]*dt+.5*(cv[i+1]-cv[i])/(T[i+1]-T[i])*dt*dt
    cache={};ranges=[];air=PhaseAwareViscousAir()
    for name,state,refresh in (('fine_constructed',a,False),('coarse_constructed',b,False),('fine_initialized',f,True),('coarse_restricted',c,True)):
        vel=state[:,6:9]/state[:,5,None];ep=state[:,9]/state[:,5]-.5*np.sum(vel*vel,axis=1)-300000
        values,inverse=np.unique(ep,return_inverse=True)
        for value in values:
            if float(value) not in cache:cache[float(value)]=brentq(lambda t:energy(t)-value,T[0],T[-1],xtol=1e-10)
        ts=np.array([cache[float(value)] for value in values])[inverse];p=state[:,5]*R*ts;validate(p,ts)
        ranges.append(dict(name=name,points=len(state),minT=float(min(ts)),maxT=float(max(ts)),minP=float(min(p)),maxP=float(max(p))))
        if refresh:
            check('primitiveTemperature',state[:,12],ts);check('primitivePressure',state[:,16],p);check('primitiveVelocity',state[:,13:16],vel)
            cpv=np.interp(ts,T,cp);cvv=cpv-R;check('Cp',state[:,23],cpv);check('physicalH',state[:,18]-300000,ep+R*ts+.5*np.sum(vel*vel,axis=1))
            j=np.minimum(len(T)-2,np.searchsorted(T,ts,side='right')-1)
            secondary=np.column_stack((R*ts,state[:,5]*R/cvv,np.zeros(len(state)),1/cvv,np.zeros(len(state)),
              np.diff(air.mu)[j]/np.diff(T)[j],np.zeros(len(state)),np.diff(air.kappa)[j]/np.diff(T)[j]))
            check('secondary',state[:,24:32],secondary)
    report=dict(passed=True,geometryPreflightPassed=not reject,completeFASCompatibilityPassed=False,productionAdmitted=False,
        finePoints=N,coarsePoints=M,fineCells=705490,originalWallFaces=19752,originalWallAreaM2=float(sum(areas)),
        originalMeshBitEqual=True,originalCoordinateBitEqual=True,allOriginalWallFacesRecovered=True,allFineChildrenExactlyOnce=True,
        childCountQuantiles=np.quantile(count,[0,.01,.5,.99,1]).tolist(),boundary=summary,rejections=reject,
        independentComparedScalars=int(checks),independentStateDomainChecks=2*(N+M),uniqueEnergyInversions=len(cache),stateRanges=ranges,maxScaledErrors=errors,
        nativePublicVolumeRestriction=True,privateFASWallProjectionExecuted=False,fullSARestrictionExecuted=False,flowUpdates=0,SAUpdates=0,
        givenTwK=300,wallFluidStateIsInitial=True,solidTwSolved=False,heatFluxSolved=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print({k:report[k] for k in ('geometryPreflightPassed','finePoints','coarsePoints','independentStateDomainChecks','rejections')})


def main():
    if CLOSED:raise SystemExit('v246 geometry audit closed; follow NEXT_DECISION')
    g=run.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    locks += [p for p in run.CASE.iterdir() if p.is_file()]
    return g.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,512,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v246 geometry audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
