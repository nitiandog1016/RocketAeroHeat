"""接口交换的独立数据/误用控制；不把合成剖面当真实求解或守恒验证。"""
from pathlib import Path
import sys, json, copy, time, shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import sequentialRegionSnapshot as api
from phaseAwareViscousAir import PhaseAwareViscousAir

def fixture():
    """模块一：四种法向传播方向及合法负内能，显式标成合成数据。"""
    n=4;gas=PhaseAwareViscousAir();p=np.array([1000.,2000.,2500.,10000.]);t=np.array([220.,250.,300.,400.])
    q=gas.evaluate(p,t);tr=gas.interpolated_transport(p,t)
    vel=np.array([[-1000.,20.,10.],[1000.,15.,20.],[-50.,0.,1.],[50.,20.,5.]])
    vals=dict(rho=q['densityKgm3'],p=p,p_e=p,a=q['soundSpeedMps'],T=t,u=q['internalEnergyJkg'],
        mu=tr['dynamicViscosityPas'],k=tr['conductivityWmK'],sigma=np.zeros(n),quality=np.ones(n),
        massf=np.ones(n),rho_s=q['densityKgm3'],velx=vel[:,0],vely=vel[:,1],velz=vel[:,2],
        turb0=np.full(n,1e-5),turb1=np.ones(n),mu_t=np.full(n,2e-6),k_t=np.full(n,.001),S=np.zeros(n),
        h0=q['enthalpyJkg']+(vel*vel).sum(axis=1)/2)
    recon=np.ones((n,10,6));recon[:,:,:3]=np.arange(n*30).reshape(n,10,3)*1e-4
    recon[:,:,4]=-1e6;recon[:,:,5]=1e6
    # n_species=1、SA一方程及rhou时的未使用槽保留原生哨兵，不给NaN替补物理值。
    for name in ('turb1','rho_s','p','T'):
        recon[:,api.RECON.index(name),:]=np.nan;recon[:,api.RECON.index(name),3]=-1.
    visc=np.arange(n*21,dtype=float).reshape(n,7,3)*.002;visc[:,3,:]=0;visc[:,6,:]=np.nan
    geom=np.ones((n,10));geom[:,:3]=np.array([[.1,0,0],[.2,.1,0],[.3,0,.1],[.4,.1,.1]])
    geom[:,3]=1e-8;geom[:,4]=1e-5;geom[:,5:8]=.01;geom[:,8]=.001;geom[:,9]=.02
    face=np.zeros((n,7));face[:,:3]=geom[:,:3]+[.005,0,0];face[:,3]=1.;face[:,6]=1e-4
    arrays=dict(faceIDs=np.array([31,12,45,19],dtype=np.int64),donorCellIDs=np.array([5,8,5,9],dtype=np.int64),
        geometry=geom,faceGeometry=face,state=np.column_stack([vals[k] for k in api.STATE]),
        reconstruction=recon,viscous=visc,wallDistance=np.full(n,.0003),
        nativeIntegratedFlux=np.arange(n*6,dtype=float).reshape(n,6)*.1)
    # 以上通量只是数据搬运哨兵，不是本脚本自行算出的或原生实际通量。
    meta={k:format(i+1,'064x') for i,k in enumerate(api.HASH_KEYS)}
    meta.update(schema=api.SCHEMA,model=api.MODEL,dataSHA256=api.DATA_SHA,sourceCommit=api.SOURCE_COMMIT,
        donorRegion='nose',receiverRegion='contraction',outerIteration=2,donorIteration=2,sweepPosition=1,
        physicalTimeS=0.,origin='synthetic-control',fieldLayout=api.layout(),frozenDonor=True,
        fluxConvention='recipient-outward-integrated-SI',thermoInterpolator='rhou')
    mapping={k:arrays[k].copy() for k in ('faceIDs','donorCellIDs','geometry','faceGeometry','wallDistance')}
    return meta,arrays,mapping

def main():
    start=time.monotonic();out=ROOT/'outputs/v212_validation/packet01';out.mkdir(exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(),out/'executed_test.py')
    report=dict(passed=False,CFDRuns=0,nativeBuilds=0,nativeAdapterTested=False,
        physicalQualification=False,newQualifiedAreaM2=0)
    try:
        meta,arrays,mapping=fixture();expected=copy.deepcopy(meta)
        receipt=api.write_snapshot(out/'synthetic.npz',meta,arrays,expected,mapping)
        frozen=api.read_snapshot(out/'synthetic.npz',receipt['sha256'],expected,mapping)
        # 不以容差比较代替原值搬运；包含未使用NaN槽的原始位模式。
        assert all(frozen[k].tobytes()==v.tobytes() for k,v in arrays.items())
        perm=np.array([2,0,3,1]);shuffled={k:v[perm] for k,v in arrays.items()}
        reordered=api.validate(meta,shuffled,expected,mapping)
        assert all(reordered[k].tobytes()==v.tobytes() for k,v in arrays.items())
        assert not any(v.flags.writeable for v in frozen.values())
        assert np.any(frozen['state'][:,api.STATE.index('u')]<0)
        arrays['state'][0,0]+=1
        assert frozen['state'].tobytes()!=arrays['state'].tobytes()
        meta,arrays,mapping=fixture()
        labels=api.characteristic_labels(frozen)
        assert labels['labels'].tolist()==['supersonic-inflow','supersonic-outflow','mixed','mixed']
        assert not labels['boundaryQualified']
        rejected=[]
        def reject(name,fn):
            try: fn()
            except (ValueError,KeyError): rejected.append(name);return
            raise AssertionError('未拒绝: '+name)
        def changed(k,index,value):
            z=copy.deepcopy(arrays);z[k][index]=value;return z
        def check(z):return api.validate(meta,z,expected,mapping)
        # 模块二：数据身份/时间顺序/真实编号/状态/梯度必须分别有拒绝证据。
        for key in ('outerIteration','donorIteration','sweepPosition','model','donorFrameSHA256','motherMeshSHA256',
                    'caseSHA256','dataSHA256','discreteOperatorSHA256','nativeBuildSHA256','wallDistanceFieldSHA256',
                    'sourceCommit','receiverRegion','donorRegion','origin','thermoInterpolator','frozenDonor','fluxConvention'):
            m=copy.deepcopy(meta);m[key]='wrong'
            reject('wrong-'+key,lambda m=m:api.validate(m,arrays,expected,mapping))
        for k in ('reconstruction','viscous','nativeIntegratedFlux','wallDistance'):
            z=copy.deepcopy(arrays);z.pop(k);reject('missing-'+k,lambda z=z:check(z))
        for name,z in [
            ('duplicate-face',changed('faceIDs',0,12)),('wrong-face',changed('faceIDs',0,30)),
            ('wrong-donor',changed('donorCellIDs',0,6)),('mirror-position',changed('geometry',(0,0),.09)),
            ('reversed-normal',changed('faceGeometry',(0,3),-1)),('wrong-volume',changed('geometry',(0,3),2e-8)),
            ('wrong-wall-distance',changed('wallDistance',0,.1)),
            ('zero-area',changed('faceGeometry',(0,6),0)),
            ('negative-pressure',changed('state',(0,api.STATE.index('p')),-1)),
            ('low-temperature',changed('state',(0,api.STATE.index('T')),59)),
            ('wrong-density',changed('state',(0,api.STATE.index('rho')),1)),
            ('wrong-total-enthalpy',changed('state',(0,api.STATE.index('h0')),1)),
            ('wrong-viscosity',changed('state',(0,api.STATE.index('mu')),1e-3)),
            ('negative-SA',changed('state',(0,api.STATE.index('turb0')),-1)),
            ('active-gradient-nan',changed('reconstruction',(0,0,0),np.nan)),
            ('active-gradient-inf',changed('reconstruction',(0,0,0),np.inf)),
            ('active-viscous-nan',changed('viscous',(0,0,0),np.nan)),
            ('disabled-active-rhou',changed('reconstruction',(0,api.RECON.index('rho'),3),-1)),
            ('activated-unused-nan',changed('reconstruction',(0,api.RECON.index('p'),3),1)),
            ('wrong-mass-gradient',changed('viscous',(0,3,0),1)),
            ('invalid-extrema',changed('reconstruction',(0,0,4),2e6))]:
            reject(name,lambda z=z:check(z))
        z=copy.deepcopy(arrays);z['state']=z['state'].astype(complex);reject('complex-donor',lambda:check(z))
        z=copy.deepcopy(arrays);z['faceIDs']=z['faceIDs'].astype(float);reject('floating-IDs',lambda:check(z))
        z=copy.deepcopy(arrays);z['state']=z['state'].astype(np.float32);reject('low-precision-state',lambda:check(z))
        m=copy.deepcopy(meta);m['origin']='native-residual-snapshot'
        reject('synthetic-as-native',lambda:api.read_snapshot(out/'synthetic.npz',receipt['sha256'],m,mapping))
        reject('file-corruption',lambda:api.read_snapshot(out/'synthetic.npz','0'*64,expected,mapping))
        # 保护已交付/并发文件；本测试不改UI或旧流动场。
        protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
        assert len(protected)==29 and all(api.sha(ROOT/k)==v for k,v in protected.items())
        report.update(passed=True,roundtripBitwise=True,permutationBitwise=True,
            frozenIndependentCopy=True,unusedNativeSentinelsPreserved=True,negativeEnergyPreserved=True,
            characteristicLabels=labels['labels'].tolist(),rejected=rejected,rejectedCount=len(rejected),
            protectedFilesUnchanged=29,snapshotReceipt=receipt,
            inputSHA256={p.relative_to(ROOT).as_posix():api.sha(p) for p in [Path(__file__).resolve(),
                ROOT/'scripts/sequentialRegionSnapshot.py',ROOT/'docs/validation/v212/PLAN.md',
                ROOT/'data/thermo/phase_aware_dry_air_v1.json']})
    except Exception as exc:
        report['error']=repr(exc);raise
    finally:
        report['seconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf8')
    assert report['seconds']<60
    print(json.dumps({k:report[k] for k in ('passed','rejectedCount','roundtripBitwise','protectedFilesUnchanged','seconds')}))

if __name__=='__main__':main()
