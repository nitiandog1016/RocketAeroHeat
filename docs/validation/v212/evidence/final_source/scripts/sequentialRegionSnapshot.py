"""顺序区域的冻结接口交换文件；只搬运并审核数据，不计算CFD通量或授热资格。

当前是离线适配层。原生导出/导入及残差等价尚未实现，不能直接接StaticProfile。
冻结的是一次外迭代内的邻域数据；受体内部复数Jacobian仍由原生求解器计算。
"""
from pathlib import Path
import hashlib, json, re, zipfile
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir

SCHEMA='RAH.FrozenMappedSnapshot.v1'
MODEL='rah-phase-aware-dry-air-eilmer-v1'
SOURCE_COMMIT='f53f4609a0331d48efee69a4e4f3c3598378cc03'
DATA_SHA='7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'
GEOMETRY=('x','y','z','volume','areaxy','iLength','jLength','kLength','L_min','L_max')
STATE=('rho','p','p_e','a','T','u','mu','k','sigma','quality','massf','rho_s',
       'velx','vely','velz','turb0','turb1','mu_t','k_t','S','h0')
RECON=('velx','vely','velz','turb0','turb1','rho_s','rho','p','T','u')
RECON_PARTS=('dx','dy','dz','phi','min','max')
VISCOUS=('velx','vely','velz','massf','T','turb0','turb1')
SHAPES={'faceIDs':(), 'donorCellIDs':(), 'geometry':(10,), 'faceGeometry':(7,),
        'state':(21,), 'reconstruction':(10,6), 'viscous':(7,3),
        'wallDistance':(), 'nativeIntegratedFlux':(6,)}
HASH_KEYS=('motherMeshSHA256','discreteOperatorSHA256','caseSHA256','donorFrameSHA256',
           'nativeBuildSHA256','wallDistanceFieldSHA256')
META_KEYS=set(HASH_KEYS)|{'schema','model','dataSHA256','sourceCommit','donorRegion','receiverRegion',
    'outerIteration','donorIteration','sweepPosition','physicalTimeS','origin','fieldLayout',
    'frozenDonor','fluxConvention','thermoInterpolator'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def layout():
    return dict(geometry=list(GEOMETRY),state=list(STATE),reconstruction=list(RECON),
        reconstructionParts=list(RECON_PARTS),viscous=list(VISCOUS),
        faceGeometry=['x','y','z','outwardNx','outwardNy','outwardNz','areaM2'],
        flux=['mass','momentumX','momentumY','momentumZ','totalEnergy','SA'])

def _require(condition, reason):
    if not condition: raise ValueError('SNAPSHOT_'+reason)

def _metadata(meta, expected):
    """模块一：绑定真正的母网格/物性/离散算子/轮次，调用者提供预期身份。"""
    _require(set(meta)==META_KEYS and set(expected)==META_KEYS,'METADATA_FIELDS')
    _require(meta==expected,'STALE_OR_WRONG_IDENTITY')
    _require(meta['schema']==SCHEMA and meta['model']==MODEL and meta['dataSHA256']==DATA_SHA
             and meta['sourceCommit']==SOURCE_COMMIT,'MODEL')
    for k in HASH_KEYS:
        _require(isinstance(meta[k],str) and re.fullmatch('[0-9a-f]{64}',meta[k]) is not None,'HASH_'+k)
    _require(meta['fieldLayout']==layout(),'LAYOUT')
    _require(meta['thermoInterpolator'] in ('pt','rhou','rhop','rhoT'),'THERMO_INTERPOLATOR')
    for k in ('outerIteration','donorIteration','sweepPosition'):
        _require(type(meta[k]) is int and meta[k]>=0,'ITERATION')
    _require(type(meta['physicalTimeS']) in (int,float) and np.isfinite(meta['physicalTimeS'])
             and meta['physicalTimeS']>=0,'TIME')
    for k in ('donorRegion','receiverRegion'):
        _require(isinstance(meta[k],str) and re.fullmatch('[A-Za-z0-9_-]+',meta[k]) is not None,'REGION')
    _require(meta['donorRegion']!=meta['receiverRegion'],'SELF_DONOR')
    _require(meta['origin'] in ('synthetic-control','native-residual-snapshot'),'ORIGIN')
    _require(meta['frozenDonor'] is True and meta['fluxConvention']=='recipient-outward-integrated-SI','CONVENTION')

def _arrays(arrays):
    _require(set(arrays)==set(SHAPES),'ARRAY_FIELDS')
    ids=np.asarray(arrays['faceIDs']); _require(ids.ndim==1 and len(ids)>0,'EMPTY_IDS')
    n=len(ids)
    for k,tail in SHAPES.items():
        a=np.asarray(arrays[k]); _require(a.shape==(n,)+tail,'SHAPE_'+k)
        # 不默默裁复数虚部、把float索引取整，或把低精度存储当原始双精度。
        dtype=np.dtype('int64') if k in ('faceIDs','donorCellIDs') else np.dtype('float64')
        _require(a.dtype==dtype,'DTYPE_'+k)
        # 原生未使用的重构槽可能保留NaN；其phi=-1语义在validate中逐项审核。
        if k not in ('reconstruction','viscous'):
            _require(np.isfinite(a).all(),'NONFINITE_'+k)
        else:
            _require(not np.isinf(a).any(),'INFINITE_'+k)
    _require(len(np.unique(ids))==n and np.all(ids>=0) and np.all(arrays['donorCellIDs']>=0),'IDS')
    g=arrays['geometry']; face=arrays['faceGeometry']
    _require(np.all(g[:,3:]>0) and np.all(g[:,8]<=g[:,9]),'GEOMETRY')
    _require(np.all(face[:,6]>0) and np.all(abs(np.linalg.norm(face[:,3:6],axis=1)-1)<1e-12),'FACE_GEOMETRY')
    _require(np.all(arrays['wallDistance']>0),'WALL_DISTANCE')
    return n

def validate(meta, arrays, expected, mapping):
    """模块二：精确ID回排并对照母网格；逐点核验相容状态，不修改输入值。"""
    _metadata(meta,expected); n=_arrays(arrays)
    keys=('faceIDs','donorCellIDs','geometry','faceGeometry','wallDistance')
    _require(set(mapping)==set(keys),'MAP_FIELDS')
    for k in keys:
        m=np.asarray(mapping[k]); _require(m.shape==arrays[k].shape and m.dtype==arrays[k].dtype,'MAP_'+k)
        _require(np.isfinite(m).all(),'MAP_NONFINITE')
    _require(len(np.unique(mapping['faceIDs']))==n,'MAP_DUPLICATE')
    _require(np.array_equal(np.sort(arrays['faceIDs']),np.sort(mapping['faceIDs'])),'MAP_IDS')
    index={int(v):i for i,v in enumerate(arrays['faceIDs'])}
    order=np.array([index[int(v)] for v in mapping['faceIDs']],dtype=np.int64)
    for k in keys:
        _require(np.array_equal(arrays[k][order],mapping[k]),'MAP_MISMATCH_'+k)
    state=arrays['state']; col={name:state[:,i] for i,name in enumerate(STATE)}
    gas=PhaseAwareViscousAir(); physical=gas.evaluate(col['p'],col['T'])
    transport=gas.interpolated_transport(col['p'],col['T'])
    wanted=dict(rho=physical['densityKgm3'], a=physical['soundSpeedMps'],u=physical['internalEnergyJkg'],
        mu=transport['dynamicViscosityPas'],k=transport['conductivityWmK'],rho_s=physical['densityKgm3'])
    speed2=sum(col[k]**2 for k in ('velx','vely','velz'))
    wanted['h0']=physical['enthalpyJkg']+speed2/2
    for k,v in wanted.items():
        scale_floor=1.0 if k in ('u','h0','a') else 1e-12
        _require(np.all(abs(col[k]-v)<=5e-10*np.maximum(abs(v),scale_floor)),'THERMO_'+k)
    _require(np.all(col['massf']==1) and np.all(col['quality']==1) and np.all(col['sigma']==0),'SINGLE_NEUTRAL_AIR')
    _require(np.all(col['p_e']>=0) and np.all(col['turb0']>=0) and np.all(col['mu_t']>=0)
             and np.all(col['k_t']>=0),'SA_OR_TRANSPORT')
    # 模块使用量必须有限；原生未使用槽按phi=-1原样保存，不能误将NaN转成已解状态。
    active={'velx','vely','velz','turb0'}|set(dict(pt=('p','T'),rhou=('rho','u'),rhop=('rho','p'),rhoT=('rho','T'))[meta['thermoInterpolator']])
    for j,name in enumerate(RECON):
        r=arrays['reconstruction'][:,j,:]
        used=np.ones(n,dtype=bool) if name in active else r[:,3]!=-1
        _require(np.isfinite(r[used]).all(),'ACTIVE_RECONSTRUCTION_'+name)
        _require(np.all(r[used,3]>=0) and np.all(r[used,4]<=r[used,5]),'LIMITER_'+name)
    _require(np.isfinite(arrays['viscous'][:,:6,:]).all(),'ACTIVE_VISCOUS')
    _require(np.all(arrays['viscous'][:,VISCOUS.index('massf'),:]==0),'SINGLE_SPECIES_GRADIENT')
    # 返回独立副本，阻止调用者之后修改缓存的数据；不是原生Jacobian验收。
    frozen={k:np.array(v[order],copy=True) for k,v in arrays.items()}
    for a in frozen.values(): a.flags.writeable=False
    return frozen

def write_snapshot(path, meta, arrays, expected, mapping):
    """模块三：只写新的交换文件；上游运行快照及轮次由调用者事前绑定。"""
    checked=validate(meta,arrays,expected,mapping)
    header=np.frombuffer(json.dumps(meta,sort_keys=True,ensure_ascii=True,allow_nan=False).encode('utf-8'),dtype=np.uint8)
    path=Path(path)
    with path.open('xb') as stream:
        np.savez(stream,metadata=header,**checked)
    return dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size,
                origin=meta['origin'],nativeAdapterTested=False,runAuthorized=False,physicalQualification=False)

def read_snapshot(path, expected_sha, expected, mapping):
    """模块四：调用方给出已绑定文件摘要；禁用pickle，损坏/陈旧/错面必须拒绝。"""
    path=Path(path); _require(path.stat().st_size<=16*1024**2,'FILE_SIZE')
    _require(sha(path)==expected_sha,'FILE_SHA')
    with zipfile.ZipFile(str(path)) as z:
        names=z.namelist(); _require(len(names)==len(set(names)),'DUPLICATE_MEMBER')
        _require(set(names)=={k+'.npy' for k in SHAPES}|{'metadata.npy'},'MEMBERS')
        _require(sum(i.file_size for i in z.infolist())<=16*1024**2,'EXPANDED_SIZE')
    with np.load(str(path),allow_pickle=False) as z:
        raw=z['metadata']; _require(raw.dtype==np.uint8 and raw.ndim==1 and len(raw)<32768,'HEADER')
        meta=json.loads(raw.tobytes().decode('utf-8'))
        arrays={k:z[k] for k in SHAPES}
    return validate(meta,arrays,expected,mapping)

def characteristic_labels(frozen):
    """模块五：只显示局部无黏法向特征符号；不把它当完整黏性/回流边界条件。"""
    state=frozen['state']; n=frozen['faceGeometry'][:,3:6]
    vel=state[:,[STATE.index(k) for k in ('velx','vely','velz')]]
    a=state[:,STATE.index('a')]; un=np.sum(vel*n,axis=1)
    eigen=np.column_stack((un-a,un,un+a))
    labels=np.full(len(un),'mixed',dtype='<U20')
    labels[np.all(eigen>0,axis=1)]='supersonic-outflow'
    labels[np.all(eigen<0,axis=1)]='supersonic-inflow'
    return dict(normalMach=un/a,characteristicSpeeds=eigen,labels=labels,
        boundaryQualified=False)

if __name__=='__main__':
    raise SystemExit('仅供离线接口审核；没有原生适配、网格或CFD运行准入。')
