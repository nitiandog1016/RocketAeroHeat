"""模块职责：独立读取 v2 原生冻结文件，检查身份、有效槽和相容物性。

此接口不兼容地替代 v1，也不把制造状态认作流动解。未使用槽按位保留。
"""
from pathlib import Path
import hashlib,json,struct,re
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir

SCHEMA='RAH.FrozenMappedSnapshot.v2'
GEOMETRY=('x','y','z','volume','areaxy','iLength','jLength','kLength','L_min','L_max')
STATE=('rho','p','p_e','a','T','u','mu','k','sigma','quality','massf','rho_s',
       'velx','vely','velz','turb0','turb1','mu_t','k_t','S','Bx','By','Bz','psi','divB')
RECON=('velx','vely','velz','rho','p','T','u','Bx','By','Bz','psi','turb0','turb1','rho_s')
VISCOUS=('velx','vely','velz','T','massf','turb0','turb1')
META_KEYS={'schema','model','dataSHA256','sourceCommit','motherMeshSHA256','geometryMappingSHA256',
    'initialFrameSHA256','sourceRecipeSHA256','discreteConfigSHA256','adapterSHA256','nativeBuildSHA256',
    'donorRegion','receiverRegion','origin','outerIteration','donorIteration','sweepPosition','physicalTimeS',
    'configuration','layout'}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(ok,code):
    if not ok:raise ValueError('RAH_FILE_'+code)
def pairs(items):
    d={}
    for k,v in items:
        require(k not in d,'DUPLICATE_KEY');d[k]=v
    return d
def load(path):return json.loads(Path(path).read_text(encoding='utf8'),object_pairs_hook=pairs)
def scalar(value,imag=0.):return struct.pack('>dd',float(value),float(imag)).hex()
def unpack(values):
    require(isinstance(values,list) and len(values)==140,'SCALAR_COUNT')
    require(all(isinstance(v,str) and re.fullmatch('[0-9a-f]{32}',v) for v in values),'SCALAR_ENCODING')
    return np.array([struct.unpack('>dd',bytes.fromhex(v)) for v in values])
def validate(doc,expected,mapping):
    """模块一：母面用母节点规范键，原生 cell ID 单独绑定；不凭行号猜编号。"""
    require(set(doc)=={'metadata','records'},'TOP_FIELDS');meta=doc['metadata']
    require(set(meta)==META_KEYS and meta==expected,'IDENTITY')
    require(meta['schema']==SCHEMA,'SCHEMA')
    require(meta['model']=='rah-phase-aware-dry-air-eilmer-v1','MODEL')
    require(meta['configuration']=='3D-single-neutral-SA-Edwards-rhou-noMHD-noElectric','CONFIG')
    require(meta['layout']=='geometry10-state25-recon14x6-viscous7x3-complexHex64','LAYOUT')
    for k in ('outerIteration','donorIteration','sweepPosition'):
        require(type(meta[k]) is int and meta[k]>=0,'ITERATION')
    require(type(meta['physicalTimeS']) is float and np.isfinite(meta['physicalTimeS']) and meta['physicalTimeS']>=0,'TIME')
    rows=doc['records'];require(isinstance(rows,list) and len(rows)==len(mapping)>0,'COUNT')
    byid={}
    for row in rows:
        require(set(row)=={'receiverFaceID','motherCellID','motherFaceNodes','values'},'ROW_FIELDS')
        i=row['receiverFaceID'];require(type(i) is int and i>=0 and i not in byid,'DUPLICATE_OR_INVALID_FACE')
        byid[i]=row
    values=[]
    for i,b in mapping.items():
        require(i in byid,'FACE_ID');r=byid[i]
        require(r['motherCellID']==b['motherCellID'] and r['motherFaceNodes']==b['motherFaceNodes'],'MOTHER_MAPPING')
        require(r['values'][:10]==b['geometry'],'SOURCE_GEOMETRY');values.append(unpack(r['values']))
    a=np.array(values);require(not np.isinf(a).any(),'INFINITY')
    require(np.all((a[:,:,1]==0)|np.isnan(a[:,:,1])),'COMPLEX')
    def active(x):require(np.isfinite(x[:,:,0]).all() and np.all(x[:,:,1]==0),'ACTIVE_NONREAL_OR_NONFINITE')
    active(a[:,[0,1,2,3,5,6,7,8,9]])
    require(np.all(a[:,[3,5,6,7,8,9],0]>0) and np.all(a[:,8,0]<=a[:,9,0]),'GEOMETRY')
    s=a[:,10:35];active(s[:,[0,1,3,4,5,6,7,8,9,10,12,13,14,15,17,18,19]])
    require(np.all(s[:,10,0]==1) and np.all(s[:,8,0]==0) and np.all(s[:,9,0]==1),'NEUTRAL_AIR')
    require(np.all(s[:,[15,17,18],0]>=0),'SA_TRANSPORT')
    gas=PhaseAwareViscousAir();props=gas.evaluate(s[:,1,0],s[:,4,0]);trans=gas.interpolated_transport(s[:,1,0],s[:,4,0])
    for i,w in [(0,props['densityKgm3']),(3,props['soundSpeedMps']),(5,props['internalEnergyJkg']),
                (6,trans['dynamicViscosityPas']),(7,trans['conductivityWmK'])]:
        require(np.all(abs(s[:,i,0]-w)<=5e-10*np.maximum(1. if i in (3,5) else 1e-12,abs(w))),'THERMO')
    r=a[:,35:119].reshape((-1,14,6,2))
    for j in range(14):
        used=np.ones(len(a),bool) if j in (0,1,2,3,6,11) else np.zeros(len(a),bool)
        if j in (4,5,13):used|=r[:,j,3,0]!=-1
        active(r[used,j]);require(np.all(r[used,j,3,0]>=0) and np.all(r[used,j,4,0]<=r[used,j,5,0]),'LIMITER')
    active(a[:,119:137]);require(np.all(a[:,131:134,0]==0),'SINGLE_SPECIES_GRADIENT')
    a.flags.writeable=False
    return a
def read_snapshot(path,expected_sha,expected,mapping):
    """模块二：读取前检查大小及调用方的已绑定摘要；不自动信任文件内的来源声明。"""
    require(Path(path).stat().st_size<=16*1024**2,'SIZE');require(sha(path)==expected_sha,'SHA')
    return validate(load(path),expected,mapping)
