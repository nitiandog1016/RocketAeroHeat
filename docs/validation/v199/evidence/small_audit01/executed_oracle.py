"""只读离散守恒校核器：重放原生一阶Kurganov面通量，不参与任何CFD推进。"""
from pathlib import Path
import re,json
import numpy as np
from openFoamFields import internal_field,read_state,audit_state
from openFoamAirThermo import evaluate,R


def array_entry(text,key,count,vector=False):
    """模块一：读取准确边界value；不把freestreamValue或缺失场当作实际值。"""
    columns=3 if vector else 1
    m=re.search(r'\b'+key+r'\s+uniform\s+([^;]+);',text)
    if m:
        v=np.fromstring(m.group(1).replace('(',' ').replace(')',' '),sep=' ')
        if len(v)!=columns:raise ValueError('边界向量/标量维度错误')
        a=np.tile(v,(count,1))
    else:
        m=re.search(r'\b'+key+r'\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
        if not m:return None
        if int(m.group(2))!=count or m.group(1)!=('vector' if vector else 'scalar'):raise ValueError('边界数量/类型错误')
        v=np.fromstring(m.group(3).replace('(',' ').replace(')',' '),sep=' ')
        if len(v)!=count*columns:raise ValueError('边界数据不完整')
        a=v.reshape(count,columns)
    if not np.isfinite(a).all():raise ValueError('边界非有限')
    return a if vector else a[:,0]


def geometry(case):
    """模块二：从原连接独立计算向量面积/体积，不用显示网格或近邻映射。"""
    case=Path(case);data=np.load(case/'exact_mesh.npz');g={k:data[k] for k in data.files};p=g['points'];c=g['cells'];xyz=p[g['faces']]
    g['sf']=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])/2;g['area']=np.linalg.norm(g['sf'],axis=1)
    g['volume']=abs(np.linalg.det((p[c][:,1:]-p[c][:,:1]).transpose(0,2,1)))/6
    g['patches']=json.loads((case/'patches.json').read_text());g['centers']=p[c].mean(axis=1)
    if np.any(g['volume']<=0) or np.any(g['area']<=0):raise ValueError('几何非正')
    return g


def boundary_values(directory,name,internal,g):
    """模块三：实际写出的边界优先；零梯度/滑移用其精确原生代数定义。"""
    text=(Path(directory)/name).read_text(encoding='utf-8');text=text[text.index('boundaryField'):]
    ni=len(g['neighbor']);n=len(g['owner'])-ni;vector=name=='U';values=np.zeros((n,3)) if vector else np.zeros(n)
    for patch in g['patches']:
        start=patch['startFace'];count=patch['nFaces'];ids=np.arange(start,start+count);out=ids-ni
        m=re.search(r'\b'+re.escape(patch['name'])+r'\s*\{([^}]+)\}',text,re.S)
        if not m:raise ValueError('缺少边界 '+patch['name'])
        body=m.group(1);kind=re.search(r'\btype\s+(\w+);',body).group(1)
        a=array_entry(body,'value',count,vector)
        if a is None:
            if kind=='zeroGradient':a=internal[g['owner'][ids]]
            elif kind=='slip' and vector:
                v=internal[g['owner'][ids]];normal=g['sf'][ids]/g['area'][ids,None]
                a=v-np.einsum('ij,ij->i',v,normal)[:,None]*normal
            else:raise ValueError('没有可证明的实际边界值 '+name+'/'+patch['name']+'/'+kind)
        values[out]=a
    return values


def conservative(fields):
    rho=fields['rho'];u=fields['U'];total=fields['e']+np.einsum('ij,ij->i',u,u)/2
    return np.column_stack([rho,rho[:,None]*u,rho*total])


def flux_and_residual(directory,g):
    """模块四：同EOS完整两侧状态→一个共享面通量→只读空间残差。"""
    f=read_state(directory,len(g['cells']));audit_state(f);owner=g['owner'];nei=g['neighbor'];ni=len(nei)
    side={}
    for name in ['rho','U','T','e']:
        boundary=boundary_values(directory,name,f[name],g)
        left=np.concatenate([f[name][owner[:ni]],boundary]);right=np.concatenate([f[name][nei],boundary]);side[name]=(left,right)
    rhoL,rhoR=side['rho'];uL,uR=side['U'];tL,tR=side['T'];eL,eR=side['e']
    pL=rhoL*R*tL;pR=rhoR*R*tR
    aL=evaluate(pL,tL)['soundSpeedMps'];aR=evaluate(pR,tR)['soundSpeedMps']
    sf=g['sf'];mag=g['area'];phiL=np.einsum('ij,ij->i',uL,sf);phiR=np.einsum('ij,ij->i',uR,sf)
    ap=np.maximum.reduce([phiL+aL*mag,phiR+aR*mag,np.zeros_like(mag)])
    am=np.minimum.reduce([phiL-aL*mag,phiR-aR*mag,np.zeros_like(mag)])
    weight=ap/(ap-am);asf=am*weight;alpha=weight*phiL-asf;beta=(1-weight)*phiR+asf
    mass=alpha*rhoL+beta*rhoR
    momentum=(alpha*rhoL)[:,None]*uL+(beta*rhoR)[:,None]*uR+(weight*pL+(1-weight)*pR)[:,None]*sf
    EL=eL+np.einsum('ij,ij->i',uL,uL)/2;ER=eR+np.einsum('ij,ij->i',uR,uR)/2
    energy=alpha*(rhoL*EL+pL)+beta*(rhoR*ER+pR)+asf*(pL-pR)
    flux=np.column_stack([mass,momentum,energy]);total=np.zeros((len(g['cells']),5))
    np.add.at(total,owner,flux);np.add.at(total,nei,-flux[:ni]);residual=total/g['volume'][:,None]
    return f,flux,residual


def reference_scales():
    gas=evaluate(5000,300);rho=float(gas['densityKgm3']);speed=3*float(gas['soundSpeedMps'])
    scale=np.array([rho,*([rho*speed]*3),rho*(abs(float(gas['internalEnergyJkg']))+speed**2)])
    return scale,speed


def native_mass_flux(directory,g):
    internal=internal_field(Path(directory)/'phi',len(g['neighbor']))
    # phi属于surface场，所有patch均必须有实际value，不用cell owner作代填。
    return np.r_[internal,boundary_values(directory,'phi',internal,g)]
