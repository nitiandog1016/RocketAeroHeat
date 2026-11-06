"""一般混合域的实际三角壁距初猜；不按层号推坐标，不复制旧流场。"""
from pathlib import Path
import argparse,json,time,re
import numpy as np
from prepareSU2CheckpointContinuation import ROOT,sha,read,config
from prepareSU2SSTCase import prepare as prepare_sst,set_value
from su2SSTInterface import state,wall_geometry,mesh

KIND='general-mesh-SST-wall-distance-initial-guess'

def project_triangle(p,t):
    """模块一：精确面内/边/顶点投影，返回凸重心权重；不用最近顶点近似壁距。"""
    a,b,c=t;ab=b-a;ac=c-a;ap=p-a
    aa=ab@ab;bb=ab@ac;cc=ac@ac;den=aa*cc-bb*bb;assert den>0
    s=(ap@ab*cc-ap@ac*bb)/den;v=(ap@ac*aa-ap@ab*bb)/den
    weights=np.column_stack([1-s-v,s,v]);inside=np.all(weights>=0,axis=1)
    dist=np.full(len(p),np.inf);best=np.zeros((len(p),3))
    dist[inside]=np.sum((p[inside]-weights[inside]@t)**2,axis=1);best[inside]=weights[inside]
    for i,j in [(0,1),(1,2),(2,0)]:
        e=t[j]-t[i];f=np.clip((p-t[i])@e/(e@e),0,1);q=t[i]+f[:,None]*e
        d=np.sum((p-q)**2,axis=1);take=d<dist
        dist[take]=d[take];best[take]=0;best[take,i]=1-f[take];best[take,j]=f[take]
    return dist,best

def surface_distance(points,vertices,faces,limit):
    """模块二：确定性AABB树精确搜索；超过初猜厚度只需证明d>=limit。"""
    assert np.isfinite(points).all() and np.isfinite(vertices).all() and np.isfinite(limit) and limit>0
    tri=vertices[faces];lo=tri.min(axis=1);hi=tri.max(axis=1);center=tri.mean(axis=1)
    def build(ids):
        lower=lo[ids].min(axis=0);upper=hi[ids].max(axis=0)
        if len(ids)<=12:return lower,upper,ids,None
        axis=np.argmax(np.ptp(center[ids],axis=0));order=ids[np.argsort(center[ids,axis],kind='mergesort')];mid=len(ids)//2
        return lower,upper,None,(build(order[:mid]),build(order[mid:]))
    tree=build(np.arange(len(faces)));distance=np.full(len(points),limit**2);which=np.full(len(points),-1,np.int64)
    weights=np.zeros((len(points),3));pairs=0;pad=64*np.finfo(float).eps*max(1.,float(abs(vertices).max()),limit)
    def visit(node,rows):
        nonlocal pairs
        lower,upper,ids,children=node
        gap=np.maximum(np.maximum(lower-points[rows],points[rows]-upper),0)
        rows=rows[np.sum(gap*gap,axis=1)<=(np.sqrt(distance[rows])+pad)**2]
        if len(rows)==0:return
        if ids is not None:
            for i in ids:
                d,w=project_triangle(points[rows],tri[i]);pairs+=len(rows);take=d<distance[rows];sel=rows[take]
                distance[sel]=d[take];which[sel]=i;weights[sel]=w[take]
        else:
            visit(children[0],rows);visit(children[1],rows)
    visit(tree,np.arange(len(points)))
    return np.sqrt(distance),which,weights,pairs

def initial_values(points,vertices,faces,wall_ids,wall_lengths,far_ids,thickness,U,pinf,Tinf,Tw,R,gamma,kinf,winf,nu_wall):
    """模块三：正温/密度、无滑移、精确远场和含k能量；这是独立初猜而非物理解。"""
    d,face,w,pairs=surface_distance(points,vertices,faces,thickness)
    assert np.all(d[wall_ids]<1e-12) and np.all(face[wall_ids]>=0)
    d[wall_ids]=0.;eta=d/thickness;b=eta**2*(3-2*eta);assert np.all(b[far_ids]==1)
    cp=gamma*R/(gamma-1);speed2=float(np.dot(U,U));velocity=b[:,None]*U
    T=Tw+(Tinf-Tw)*b+speed2/(2*cp)*b*(1-b);rho=pinf/(R*T);k=kinf*b
    omega=np.full(len(points),winf);inside=face>=0
    lengths=np.full(len(vertices),np.nan);lengths[wall_ids]=wall_lengths
    ell=np.sum(lengths[faces[face[inside]]]*w[inside],axis=1);assert np.all(ell>0)
    near=60*nu_wall/(.075*(d[inside]+ell)**2)
    omega[inside]=np.exp((1-b[inside])*np.log(near)+b[inside]*np.log(winf))
    omega[wall_ids]=60*nu_wall/(.075*wall_lengths**2)
    energy=rho*(R*T/(gamma-1)+.5*np.sum(velocity**2,axis=1)+k)
    fields=np.column_stack([rho,rho[:,None]*velocity,energy,k,omega])
    assert np.isfinite(fields).all() and min(T)>0 and min(rho)>0 and min(omega)>0 and min(k)>=0
    assert np.all(velocity[wall_ids]==0) and np.all(T[wall_ids]==Tw) and np.all(k[wall_ids]==0)
    return fields,dict(distance=d,face=face,weights=w,blend=b,temperature=T,pairs=pairs)

def boundary_faces(path):
    """读取具名壁面/远场，节点坐标与SU2实际网格共同作为初态输入。"""
    with Path(path).open() as f:
        for line in f:
            if line.startswith('NMARK='):break
        count=int(line.split('=')[1]);out={}
        for _ in range(count):
            name=f.readline().split('=')[1].strip();n=int(f.readline().split('=')[1])
            a=np.array([list(map(int,f.readline().split())) for _ in range(n)],np.int64);assert np.all(a[:,0]==5)
            out[name]=a[:,1:]
    assert set(out)=={'wall','farfield'};return out

def prepare(folder,thickness,intensity,viscosity_ratio,prt,provenance):
    start=time.monotonic();folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    assert not (folder/'initial_flow.csv').exists() and not (folder/'engine_process.json').exists()
    prepare_sst(folder,None,intensity,viscosity_ratio,prt,provenance)
    c=config(folder/'case.cfg');manifest=read(folder/'case_manifest.json');contract=read(folder/'sst_contract.json')
    assert manifest['thermalBoundary']['mode']=='isothermal'
    g=wall_geometry(folder/'mesh.su2');p=g['points'];assert p.shape[1]==3
    boundary=boundary_faces(folder/'mesh.su2');F=boundary['wall'];far=np.unique(boundary['farfield']);ids=g['wallIDs']
    R=float(c['GAS_CONSTANT']);gamma=float(c['GAMMA_VALUE']);pinf=float(c['FREESTREAM_PRESSURE']);Tinf=float(c['FREESTREAM_TEMPERATURE'])
    Tw=float(re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)',c['MARKER_ISOTHERMAL'])[1])
    mu=float(c['MU_REF'])*(Tw/float(c['MU_T_REF']))**1.5*(float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(Tw+float(c['SUTHERLAND_CONSTANT']))
    U=np.asarray(manifest['expectedVelocity'],float).ravel();assert abs(np.linalg.norm(U)/(float(c['MACH_NUMBER'])*np.sqrt(gamma*R*Tinf))-1)<1e-12
    fields,detail=initial_values(p,p,F,ids,g['nearestWallDistance'],far,thickness,U,pinf,Tinf,Tw,R,gamma,contract['kInfinity'],contract['omegaInfinity'],mu/(pinf/(R*Tw)))
    keys=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy','Turb_Kin_Energy','Omega']
    path=folder/'initial_flow.csv';np.savetxt(str(path),np.column_stack([np.arange(len(p)),p,fields]),delimiter=',',fmt='%.17g',header=','.join('"'+s+'"' for s in keys),comments='')
    _,_,_,rho,u,e=state(path,True);err=float(max(abs(e*(gamma-1)/R-detail['temperature'])/detail['temperature']));assert err<1e-12
    geometry=folder/'wall_distance_initialization.npz'
    np.savez_compressed(str(geometry),distance=detail['distance'],nearestFace=detail['face'],weights=detail['weights'],wallIDs=ids,wallLengths=g['nearestWallDistance'],farIDs=far)
    locks=[Path(__file__).resolve(),ROOT/'scripts/su2SSTInterface.py',ROOT/'scripts/prepareSU2SSTCase.py',folder/'mesh.su2',geometry]
    init=dict(kind=KIND,method='exact triangular wall distance; C1 velocity/total enthalpy; positive SST initial guess',
        thicknessM=thickness,nodes=len(p),wallNodes=len(ids),localIterationStartsAt=0,sourceLocalIteration=None,
        sourceDirectory=None,initialGuessOnly=True,oldFlowInterpolated=False,physicsQualified=False,
        initialStateSHA256=sha(path),sourceStateHasNativePrimitive=False,thermalEnergyRelativeError=err,
        distanceTrianglePairs=detail['pairs'],temperatureRangeK=[float(min(detail['temperature'])),float(max(detail['temperature']))],
        nearestWallDistanceRangeM=[float(min(g['nearestWallDistance'])),float(max(g['nearestWallDistance']))],
        exactFarfield=True,wallNoSlip=True,givenWallTemperatureK=Tw,preparationSeconds=time.monotonic()-start,
        sourceSHA256={str(q):sha(q) for q in locks})
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n')
    text=(folder/'case.cfg').read_text()
    for key,val in dict(RESTART_SOL='YES',SOLUTION_FILENAME='initial_flow',READ_BINARY_RESTART='NO').items():text=set_value(text,key,val)
    (folder/'case.cfg').write_text(text)
    contract.pop('initialization',None);contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'))
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    manifest.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');return init

def validate_source(folder,c,init):
    """模块四：运行前核对初态全部坐标、热状态和边界，不接受仅有标签的重启动。"""
    assert init['kind']==KIND and init['localIterationStartsAt']==0 and init['sourceLocalIteration'] is None
    assert init['initialGuessOnly'] and not init['oldFlowInterpolated'] and not init['physicsQualified']
    for n,h in init['sourceSHA256'].items():assert sha(Path(n))==h
    assert sha(folder/'initial_flow.csv')==init['initialStateSHA256']
    a,keys,nd,rho,u,e=state(folder/'initial_flow.csv',True);assert nd==3 and len(a)==init['nodes']
    p,_,wall=mesh(folder/'mesh.su2');assert np.array_equal(a[:,1:4],p)
    d=np.load(str(folder/'wall_distance_initialization.npz'));b=(d['distance']/init['thicknessM'])**2*(3-2*d['distance']/init['thicknessM'])
    manifest=read(folder/'case_manifest.json');U=np.asarray(manifest['expectedVelocity']).ravel();R=float(c['GAS_CONSTANT']);gamma=float(c['GAMMA_VALUE'])
    T=init['givenWallTemperatureK']+(float(c['FREESTREAM_TEMPERATURE'])-init['givenWallTemperatureK'])*b+float(U@U)/(2*gamma*R/(gamma-1))*b*(1-b)
    assert max(abs(e*(gamma-1)/R-T)/T)<1e-12 and max(abs(rho*R*T/float(c['FREESTREAM_PRESSURE'])-1))<1e-12
    assert np.allclose(u,b[:,None]*U,rtol=2e-15,atol=1e-14)
    assert np.all(a[d['wallIDs'],keys.index('Turb_Kin_Energy')]==0) and np.all(b[d['farIDs']]==1)
    model=manifest['sst'];k=a[:,keys.index('Turb_Kin_Energy')];omega=a[:,keys.index('Omega')]
    assert np.allclose(k,model['kInfinity']*b,rtol=2e-15,atol=0)
    F=np.array([row[1:] for row in wall],np.int64);inside=d['nearestFace']>=0
    lengths=np.full(len(p),np.nan);lengths[d['wallIDs']]=d['wallLengths']
    ell=np.sum(lengths[F[d['nearestFace'][inside]]]*d['weights'][inside],axis=1)
    Tw=init['givenWallTemperatureK'];mu=float(c['MU_REF'])*(Tw/float(c['MU_T_REF']))**1.5*(float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(Tw+float(c['SUTHERLAND_CONSTANT']))
    nu=mu/(float(c['FREESTREAM_PRESSURE'])/(R*Tw));expected=np.full(len(p),model['omegaInfinity'])
    expected[inside]=np.exp((1-b[inside])*np.log(60*nu/(.075*(d['distance'][inside]+ell)**2))+b[inside]*np.log(model['omegaInfinity']))
    expected[d['wallIDs']]=60*nu/(.075*d['wallLengths']**2)
    assert max(abs(omega/expected-1))<1e-12
    return init['sourceSHA256']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--thickness',type=float,required=True)
    p.add_argument('--intensity',type=float,required=True);p.add_argument('--viscosity-ratio',type=float,required=True)
    p.add_argument('--prt',type=float,required=True);p.add_argument('--provenance',required=True);a=p.parse_args()
    print(json.dumps(prepare(a.folder,a.thickness,a.intensity,a.viscosity_ratio,a.prt,a.provenance),indent=2))
