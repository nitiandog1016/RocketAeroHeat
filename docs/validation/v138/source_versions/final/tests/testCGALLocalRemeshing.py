"""独立核对小型局部重网格：准确导入、固定边界、冻结子域、质量与有孔域拒绝。"""
from pathlib import Path
import itertools,json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareCGALRemeshingRuntime import ROOT,DEST,sha
from buildCGALRemeshingRuntime import run
from cgalTetrahedralContract import condition,faces_and_neighbors
from reviewSU2ExplicitQuality import quality

def orient(P,T):
    T=np.array(T,int);q=P[T];bad=np.linalg.det(np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0],q[:,3]-q[:,0]],axis=2))<0
    T[bad]=T[bad][:,[0,2,1,3]];return T

def write_medit(path,P,T,labels):
    with path.open('w') as f:
        f.write('MeshVersionFormatted 2\nDimension 3\nVertices\n%d\n'%len(P))
        for p in P:f.write('%.17g %.17g %.17g 0\n'%tuple(p))
        f.write('Tetrahedra\n%d\n'%len(T))
        for t,label in zip(T,labels):f.write('%d %d %d %d %d\n'%tuple(list(t+1)+[label]))
        f.write('End\n')

def read_native(path):
    with path.open() as f:
        nv,nt=map(int,f.readline().split());P=np.array([list(map(float,f.readline().split())) for _ in range(nv)])
        T=np.array([list(map(int,f.readline().split())) for _ in range(nt)])
        assert not f.read().strip()
    return P,T[:,:4],T[:,4]

def keys(P,T):return sorted(tuple(sorted(tuple(p) for p in P[t])) for t in T)

def fixtures():
    a=np.array([[0.,0,0],[1,0,0],[.5,np.sqrt(3)/2,0],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
    P=np.r_[a,[[.5,.2,.001]]];T=orient(P,[[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]])
    yield 'convex_interior',P,T,np.ones(4,int),.85
    P=np.r_[a,[[.5,np.sqrt(3)/6,-np.sqrt(2/3)],[.5,.2,.001]]]
    T=orient(P,[[5,1,2,3],[0,5,2,3],[0,1,5,3],[0,1,2,5],[0,2,1,4]])
    yield 'selected_subdomain',P,T,np.array([1,1,1,1,2]),.85

def hollow():
    pts=list(itertools.product(range(4),repeat=3));P=np.array(pts,float);idx={p:i for i,p in enumerate(pts)};T=[]
    for origin in itertools.product(range(3),repeat=3):
        if origin==(1,1,1):continue
        for perm in itertools.permutations(range(3)):
            q=np.array(origin);tet=[idx[tuple(q)]]
            for axis in perm:q=q.copy();q[axis]+=1;tet.append(idx[tuple(q)])
            T.append(tet)
    return P,orient(P,T)

def main():
    out=ROOT/'docs/validation/v137/local_tests';out.mkdir(parents=True,exist_ok=False)
    runtime=json.loads((DEST/'runtime.json').read_text());exe=Path(runtime['executable']);assert sha(exe)==runtime['executableSHA256']
    locked={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__).resolve(),
        ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',
        ROOT/'scripts/buildCGALRemeshingRuntime.py',ROOT/'third_party/cgal/local_tetra_remesh.cpp',exe,DEST/'runtime.json']}
    reports=[]
    for name,P,T,label,h in fixtures():
        folder=out/name;folder.mkdir();inp=folder/'input.mesh';write_medit(inp,P,T,label)
        before,_,_=condition(P,T);assert before['directImportNecessaryConditions']
        r=run([str(exe),str(inp),str(folder/'native'),str(h),'1'],folder/'process',60,
              dict(locked,**{inp.relative_to(ROOT).as_posix():sha(inp)}))
        assert r['exitCode']==0 and all(r['sourceStableAfter'].values()),name
        ip,it,il=read_native(folder/'native_import.txt');np_,nt,nl=read_native(folder/'native_result.txt')
        assert sorted(map(tuple,ip))==sorted(map(tuple,P)) and keys(ip,it)==keys(P,T)
        for sub in np.unique(label):assert keys(ip,it[il==sub])==keys(P,T[label==sub])
        after,_,_=condition(np_,nt);assert after['directImportNecessaryConditions']
        bf,_=faces_and_neighbors(T);af,_=faces_and_neighbors(nt);assert keys(P,bf)==keys(np_,af)
        assert abs(after['tetVolume']/before['tetVolume']-1)<2e-12
        assert keys(P,T[label==2])==keys(np_,nt[nl==2])
        # 边界循环和内面反向独立检查，防止只比排序集合掩盖方向错误。
        allfaces=nt[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
        from auditSU2HybridTopology import paired
        paired(allfaces,af)
        q0=float(min(quality(P[T])));q1=float(min(quality(np_[nt])))
        assert keys(P,T)!=keys(np_,nt) and q1>q0 and q1>=.1
        item=dict(case=name,passed=True,exactConnectivityImport=True,fixedBoundaryExact=True,frozenSubdomainExact=True,
            before=before,after=after,minimumSICNBefore=q0,minimumSICNAfter=q1,nativeSeconds=r['seconds'],
            filesSHA256={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()})
        reports.append(item);print(name,q0,q1,flush=True)
    P,T=hollow();path=out/'hollow_cube.mesh';write_medit(path,P,T,np.ones(len(T),int));h,_,_=condition(P,T)
    assert abs(h['tetVolume']-26)<1e-12 and abs(h['convexHullVolume']-27)<1e-12 and not h['directImportNecessaryConditions']
    reports.append(dict(case='hollow_cube',passed=True,expectedVolume=26,expectedHullVolume=27,
        preflightRejected=True,nativeInvoked=False,report=h,inputSHA256=sha(path)))
    result=dict(passed=True,cases=reports,physicsQualified=False,newQualifiedArea=0,
        lockedSourceSHA256=locked,sourceStableAfter={n:sha(ROOT/n)==h for n,h in locked.items()})
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('three discriminating checks passed')
if __name__=='__main__':main()
