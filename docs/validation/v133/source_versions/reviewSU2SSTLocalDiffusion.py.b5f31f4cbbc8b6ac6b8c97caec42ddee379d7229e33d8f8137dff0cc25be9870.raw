"""v132只读重建尾缘局部扩散：分离两点、非正交和交叉项，不冒充原历史增量。"""
from pathlib import Path
import json,time,zipfile,hashlib
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha
from su2SSTInterface import state,triangle_distances
BASE=ROOT/'outputs/v132_validation'
FACES={6:[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]],
       4:[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]}
EDGES={n:sorted({tuple(sorted((a,b))) for f in faces for a,b in zip(f,f[1:]+f[:1])}) for n,faces in FACES.items()}
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def selected_cells(cells,selected,n):
    mark=np.zeros(n,bool);mark[list(selected)]=True
    return cells[np.any(mark[cells],axis=1)]
def neighbors_of(cells,selected,n):
    result={int(i):set() for i in selected}
    for group in cells:
        for row in selected_cells(group,selected,n):
            for a,b in EDGES[len(row)]:
                i,j=int(row[a]),int(row[b])
                if i in result:result[i].add(j)
                if j in result:result[j].add(i)
    return result
def local_geometry(points,cells,core):
    """按原生中位对偶的边中点—面心—单元心三角片重建面积向量和正子体积。"""
    volume={i:0. for i in core};normal={}
    for group in cells:
        for row in selected_cells(group,core,len(points)):
            q=points[row];center=q.mean(axis=0)
            for face in FACES[len(row)]:
                fc=q[face].mean(axis=0)
                for a,b in zip(face,face[1:]+face[:1]):
                    ia,ib=int(row[a]),int(row[b])
                    if ia not in core and ib not in core:continue
                    lo,hi=sorted((ia,ib));mid=(points[lo]+points[hi])/2
                    av=.5*np.cross(center-mid,fc-mid);dot=float(av@(points[hi]-points[lo]))
                    assert dot!=0
                    av*=1 if dot>0 else -1
                    normal[(lo,hi)]=normal.get((lo,hi),np.zeros(3))+av
                    for i in [lo,hi]:
                        if i in core:
                            sub=abs(np.cross(mid-points[i],fc-points[i])@(center-points[i]))/6
                            assert sub>0;volume[i]+=sub
    adj={i:[] for i in core}
    for (i,j),v in normal.items():
        if i in core:adj[i].append((j,v))
        if j in core:adj[j].append((i,-v))
    result={}
    for i in core:
        js=np.array([j for j,v in adj[i]],int);ns=np.array([v for j,v in adj[i]])
        e=points[js]-points[i];w=1/np.sum(e*e,axis=1);A=(e*w[:,None]).T@e
        condition=float(np.linalg.cond(A));inverse=np.linalg.inv(A)
        moment=inverse@((e*w[:,None]).T@e)
        closure=float(np.linalg.norm(ns.sum(axis=0))/np.linalg.norm(ns,axis=1).sum())
        gg=e.T@ns/(2*volume[i])
        result[i]=dict(ids=js,normals=ns,volume=volume[i],edges=e,inverse=inverse,weights=w,
            condition=condition,closure=closure,linearError=float(np.max(abs(moment-np.eye(3)))),
            ggLinearError=float(np.max(abs(gg-np.eye(3)))))
    return result
def gradients(g,values):
    """等价使用差值形式减少常量抵消；未验证为原生历史梯度逐字重放。"""
    result={}
    for i,row in g.items():
        delta=values[row['ids']]-values[i]
        ls=(row['inverse']@((row['edges']*row['weights'][:,None]).T@delta)).T
        gg=delta.T@row['normals']/(2*row['volume'])
        result[i]=dict(wls=ls,gg=gg)
    return result
def projections(e,n,delta,avg):
    h=float(n@e/(e@e));orth=delta*h
    correction=avg@(n-h*e)
    return orth,correction
def main():
    started=time.monotonic();assert not BASE.exists();BASE.mkdir()
    meshfile=ROOT/'outputs/v120_validation/mesh32/hybrid_volume.mat'
    d=loadmat(str(meshfile),squeeze_me=True);P=d['vertices'];F=d['sourceFaces'].astype(int)-1
    cells=[d['wedges'].astype(int)-1,d['tetrahedra'].astype(int)-1];nw=9878
    assert len(P)==349206 and len(F)==19752 and len(d['sourceVertices'])==nw
    targets=[29650,29651,29652,39527,39528,39529,39530,148188,148189,158067]
    controls=[layer*nw+parent for layer in [3,4,15,16] for parent in [0,6,23,30,38,45]]
    focus=set(targets+controls);adj=neighbors_of(cells,focus,len(P))
    core=focus|set(j for v in adj.values() for j in v)
    assert min(core)>=nw and max(core)<32*nw
    g=local_geometry(P,cells,core)
    assert all(g[i]['closure']<=1e-10 and g[i]['linearError']<=1e-10 for i in core)
    # 完整三角搜索半径上界，保证包含最短壁距候选，而非只取父节点距离。
    tri=P[F];centers=tri.mean(axis=1);radius=np.linalg.norm(tri-centers[:,None],axis=2).max()
    vt=cKDTree(P[:nw]);ct=cKDTree(centers);distance={}
    for i in core:
        upper=vt.query(P[i])[0];near=ct.query_ball_point(P[i],radius+upper+1e-12)
        distance[i]=float(min(triangle_distances(P[i],tri[near,0],tri[near,1],tri[near,2])))
        assert distance[i]>0
    srcbase=ROOT/'tools/runtime/su2-thermal-build/source'
    paths=['Common/src/geometry/dual_grid/CEdge.cpp','SU2_CFD/include/gradients/computeGradientsLeastSquares.hpp',
        'SU2_CFD/include/gradients/computeGradientsGreenGauss.hpp','SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp',
        'SU2_CFD/include/numerics/turbulent/turb_sst_edge_flux.hpp','SU2_CFD/src/variables/CTurbSSTVariable.cpp']
    source={}
    with zipfile.ZipFile(str(ROOT/'tools/runtime/su2-thermal-build/su2.zip')) as z:
        for rel in paths:
            raw=z.read(next(n for n in z.namelist() if n.endswith('/'+rel)))
            assert (srcbase/rel).read_bytes()==raw
            source[str(srcbase/rel)]=hashlib.sha256(raw).hexdigest()
    report=dict(targetNodeIDs=targets,controlNodeIDs=controls,gradientNodes=len(core),
        sourceSHA256=source,geometry={},variants={},nativeRuns=0,originalHistoricalRHSReproduced=False,
        originalHistoricalIncrementKnown=False,physicsQualified=False,newQualifiedArea=0)
    for i in sorted(focus):
        row=g[i];e=row['edges'];n=row['normals']
        cosine=np.sum(e*n,axis=1)/(np.linalg.norm(e,axis=1)*np.linalg.norm(n,axis=1))
        wall=i%nw;faces=np.where(np.any(F==wall,axis=1))[0];q=P[F[faces]]
        normals=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
        report['geometry'][str(i)]=dict(xyzM=P[i].tolist(),wallParentID=wall,normalLayer=i//nw,
            wallDistanceM=distance[i],volumeM3=row['volume'],wlsCondition=row['condition'],
            wlsLinearGeometryError=row['linearError'],ggLinearGeometryError=row['ggLinearError'],
            dualAreaClosureRelativeError=row['closure'],minimumDualAlignmentCosine=float(min(cosine)),
            maximumIncidentWallNormalAngleDegrees=float(np.rad2deg(np.arccos(np.clip((normals@normals.T).min(),-1,1)))))
    for name,folder in [('official_v130',ROOT/'outputs/v130_validation/column_final'),('bounded_v131',ROOT/'outputs/v131_validation/actual')]:
        receipt=read(folder/'engine_process.json');path=folder/'restart_flow.csv'
        assert receipt['exitCode']==0 and all(receipt['inputsStable'].values())
        assert sha(path)==receipt['outputSHA256']['restart_flow.csv']
        a,keys,nd,rho,uv,energy=state(path,True,True);assert nd==3
        assert np.max(abs(a[:,1:4]-P))<1e-12
        kw=a[:,[keys.index('Turb_Kin_Energy'),keys.index('Omega')]]
        mu=a[:,keys.index('Laminar_Viscosity')];mut=a[:,keys.index('Eddy_Viscosity')]
        grad=gradients(g,kw);F1={}
        for i in core:
            k,w=kw[i];y=distance[i];cd=max(2*rho[i]*.856/w*(grad[i]['wls'][0]@grad[i]['wls'][1]),1e-10)
            arg=min(max(np.sqrt(k)/(.09*w*y+1e-32),500*mu[i]/(rho[i]*y*y*w+1e-32)),4*rho[i]*.856*k/(cd*y*y+1e-32))
            F1[i]=float(np.tanh(arg**4))
        rows={}
        for i in sorted(focus):
            row=g[i];js=row['ids'];sigma_i=.5*F1[i]+.856*(1-F1[i])
            result=dict(k=float(kw[i,0]),omega=float(kw[i,1]),F1FromSavedStateWLS=F1[i],
                localOmegaMinimum=bool(kw[i,1]<=min(kw[js,1])*(1+1e-12)),
                minimumNeighborOmega=float(min(kw[js,1])),maximumNeighborOmega=float(max(kw[js,1])),
                frozenSavedMu= float(mu[i]),frozenSavedMuT=float(mut[i]),methods={})
            for method in ['wls','gg']:
                edges=[]
                for j,e,n in zip(js,row['edges'],row['normals']):
                    j=int(j);sigma_j=.5*F1[j]+.856*(1-F1[j])
                    D=.5*(mu[i]+sigma_i*mut[i]+mu[j]+sigma_j*mut[j]);assert D>0
                    lam=(1-F1[i])*rho[i]*sigma_i+(1-F1[j])*rho[j]*sigma_j
                    cross=lam*(1-2*kw[i,1]/(kw[i,1]+kw[j,1]))
                    orth,corr=projections(e,n,kw[j]-kw[i],.5*(grad[i][method]+grad[j][method]))
                    diag=D*(orth[1]+corr[1]);crossFlux=cross*(orth[0]+corr[0])
                    edges.append(dict(neighbor=j,orthogonalOmegaDiffusion=float(D*orth[1]),
                        nonorthogonalOmegaCorrection=float(D*corr[1]),totalDiagonalOmegaDiffusion=float(diag),
                        crossDiffusion=float(crossFlux),totalOmegaDiffusion=float(diag+crossFlux),
                        neighborOmega=float(kw[j,1]),neighborK=float(kw[j,0])))
                totals={key:float(sum(e[key] for e in edges)) for key in ['orthogonalOmegaDiffusion','nonorthogonalOmegaCorrection',
                    'totalDiagonalOmegaDiffusion','crossDiffusion','totalOmegaDiffusion']}
                result['methods'][method]=dict(totals=totals,edges=edges)
            rows[str(i)]=result
        report['variants'][name]=dict(nodes=rows,coefficientConvention='saved mu/muT frozen; F1 recomputed from saved k/omega using WLS for both gradient comparisons; not original step reconstruction',
            sourceLocalIteration=99,sourceSHA256={str(folder/'engine_process.json'):sha(folder/'engine_process.json'),str(path):sha(path)})
    report['elapsedSeconds']=time.monotonic()-started
    report['sourceSHA256'].update({str(p):sha(p) for p in [Path(__file__).resolve(),meshfile,ROOT/'docs/validation/v132/PLAN.md']})
    report['summary']={}
    for name,variant in report['variants'].items():
        minima=[variant['nodes'][str(i)] for i in targets if variant['nodes'][str(i)]['localOmegaMinimum']]
        report['summary'][name]=dict(targetLocalMinima=len(minima),methods={method:dict(
            negativeDiagonalDiffusionCount=sum(r['methods'][method]['totals']['totalDiagonalOmegaDiffusion']<0 for r in minima),
            negativeTotalDiffusionCount=sum(r['methods'][method]['totals']['totalOmegaDiffusion']<0 for r in minima),
            diagonalDiffusionSum=sum(r['methods'][method]['totals']['totalDiagonalOmegaDiffusion'] for r in minima),
            crossDiffusionSum=sum(r['methods'][method]['totals']['crossDiffusion'] for r in minima)) for method in ['wls','gg']})
    dump(BASE/'local_diffusion.json',report)
    print(json.dumps(dict(elapsedSeconds=report['elapsedSeconds'],gradientNodes=len(core),summary=report['summary'],
        maximumGGLinearError=max(v['ggLinearGeometryError'] for v in report['geometry'].values())),indent=2))
if __name__=='__main__':main()
