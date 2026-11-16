"""模块职责：从原箭体导出主导残差点的真实单元星与状态，用于具名成熟方法控制。"""
from pathlib import Path
import csv,struct,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from fasRecoveryFieldsV248 import state,residual,coarse

BASE=ROOT/'outputs/v249_validation';DOC=ROOT/'docs/validation/v249';CLOSED=False
PREV=ROOT/'outputs/v248_validation';CASE=PREV/'rocket01'
MG=ROOT/'outputs/v246_validation/initial01'
EDGES={10:((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)),
       13:((0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5))}


def extract(path,seeds,cell_count=705490,point_count=335726):
    """模块一：读取原SU2文本；只保存触及种子的母单元，保留真实边而非单元完全图。"""
    selected=[];markers={};stars={int(i):[] for i in seeds};neighbors={int(i):set() for i in seeds}
    seeds=set(stars)
    with Path(path).open(encoding='utf8') as f:
        assert next(f).strip()=='NDIME= 3';assert int(next(f).split('=')[1])==cell_count
        for cid in range(cell_count):
            row=list(map(int,next(f).split()));kind=row[0];nodes=row[1:];assert kind in EDGES
            assert len(nodes)==(4 if kind==10 else 6) and len(set(nodes))==len(nodes)
            touching=seeds.intersection(nodes)
            if touching:
                selected.append([cid,kind]+nodes+[-1]*(6-len(nodes)))
                for seed in touching:stars[seed].append(cid)
                for u,v in EDGES[kind]:
                    if nodes[u] in touching:neighbors[nodes[u]].add(nodes[v])
                    if nodes[v] in touching:neighbors[nodes[v]].add(nodes[u])
        assert int(next(f).split('=')[1])==point_count;xyz=np.empty((point_count,3))
        for i in range(point_count):
            row=next(f).split();assert int(row[3])==i;xyz[i]=list(map(float,row[:3]))
        for _ in range(int(next(f).split('=')[1])):
            name=next(f).split('=')[1].strip();count=int(next(f).split('=')[1]);marker=set()
            for _ in range(count):
                row=list(map(int,next(f).split()));assert row[0]==5;marker.update(row[1:])
            markers[name]=marker
        assert not f.read().strip()
    assert set(markers)=={'wall','farfield'} and all(stars.values())
    return xyz,np.asarray(selected,dtype=np.int64),stars,neighbors,markers


def prepare():
    out=BASE/'stencil01';out.mkdir(parents=True,exist_ok=False)
    diagnosis=read(PREV/'diagnosis01/report.json');audit=read(PREV/'audit02/report.json')
    assert diagnosis['passed'] and audit['passed'];seeds=[r['globalID'] for r in diagnosis['topEnergyRHSPoints']]
    assert len(seeds)==len(set(seeds))==20 and seeds[0]==333260
    xyz,cells,stars,neighbors,markers=extract(CASE/'mesh.su2',seeds)
    current=state(CASE/'state_000200.bin');global_to_local=np.argsort(current[:,0]);a=current[global_to_local]
    assert np.array_equal(a[:,1:4],xyz)
    raw,r=residual(CASE/'terminal_spatial.bin');r=r[np.argsort(raw[:,0])]
    with (MG/'children.bin').open('rb') as f:
        assert struct.unpack('<8sQQ',f.read(24))==(b'RAH246ID',335726,2)
        pairs=np.fromfile(f,dtype='<u8').reshape(-1,2).astype(np.int64)
    parents_local=np.empty(335726,dtype=np.int64);parents_local[pairs[:,1]]=pairs[:,0]
    parents=parents_local[global_to_local];ca=coarse(CASE/'coarse_final.bin')
    boundary=np.loadtxt(str(MG/'boundary.csv'),delimiter=',',skiprows=1)
    nativefine=boundary[boundary[:,0]==0];nativecoarse=boundary[boundary[:,0]==1]
    assert set(nativefine[nativefine[:,15]==1,4].astype(int))==markers['wall']
    assert set(nativefine[nativefine[:,15]==0,4].astype(int))==markers['farfield']
    ids=np.unique(cells[:,2:][cells[:,2:]>=0]);selected=a[ids]
    history=[]
    for step in (0,1,50,100,150,200):
        f=state(CASE/('state_%06d.bin'%step),initial=step==0);f=f[np.argsort(f[:,0])]
        history.append(f[ids])
    records=[]
    for seed in seeds:
        nb=np.array(sorted(neighbors[seed]),dtype=np.int64);dx=xyz[nb]-xyz[seed];dist=np.linalg.norm(dx,axis=1)
        assert np.all(dist>0);design=dx/dist[:,None];sing=np.linalg.svd(design,compute_uv=False)
        cs=cells[np.isin(cells[:,0],stars[seed])];tet=cs[cs[:,1]==10,2:6]
        tetvol=abs(np.linalg.det(xyz[tet[:,1:]]-xyz[tet[:,0],None]))/6 if len(tet) else np.array([])
        parent=int(parents[seed]);parentrow=ca[parent]
        br=nativefine[nativefine[:,4]==seed];cp=nativecoarse[nativecoarse[:,3]==parent]
        records.append(dict(globalID=seed,localID=int(global_to_local[seed]),xyzM=xyz[seed].tolist(),
          boundaryLabels=[k for k,v in markers.items() if seed in v],coarseLocalParent=parent,
          coarseGlobalRepresentative=int(parentrow[0]),coarseBoundaryKinds=cp[:,2].astype(int).tolist(),
          nativeFineBoundaryRows=br.tolist(),cellIDs=stars[seed],edgeNeighbors=nb.tolist(),
          incidentCellTypes=np.unique(cs[:,1]).tolist(),minimumTetVolumeM3=float(min(tetvol)) if len(tetvol) else None,
          edgeLengthRangeM=[float(min(dist)),float(max(dist))],normalizedWLSDesignCondition=float(sing[0]/sing[-1]),
          conditionIsGeometricNotNativeGradientError=True,
          TRangeK=[float(min(a[np.r_[seed,nb],11])),float(max(a[np.r_[seed,nb],11]))],
          pRangePa=[float(min(a[np.r_[seed,nb],15])),float(max(a[np.r_[seed,nb],15]))],
          machRange=[float(min(np.linalg.norm(a[np.r_[seed,nb],12:15],axis=1)/a[np.r_[seed,nb],18])),
                     float(max(np.linalg.norm(a[np.r_[seed,nb],12:15],axis=1)/a[np.r_[seed,nb],18)))],
          physicalRHS=r[seed].tolist(),parentChildCount=int(sum(parents==parent)),
          sampledHistory=[dict(update=s,T=float(history[k][np.searchsorted(ids,seed),11]),p=float(history[k][np.searchsorted(ids,seed),15])) for k,s in enumerate((0,1,50,100,150,200))]))
    np.savez_compressed(str(out/'actual_stencils.npz'),globalIDs=ids,xyz=xyz[ids],motherCells=cells,
      fineState=selected,physicalRHS=r[ids],coarseParentLocal=parents[ids],history=np.asarray(history),seedIDs=np.asarray(seeds))
    save(out/'report.json',dict(passed=True,seeds=records,seedCount=20,patchNodes=len(ids),motherCellCount=len(cells),
      originalWallFaces=19752,originalAreaM2=1.9646292565715553,originalMeshSHA256=sha(CASE/'mesh.su2'),
      artifactSHA256=sha(out/'actual_stencils.npz'),nativeBoundaryMembershipMatched=True,
      fineStateAndMeshCoordinatesBitEqual=True,physicalCFDSubdomain=False,newFlowUpdates=0,physicalQualification=False))
    print(dict(passed=True,seeds=20,nodes=len(ids),cells=len(cells),leading=records[0]))


def setup_guard():
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;original=guard.size_bytes
    guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    return guard


def main():
    if CLOSED:raise SystemExit('v249 stencil closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),ROOT/'scripts/fasRecoveryFieldsV248.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
      CASE/'mesh.su2',CASE/'terminal_spatial.bin',CASE/'coarse_final.bin',MG/'children.bin',MG/'boundary.csv',
      PREV/'audit02/report.json',PREV/'diagnosis01/report.json']+[CASE/('state_%06d.bin'%i) for i in (0,1,50,100,150,200)]
    return setup_guard().supervise('stencil',[sys.executable,str(Path(__file__).resolve()),'task'],60,256,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v249 stencil closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':prepare()
    else:raise SystemExit(0 if main() else 1)
