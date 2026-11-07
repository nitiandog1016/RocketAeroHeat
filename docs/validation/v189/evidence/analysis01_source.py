# -*- coding: utf-8 -*-
"""只读原面失效原因，为混合路线选区；不推定分离或补热流。"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analyze(out):
    """模块一：绑定既有工程快照，重算原面积并用独立CSV核对身份。"""
    start=time.monotonic(); out=Path(out); out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'outputs/v104_validation/baseline02'
    paths=[base/'current_case.mat',base/'faces.csv',base/'summary.json',Path(__file__)]
    before={str(p):sha(p) for p in paths}
    result=dict(passed=False,newPhysicalQualification=False,newSolvedAreaM2=0)
    try:
        data=loadmat(paths[0],simplify_cells=True);m=data['current'];v=m['V'];f=m['F'].astype(int)-1
        area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)/2
        centers=v[f].mean(axis=1);n=len(f);code=m['heatingCoverage']['reasonCode'].astype(int)
        flags={k:np.asarray(a,dtype=bool) for k,a in m['heatingCoverage']['flags'].items()}
        solved=np.asarray(m['heating']['eligibleFace'],dtype=bool);missing=~solved
        with paths[1].open(encoding='utf-8-sig',newline='') as stream:
            rows=list(csv.DictReader(stream))
        assert n==19752==len(rows) and abs(area.sum()-1.9646292565715551)<2e-13
        assert np.array_equal(np.array([int(r['face']) for r in rows]),np.arange(n)+1)
        assert np.array_equal(code,np.array([int(r['heat_reason_code']) for r in rows]))
        assert np.max(np.abs(area-np.array([float(r['area_m2']) for r in rows])))<1e-14
        assert np.array_equal(solved,np.isfinite([float(r['heat_flux_W_m2']) for r in rows]))
        def summary(mask):
            ids=np.flatnonzero(mask)
            return dict(faces=len(ids),areaM2=float(area[ids].sum()),
                        xRangeM=[float(centers[ids,0].min()),float(centers[ids,0].max())] if len(ids) else None)
        result['primaryReasons']=[dict(code=i,**summary(code==i)) for i in np.unique(code)]
        result['flagOverlapOnMissing']={k:summary(a&missing) for k,a in flags.items()}
        edge=missing & np.isin(code,[4,14])
        result['edgePrimary']=summary(edge)
        result['edgePrimary']['fractionOfMissingArea']=float(area[edge].sum()/area[missing].sum())
        other=np.zeros(n,dtype=bool)
        for k,a in flags.items():
            if k not in ('externalFlow','upstreamExternalState'):other|=a
        result['edgeFlagWithoutOtherExistingFlags']=summary(missing & (flags['externalFlow']|flags['upstreamExternalState']) & ~other)
        # 模块二：几何连通只用于定位，不能把连通缺值块称为分离泡。
        adjacency=[set() for _ in range(n)];owners={}
        for i,tri in enumerate(f):
            for j in range(3):
                key=tuple(sorted((int(tri[j]),int(tri[(j+1)%3]))))
                if key in owners:
                    other_id=owners[key];adjacency[i].add(other_id);adjacency[other_id].add(i)
                else: owners[key]=i
        def components(mask):
            unseen=set(np.flatnonzero(mask).tolist());groups=[]
            while unseen:
                seed=min(unseen);unseen.remove(seed);stack=[seed];group=[seed]
                while stack:
                    for neighbour in adjacency[stack.pop()]:
                        if neighbour in unseen:
                            unseen.remove(neighbour);stack.append(neighbour);group.append(neighbour)
                ids=np.array(sorted(group));mask_group=np.zeros(n,dtype=bool);mask_group[ids]=True
                groups.append(dict(**summary(mask_group),firstFace=int(ids[0])+1,faceIDs=(ids+1).tolist()))
            return sorted(groups,key=lambda a:a['areaM2'],reverse=True)
        groups=components(edge);result['edgeConnectedGroups']=[{k:a[k] for k in ('faces','areaM2','xRangeM','firstFace')} for a in groups]
        # 模块三：沿既有前驱寻找热力状态拒绝根，仅作旧算法传播审计。
        predecessor=m['surfaceStreamlines']['predecessorFace'].astype(int)-1
        direct=np.asarray(m['continuumRegime']['directStateUnresolvedFace'],dtype=bool)
        root=np.full(n,-1,dtype=int);cycles=[]
        for i in np.flatnonzero(edge):
            j=int(i);seen=set()
            while j>=0 and j<n:
                if j in seen:cycles.append(int(i)+1);break
                if direct[j]:root[i]=j;break
                seen.add(j);j=int(predecessor[j])
        result['propagation']=dict(directRootsFoundFaces=int(np.sum(root[edge]>=0)),
            noDirectRootFaces=int(np.sum(root[edge]<0)),cycleCount=len(cycles),
            uniqueDirectRootFaces=int(len(np.unique(root[edge & (root>=0)]))))
        xs=np.unique(np.round(v[:,0],9));profile=[]
        for x in xs:
            ids=np.abs(v[:,0]-x)<1e-8
            profile.append(dict(xM=float(x),radiusM=float(np.hypot(v[ids,1],v[ids,2]).max())))
        result['axialGeometry']=profile
        bins=np.linspace(float(v[:,0].min()),float(v[:,0].max()),25)
        result['axialMissing']=[dict(xLowM=float(a),xHighM=float(b),
            **summary(missing&(centers[:,0]>=a)&(centers[:,0]<b)),
            edgePrimaryAreaM2=float(area[edge&(centers[:,0]>=a)&(centers[:,0]<b)].sum())) for a,b in zip(bins[:-1],bins[1:])]
        with (out/'original_face_diagnosis.csv').open('w',encoding='utf-8',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['faceID','xM','yM','zM','areaM2','primaryReason',
                'edgePrimaryCandidate','directStateRootFace','existingNonEdgeFlag'])
            for i in range(n):writer.writerow([i+1,*centers[i],area[i],int(code[i]),int(edge[i]),int(root[i]+1),int(other[i])])
        (out/'edge_components.json').write_text(json.dumps(groups,indent=2)+'\n',encoding='utf-8')
        result.update(passed=True,originalFaces=n,totalAreaM2=float(area.sum()),missing=summary(missing),
                      sourceSHA256=before,inputStable=all(sha(p)==h for p,h in before.items()),
                      independentCSVCompared=True,sourceCondition=json.loads(paths[2].read_text(encoding='utf-8'))['condition'],
                      interpretation='Existing engineering rejection/propagation; neither proven separation nor recoverable/qualified heat area.')
        assert result['inputStable']
    except Exception as exc:
        result['error']=repr(exc);raise
    finally:
        result['elapsedSeconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('passed','missing','edgePrimary','propagation','elapsedSeconds')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    analyze(parser.parse_args().out)
