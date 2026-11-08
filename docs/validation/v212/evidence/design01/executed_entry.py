"""按原面和原层列设计顺序重叠区域；只做组合选择，不生成新体网格。"""
from pathlib import Path
import json, time, shutil
import numpy as np
from scipy.io import loadmat
from eilmerNKSAResourceBudget import ROOT, VOLUME_SHA, EXPECTED_SIZES, lower_bound, extruded_counts, sha, save

BASE = ROOT / 'outputs/v212_validation'

def design(out):
    """模块一：按预定物理位置取最近原环，原鼻部/收缩根面/原扩域全部保留。"""
    start = time.monotonic(); out = Path(out); out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(), out / 'executed_entry.py')
    report = dict(passed=False, CFDRuns=0, volumeMeshesGenerated=0, runAuthorized=False)
    try:
        source = ROOT / 'outputs/v144_validation/actual/hybrid_volume.mat'
        assert sha(source) == VOLUME_SHA
        d = loadmat(str(source), variable_names=['sourceVertices', 'sourceFaces', 'wedges'])
        v = d['sourceVertices']; f = d['sourceFaces'].astype(np.int64)-1; w = d['wedges'].astype(np.int64)-1
        area = np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]], v[f[:,2]]-v[f[:,0]]), axis=1)/2
        xmin = v[f,0].min(axis=1); xmax = v[f,0].max(axis=1); rings = np.unique(v[:,0])
        ring = lambda x: float(rings[np.argmin(abs(rings-x))])
        old = json.loads((ROOT/'docs/validation/v189/region_contract.json').read_text())['local']
        roots = json.loads((ROOT/'docs/validation/v189/evidence/selection01/report.json').read_text())['thermalRoots']['faceIDs']
        # 几何设计点不是流动附着/超声速判断；重叠宽度仍须实际误差验证。
        cut_a, start_b, cut_b, start_c = map(ring, [-.40, -.28, -.88, -.70])
        owner_ab, owner_bc = map(ring, [-.34, -.79])
        end = old['extensionTestCutXM']; tol = 1e-10
        cover = xmin >= end-tol
        solves = dict(nose=xmin>=cut_a-tol,
                      contraction=(xmin>=cut_b-tol)&(xmax<=start_b+tol),
                      afterbody=(xmin>=end-tol)&(xmax<=start_c+tol))
        owners = dict(nose=xmin>=owner_ab-tol,
                      contraction=(xmin>=owner_bc-tol)&(xmax<=owner_ab+tol),
                      afterbody=cover&(xmax<=owner_bc+tol), remainder=~cover)
        assert np.all(sum(owners.values())==1) and np.all(np.logical_or.reduce(list(solves.values()))==cover)
        assert np.all(owners['contraction'][np.array(roots)-1]) and len(roots)==336
        assert abs(float(area.sum())-1.9646292565715553)<1e-12 and len(f)==19752
        regions = {}
        # 模块二：读取已有32层柱，不改坐标/首层/连接，必要载荷与独立整数公式交叉核对。
        for name,mask in solves.items():
            assert np.all(~owners[name] | mask)
            counts, topo = extruded_counts(f[mask], w, len(v))
            bound = lower_bound(counts, EXPECTED_SIZES)
            n=counts['cells']; c=n+counts['ghosts']
            independent=(13217*c+3003*counts['faces']+745*counts['vertices']+3170*n+
                32*counts['cellVertexLinks']+124*counts['cellFaceLinks']+24*counts['faceVertexLinks']+
                584*(n+2*counts['internalFaces'])+344464)
            assert independent==bound['payloadLowerBoundBytes']
            # 每个子域的表面都必须是完整的圆周片；检查边界环、连通和欧拉数。
            edges=np.sort(f[mask][:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
            unique, freq=np.unique(edges,axis=0,return_counts=True)
            boundary=unique[freq==1]; bv,degree=np.unique(boundary,return_counts=True)
            assert np.all(degree==2) and np.all((freq==1)|(freq==2))
            bx=np.unique(v[bv,0]); assert len(bx)==(1 if name=='nose' else 2)
            assert len(np.unique(f[mask]))-len(unique)+int(mask.sum())==(1 if name=='nose' else 0)
            assert bound['payloadLowerBoundBytes'] < 8*1024**3
            regions[name]=dict(originalFaceIDs=(np.flatnonzero(mask)+1).tolist(),
                ownedFaceIDs=(np.flatnonzero(owners[name])+1).tolist(), faceCount=int(mask.sum()),
                ownedFaceCount=int(owners[name].sum()), areaM2=float(area[mask].sum()),
                ownedAreaM2=float(area[owners[name]].sum()), boundaryRingXM=bx.tolist(),
                topology=topo, counts=counts, layerPayloadLowerBoundBytes=independent,
                layerPayloadLowerBoundGiB=independent/1024**3,
                peakBytes=None, fitProven=False, runAuthorized=False, state='design-only-unqualified')
        overlap={}
        for a,b in [('nose','contraction'),('contraction','afterbody')]:
            ids=np.flatnonzero(solves[a]&solves[b]);assert len(ids)>0
            overlap[a+'--'+b]=dict(originalFaceIDs=(ids+1).tolist(), faceCount=len(ids),
                areaM2=float(area[ids].sum()), xminM=float(xmin[ids].min()), xmaxM=float(xmax[ids].max()))
        # 全箭责任台账：未解决的其余面明确保留，绝不把选择范围当全面积结果。
        remainder=np.flatnonzero(owners['remainder'])
        result=dict(schema='RAH.SequentialRegionDesign.v1',status='offline-design-not-native-admitted',
            motherVolumeSHA256=VOLUME_SHA, originalFaceCount=len(f), originalAreaM2=float(area.sum()),
            originalMinimumHistoricalLayers=19, selectedLayers=32, selectedFirstHeightM=2e-6,
            ownershipCutsXM=[owner_ab,owner_bc], baselineEndXM=old['cutXM'], extensionEndXM=end,
            regions=regions, overlap=overlap, rootsFaceIDs=roots,
            remainder=dict(faceIDs=(remainder+1).tolist(), faceCount=len(remainder),
                areaM2=float(area[remainder].sum()), responsibility='existing-engineering-where-applicable; unresolved-aft/base-needs-local-solution',
                newQualifiedAreaM2=0),
            newQualifiedAreaM2=0, physicalQualification=False, runAuthorized=False,
            pending=['native full residual and frozen-donor derivative equivalence',
                'complete exterior mesh and memory peak, including snapshots and all stages',
                'physical terminal pressure/backflow boundary or coupling to downstream full domain',
                'interface/overlap/extension/global residual and heat verification'])
        save(out/'regions.json', result)
        report.update(passed=True, originalFaces=len(f), coveredFaces=int(cover.sum()), remainderFaces=len(remainder),
            ownedAreasM2={k:float(area[m].sum()) for k,m in owners.items()},
            layerLowerBoundGiB={k:r['layerPayloadLowerBoundGiB'] for k,r in regions.items()},
            overlapFaces={k:r['faceCount'] for k,r in overlap.items()}, rootsInCore=336,
            inputSHA256={str(source.relative_to(ROOT)):VOLUME_SHA,
                'docs/validation/v212/PLAN.md':sha(ROOT/'docs/validation/v212/PLAN.md'),
                'scripts/designSequentialRegions.py':sha(Path(__file__).resolve())},
            regionsSHA256=sha(out/'regions.json'),newQualifiedAreaM2=0)
        assert sha(source)==VOLUME_SHA
    except Exception as exc:
        report['error']=repr(exc);raise
    finally:
        report['seconds']=time.monotonic()-start;save(out/'report.json',report)
    assert report['seconds']<60
    print(json.dumps(report,ensure_ascii=True))

if __name__=='__main__':
    design(BASE/'design01')
