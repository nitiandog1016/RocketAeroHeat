"""只读恢复资源报告并独立审核整数账；不重复原生布局或生成体网格。"""
from pathlib import Path
import copy,json,time,shutil,ast
import numpy as np
from scipy.io import loadmat
from eilmerNKSAResourceBudget import ROOT,BASE,PROFILE,MESH_SHA,VOLUME_SHA,EXPECTED_SIZES,sha,save,assess

def owner_sum(r):
    """模块一：以所有者独立重组，不调用生产逐行求和函数。"""
    d=r['counts'];n=d['cells'];c=n+d['ghosts'] if r['includeGhosts'] else n
    return (13217*c+3003*d['faces']+745*d['vertices']+3170*n+32*d['cellVertexLinks']+
        124*d['cellFaceLinks']+24*d['faceVertexLinks']+584*(n+2*d['internalFaces'])+344464)

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    report=dict(passed=False,engineBuilds=0,CFDRuns=0,newQualifiedHeatAreaM2=0,physicalQualification=False)
    try:
        original=BASE/'resource01';old=json.loads((original/'report.json').read_text())
        assert not old['passed'] and 'does not start with' in old['error']
        names=['baseLayer','baseWithoutGhosts','extensionLayer','extensionWithoutGhosts']
        rows={n:json.loads((original/(n+'.json')).read_text()) for n in names}
        hashes={n:sha(original/(n+'.json')) for n in names}
        for r in rows.values():
            assert r['payloadLowerBoundBytes']==owner_sum(r)==sum(x['count']*x['unitBytes'] for x in r['rows'])
            assert all(x['bytes']==x['count']*x['unitBytes'] and x['bytes']>=0 for x in r['rows'])
            assert not r['runAuthorized'] and not r['fitProven'] and r['fullPeakBytes'] is None
        layout=json.loads((BASE/'layout01/report.json').read_text());proc=layout['process']
        assert layout['passed'] and layout['sizes']==EXPECTED_SIZES
        assert proc['passed'] and proc['inputsStable'] and not proc['liveGroupAfterCleanup']
        sr=json.loads((BASE/'source01/report.json').read_text())
        assert all(sha(BASE/'source01/source'/p['relativePath'])==p['sha256'] for p in sr['files'])
        addon=json.loads((BASE/'source02/report.json').read_text())
        assert sha(BASE/'source02/jacobian.d')==addon['sha256']
        volume=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(volume)==VOLUME_SHA
        mesh=ROOT/'outputs/v208_validation/layer_mesh02/exact_mesh.npz';assert sha(mesh)==MESH_SHA
        # 模块二：独立按原表面连边和欧拉关系核对层域计数，不重跑原体单元选择。
        d=loadmat(str(volume),variable_names=['sourceVertices','sourceFaces']);v=d['sourceVertices'];f=d['sourceFaces'].astype(np.int64)-1
        contract=json.loads((ROOT/'docs/validation/v189/region_contract.json').read_text())['local']
        geometry=[]
        for label,cut,want_area in [('baseLayer',contract['cutXM'],contract['wallAreaM2']),
                ('extensionLayer',contract['extensionTestCutXM'],contract['extensionWallAreaM2'])]:
            mask=np.min(v[f,0],axis=1)>=cut-1e-10;faces=f[mask];count=len(faces)
            edges=np.sort(faces[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
            _,freq=np.unique(edges,axis=0,return_counts=True);assert np.all((freq==1)|(freq==2))
            assert len(np.unique(np.sort(faces,axis=1),axis=0))==count
            nv=len(np.unique(faces));nb=int(np.sum(freq==1));ne=len(freq);assert nv-ne+count==1
            expected=dict(cells=32*count,ghosts=2*count+32*nb,faces=33*count+32*ne,
                vertices=33*nv,internalFaces=31*count+32*(ne-nb),cellVertexLinks=192*count,
                cellFaceLinks=160*count,faceVertexLinks=99*count+128*ne)
            assert expected==rows[label]['counts']
            area=float(np.linalg.norm(np.cross(v[faces[:,1]]-v[faces[:,0]],v[faces[:,2]]-v[faces[:,0]]),axis=1).sum()/2)
            assert abs(area-want_area)<1e-12
            geometry.append(dict(name=label,wallFaces=count,areaM2=area,boundaryEdges=nb,cutXM=cut,
                originalFaceIDs=(np.flatnonzero(mask)+1).tolist(),counts=expected))
        save(out/'geometry_counts.json',geometry)
        # 模块三：身份/计数/布局/预算拒绝，不因预算加大而开放运行。
        args=[PROFILE,MESH_SHA,rows['baseLayer']['counts'],layout['sizes'],8192]
        assert assess(*args)['status']=='not-excluded-unproven'
        assert assess(*args[:-1],4096)['status']=='impossible'
        assert not assess(*args[:-1],16384)['runAuthorized']
        bad=[]
        for index,value in [(0,'other'),(1,'0'*64),(4,0),(4,-1),(4,8192.0),(4,True),(4,float('nan'))]:
            a=copy.deepcopy(args);a[index]=value;bad.append(a)
        for k in args[2]:
            a=copy.deepcopy(args);a[2][k]+=1;bad.append(a)
        for k in args[3]:
            a=copy.deepcopy(args);a[3][k]+=1;bad.append(a)
        for a in bad:
            try:assess(*a)
            except ValueError:continue
            raise AssertionError('错误资源身份未拒绝')
        assert len(bad)==30
        protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
        assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
        assert all(sha(original/(n+'.json'))==h for n,h in hashes.items())
        for name in ['inspectEilmerNKSAResources.py','eilmerNKSAResourceBudget.py']:
            tree=ast.parse((ROOT/'scripts'/name).read_text(encoding='utf-8'))
            assert isinstance(tree.body[-1].body[0],ast.Raise)
        bound=rows['extensionWithoutGhosts']['payloadLowerBoundBytes'];budget=8*1024**3
        assert bound>budget
        report.update(passed=True,profile=PROFILE,sizes=layout['sizes'],originalReportPassed=False,
            originalFailure=old['error'],rowFilesUnchanged=True,rowSHA256=hashes,
            lowerBoundsBytes={n:r['payloadLowerBoundBytes'] for n,r in rows.items()},
            lowerBoundsGiB={n:r['payloadLowerBoundGiB'] for n,r in rows.items()},
            baseLayerDecision=assess(*args),extensionWithoutGhostsExcessBytes=bound-budget,
            scope='Original base plus prescribed extension, complex NK / SA / ILU0 / GMRES50; not all Eilmer configurations.',
            completeOriginalResourcePlan='closed',runAuthorized=False,peakMeasured=False,peakBytes=None,
            nativeLayoutSeconds=proc['elapsedSeconds'],nativeLayoutPeakRSSBytes=proc['peakObservedGroupRSSBytes'],
            nativeLayoutInputs=len(proc['inputSHA256']),nativeLayoutGroupEmpty=True,
            sourceSeconds=sr['seconds']+addon['seconds'],sourceFiles=len(sr['files'])+1,
            rejectionTests=len(bad),budgetDecisionChecks=3,protectedFilesUnchanged=29,
            volumeMeshesGenerated=0,layerSelectionReruns=0,fullEngineBuilds=0,goalComplete=False,
            inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [volume,mesh,BASE/'layout01/report.json',
                original/'executed_entry.py',ROOT/'docs/validation/v211/READ_ONLY_RECOVERY.md',Path(__file__).resolve()]})
    except Exception as exc:
        report['error']=repr(exc);raise
    finally:
        report['seconds']=time.monotonic()-start;save(out/'report.json',report)
    assert report['seconds']<60
    print(json.dumps({k:v for k,v in report.items() if k not in ['inputSHA256','sizes','rowSHA256']}))

if __name__=='__main__':main()
