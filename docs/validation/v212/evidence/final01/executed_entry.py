"""只读独立审核v212原面职责/拓扑和现成回执，归档小证据；不重跑流动。"""
from pathlib import Path
import json,hashlib,shutil,ast,time
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v212_validation';DEST=ROOT/'docs/validation/v212/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    shutil.copyfile(Path(__file__).resolve(),out/'executed_entry.py')
    design=json.loads((BASE/'design01/regions.json').read_text());sr=json.loads((BASE/'source01/report.json').read_text())
    packet=json.loads((BASE/'packet_work02/audit/report.json').read_text());proc=json.loads((BASE/'packet_process02/report.json').read_text())
    first=json.loads((BASE/'packet01/report.json').read_text());firstproc=json.loads((BASE/'packet_process01/report.json').read_text())
    assert packet['passed'] and packet['rejectedCount']==48 and proc['passed'] and proc['nativeExitCode']==0
    assert proc['inputsStable'] and not proc['liveGroupAfterCleanup'] and len(proc['inputSHA256'])==8
    assert not first['passed'] and 'FileNotFoundError' in first['error'] and firstproc['status']=='running'
    assert sr['passed'] and len(sr['files'])==46
    assert all(sha(BASE/'source01/source'/r['path'])==r['sha256'] for r in sr['files'])
    # 模块一：独立从原面连接检查连通、无重叠责任、完整分母及真实层拓扑公式。
    volume=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(volume)==design['motherVolumeSHA256']
    d=loadmat(str(volume),variable_names=['sourceVertices','sourceFaces']);v=d['sourceVertices'];f=d['sourceFaces'].astype(int)-1
    area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)/2
    owners=[];checked={}
    for name,r in design['regions'].items():
        faces=f[np.array(r['originalFaceIDs'])-1];ids=np.array(r['ownedFaceIDs'])-1;owners.extend(ids.tolist())
        assert set(r['ownedFaceIDs']).issubset(set(r['originalFaceIDs']))
        assert abs(float(area[ids].sum())-r['ownedAreaM2'])<1e-12
        by_edge={};adj=[set() for _ in faces]
        for i,tri in enumerate(faces):
            for a,b in [(0,1),(1,2),(2,0)]:
                edge=tuple(sorted((int(tri[a]),int(tri[b]))));by_edge.setdefault(edge,[]).append(i)
        for group in by_edge.values():
            assert len(group)<=2
            if len(group)==2:adj[group[0]].add(group[1]);adj[group[1]].add(group[0])
        visited={0};todo=[0]
        while todo:
            new=adj[todo.pop()]-visited;visited.update(new);todo.extend(new)
        assert len(visited)==len(faces)
        nf=len(faces);ne=len(by_edge);nb=sum(len(x)==1 for x in by_edge.values());nv=len(np.unique(faces))
        c=dict(cells=32*nf,ghosts=2*nf+32*nb,faces=33*nf+32*ne,vertices=33*nv,
            internalFaces=31*nf+32*(ne-nb),cellVertexLinks=192*nf,cellFaceLinks=160*nf,
            faceVertexLinks=99*nf+128*ne)
        assert c==r['counts']
        n=c['cells'];total=(13217*(n+c['ghosts'])+3003*c['faces']+745*c['vertices']+3170*n+
            32*c['cellVertexLinks']+124*c['cellFaceLinks']+24*c['faceVertexLinks']+584*(n+2*c['internalFaces'])+344464)
        assert total==r['layerPayloadLowerBoundBytes']
        checked[name]=dict(edgeConnected=True,counts=c,payloadBytes=total)
    owners.extend((np.array(design['remainder']['faceIDs'])-1).tolist())
    assert sorted(owners)==list(range(19752)) and abs(float(area.sum())-design['originalAreaM2'])<1e-12
    assert set(design['rootsFaceIDs']).issubset(set(design['regions']['contraction']['ownedFaceIDs']))
    # 模块二：执行时计划与成功入口只在SHA完全匹配时恢复，不凭记忆冒充执行快照。
    plan=(ROOT/'docs/validation/v212/PLAN.md').read_text(encoding='utf8')
    want=json.loads((BASE/'design01/report.json').read_text())['inputSHA256']['docs/validation/v212/PLAN.md']
    options=[plan.split('\n首次测试')[0].rstrip()+'\n'*n for n in (1,2,3)]
    exact=[s for s in options if hashlib.sha256(s.encode('utf8')).hexdigest()==want];assert len(exact)==1
    (BASE/'design01/plan_at_design.md').write_text(exact[0],encoding='utf8')
    wrapper=ROOT/'scripts/runSequentialSnapshotAuditV212.py'
    text=wrapper.read_text(encoding='utf8');candidate='\n'.join(line for line in text.split('\n') if "raise ValueError('v212唯一恢复测试" not in line)
    assert hashlib.sha256(candidate.encode('utf8')).hexdigest()==proc['inputSHA256']['/mnt/e/RocketAeroHeat/scripts/runSequentialSnapshotAuditV212.py']
    (BASE/'packet_process02/entry_at_success.py').write_text(candidate,encoding='utf8')
    assert sha(ROOT/'scripts/sequentialRegionSnapshot.py')==packet['inputSHA256']['scripts/sequentialRegionSnapshot.py']
    assert sha(ROOT/'tests/testSequentialRegionSnapshot.py')==packet['inputSHA256']['tests/testSequentialRegionSnapshot.py']
    assert sha(BASE/'packet_work02/audit/synthetic.npz')==packet['snapshotReceipt']['sha256']
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/k.replace('\\','/'))==v for k,v in protected.items())
    for name in ['inspectSequentialRegionCoupling.py','designSequentialRegions.py','runSequentialSnapshotAuditV212.py']:
        tree=ast.parse((ROOT/'scripts'/name).read_text(encoding='utf8'));assert isinstance(tree.body[-1].body[0],ast.Raise)
    save(out/'report.json',dict(passed=True,originalSurfaceExactCover=True,independentGeometry=checked,
        nativeSourcesUnchanged=46,packetRejections=48,protectedFilesUnchanged=29,
        originalPacketFailed=True,originalSupervisorAcceptance=False,firstSupervisorEndedBy='agent SIGKILL after exact command and zombie child check',
        firstSupervisorFinalMetricsKnown=False,geometryRSSMeasured=False,geometryResourceAcceptance=False,
        recoveredPacketProcess=proc,fullEngineBuilds=0,volumeMeshesGenerated=0,CFDRuns=0,
        nativeOperatorEquivalence=False,newQualifiedAreaM2=0,goalComplete=False,seconds=time.monotonic()-start))
    # 模块三：只摘录所需原生段落，原上游树仍留E盘，不上传全部源码副本。
    selections={
        'src/lmr/lua-modules/bc.lua':[(1292,1315),(1700,1731)],
        'src/lmr/flowstate.d':[(436,490),(581,677)],
        'src/lmr/ufluidblock.d':[(669,705),(878,905),(998,1020)],
        'src/lmr/bc/ghost_cell_effect/mapped_cell_copy.d':[(43,65),(160,184),(444,454),(614,624),(720,731),(928,935),(1010,1020),(1357,1365),(1550,1557),(1640,1670)],
        'src/lmr/simcore_exchange.d':[(107,135),(152,180)]}
    extra={'lmr/newtonkrylovsolver.d':[(3359,3380),(3425,3478)],'lmr/fluidblock.d':[(1324,1342),(1627,1641)],
        'lmr/lsqinterp.d':[(252,290),(316,331),(475,520)],'lmr/flowgradients.d':[(112,136)],
        'lmr/fluidfvcell.d':[(313,324)]}
    excerpts=[]
    for path,ranges in list(selections.items())+list(extra.items()):
        source=BASE/'source01/source'/path if path.startswith('src/') else ROOT/'outputs/v211_validation/source01/source'/path
        lines=source.read_text(encoding='utf8').splitlines()
        excerpts.append(dict(path=path if path.startswith('src/') else 'src/'+path,sha256=sha(source),
            ranges=[dict(start=a,end=b,text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
    save(BASE/'source_excerpts.json',excerpts)
    save(BASE/'manual_tool_notes.json',dict(
        firstSupervisor=dict(originalReportPreserved=True,observedCommand='python3 /mnt/e/RocketAeroHeat/scripts/runSequentialSnapshotAuditV212.py',
            supervisorPID=12,childPID=26,childState='Zs',observedElapsedAtCheck='01:16',termination='SIGKILL to verified supervisor only',
            finalWaitCodeUnknown=True,finalBudgetMetricsUnknown=True,taskToolExitCode=1,laterNoLiveTask=True),
        readFailures=['Windows rg glob expressions/guessed src path failed; corrected to actual source/lmr path',
            'initial search term InFlowBC_FromProfile differs from real StaticProfile; no absence claim'],
        plot=dict(initialAxisOverlapPreserved=True,correctedWithoutNumericChange=True,renderSeconds=[.516,.437,.453]),
        sourceWeb=[dict(url='https://arxiv.org/html/2410.04668v1',read='sections 2.2 and 2.3; no target-SA proof'),
            dict(url='https://ntrs.nasa.gov/citations/19940028614',read='abstract only'),
            dict(url='https://www.sciencedirect.com/science/article/pii/S0045793020300256',read='403; not used as full text')]))
    DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        if 'source01/source/' in p.relative_to(BASE).as_posix():
            local.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size));continue
        assert p.stat().st_size<2*1024**2
        q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(p)==sha(q);rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    names=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_SEQUENTIAL_REGIONS_V212.md',
        'docs/validation/v212/PLAN.md','docs/validation/v212/NEXT_DECISION.md',
        'scripts/inspectSequentialRegionCoupling.py','scripts/designSequentialRegions.py','scripts/sequentialRegionSnapshot.py',
        'scripts/runSequentialSnapshotAuditV212.py','scripts/plotSequentialRegionsV212.py',
        'scripts/finalizeSequentialRegionsV212.py','tests/testSequentialRegionSnapshot.py']
    for name in names:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name)
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,sha256=sha(p),bytes=p.stat().st_size))
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    size=sum(r['bytes'] for r in rows);assert size<6*1024**2
    save(ROOT/'docs/validation/v212/evidence_manifest.json',dict(schema='RAH.V212Evidence.v1',files=rows,localSources=local,
        count=len(rows),bytes=size,sourceFiles=names,goalComplete=False,protectedFilesUnchanged=29,uniqueDeleted=0,
        CFDRuns=0,fullEngineBuilds=0,newQualifiedAreaM2=0,seconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True,count=len(rows),bytes=size,seconds=time.monotonic()-start)))

if __name__=='__main__':main()
