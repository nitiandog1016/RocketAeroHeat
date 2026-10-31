"""针对已定位ω下限停滞口袋的显式初场研究；保留原生模型与所有其他场量。"""
from pathlib import Path
import argparse,copy,json,os
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import OUTPUT_CONTROL,validate_source_artifacts
from prepareSU2SSTCase import set_value
from su2SSTInterface import state,mesh

KIND='same-mesh-SST-omega-pocket-initialization-study'
STUDY='v128-base-omega-floor-logharmonic-initial-guess'


def floor_nodes(omega):
    """15位科学记数及除密度的舍入界；不把数值下限提高为新的物理边界。"""
    return np.flatnonzero(omega<=1e-4*(1+8*np.finfo(float).eps))


def reconstruct(points,cells,wallIDs,omega):
    """只解下限节点的log场图Dirichlet问题；边界固定，不算物理求解。"""
    ids=floor_nodes(omega);assert len(ids)>0 and not set(ids)&set(wallIDs)
    lookup={int(i):j for j,i in enumerate(ids)};neighbors={int(i):set() for i in ids}
    edges={10:[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
        13:[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]}
    for row in cells:
        typ=row[0];v=row[1:];assert typ in edges
        if not any(i in lookup for i in v):continue
        for ia,ib in edges[typ]:
            i,j=v[ia],v[ib]
            if i in lookup:neighbors[i].add(j)
            if j in lookup:neighbors[j].add(i)
    rr=[];cc=[];vv=[];b=np.zeros(len(ids));anchors=np.zeros(len(ids),bool);boundary=[]
    for i,js in neighbors.items():
        ii=lookup[i];weightSum=0
        for j in sorted(js):
            w=1/np.linalg.norm(points[i]-points[j]);assert np.isfinite(w) and w>0;weightSum+=w
            if j in lookup:rr.append(ii);cc.append(lookup[j]);vv.append(-w)
            else:b[ii]+=w*np.log(omega[j]);anchors[ii]=True;boundary.append(omega[j])
        rr.append(ii);cc.append(ii);vv.append(weightSum)
    A=csr_matrix((vv,(rr,cc)),shape=(len(ids),len(ids)))
    nc,labels=connected_components(A,directed=False)
    assert all(any(anchors[labels==k]) for k in range(nc))
    logw=spsolve(A,b);values=np.exp(logw)
    err=float(np.max(abs(A@logw-b))/max(np.max(abs(b)),1));assert err<1e-12
    assert min(values)>=min(boundary)*(1-1e-12) and max(values)<=max(boundary)*(1+1e-12)
    assert np.all(values>1e-4)
    return ids,values,dict(nodes=len(ids),components=int(nc),linearEquationRelativeError=err,
        boundaryRange=[float(min(boundary)),float(max(boundary))],newOmegaRange=[float(min(values)),float(max(values))])


def validate_initial(folder,source,init):
    """逐字段证明只改已列出的ω节点；未改字段保持解析值精确相等。"""
    a,keys,nd,_,_,_=state(source/'restart_flow.csv',True,True)
    b,newkeys,_,_,_,_=state(folder/'initial_flow.csv',True,True);assert newkeys==keys and nd==3
    col=keys.index('Omega');indices=np.flatnonzero(a[:,col]!=b[:,col]);rep=read(folder/'pocket_initialization.json')
    assert sha(folder/'pocket_initialization.json')==init['pocketReportSHA256']
    assert indices.tolist()==rep['nodeIDs'] and len(indices)==146
    assert np.array_equal(np.delete(a,col,axis=1),np.delete(b,col,axis=1))
    assert np.array_equal(b[indices,col],np.array(rep['newOmega']))
    assert np.array_equal(indices,floor_nodes(a[:,col]))
    assert np.all(b[indices,col]>1e-4)
    assert rep['sourceSHA256']==sha(source/'restart_flow.csv') and rep['meshSHA256']==sha(folder/'mesh.su2')
    assert rep['initialGuessOnly'] and not rep['physicsQualified']


def validate_source(folder,c,init):
    source=Path(init['sourceDirectory']).resolve();a=config(source/'case.cfg')
    assert init['kind']==KIND and init['studyID']==STUDY and init['variant'] in ['control','repaired']
    assert init['sourceLocalIteration']==149 and init['localIterationStartsAt']==0
    assert init['sameMesh'] and init['samePhysicsAndNumerics'] and not init['physicsQualified']
    assert init['changedField']==('Omega-at-146-interior-floor-nodes' if init['variant']=='repaired' else 'none')
    assert {k:v for k,v in c.items() if k not in OUTPUT_CONTROL}=={k:v for k,v in a.items() if k not in OUTPUT_CONTROL}
    assert c['ITER']=='10' and c['OUTPUT_WRT_FREQ']=='10'
    assert c['VOLUME_OUTPUT']=='(COORDINATES, SOLUTION, PRIMITIVE)'
    locked=validate_source_artifacts(folder,init,validate_initial if init['variant']=='repaired' else None)
    planName=init.get('planRelativePath','docs/validation/v128/POCKET_PLAN.md')
    assert planName in ['docs/validation/v128/POCKET_PLAN.md','docs/validation/v128/POCKET_ROUNDOFF.md']
    plan=ROOT/planName;assert sha(plan)==init['planSHA256']
    locked.update({str(plan):sha(plan),str(Path(__file__).resolve()):sha(Path(__file__).resolve())})
    if init['variant']=='repaired':locked[str(folder/'pocket_initialization.json')]=init['pocketReportSHA256']
    return locked


def prepare(variant):
    source=ROOT/'outputs/v127_validation/control';folder=ROOT/'outputs/v128_validation'/('pocket_control' if variant=='control' else 'pocket_repaired_roundoff')
    assert not folder.exists();folder.mkdir(parents=True)
    os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'));text=(source/'case.cfg').read_text()
    for key,value in dict(ITER=10,OUTPUT_WRT_FREQ=10,SCREEN_WRT_FREQ_INNER=1).items():text=set_value(text,key,value)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    if variant=='control':os.link(str(source/'restart_flow.csv'),str(folder/'initial_flow.csv'))
    else:
        a,keys,_,_,_,_=state(source/'restart_flow.csv',True,True);p,cells,wall=mesh(source/'mesh.su2')
        ids,omega,report=reconstruct(p,cells,np.unique([i for r in wall for i in r[1:]]),a[:,keys.index('Omega')])
        assert len(ids)==146;mapping=dict(zip(ids,omega));col=keys.index('Omega')
        # 每个未改行原字节保留；改行只替换omega一项，其他字段不重新格式化。
        with (source/'restart_flow.csv').open('rb') as old,(folder/'initial_flow.csv').open('wb') as new:
            new.write(next(old))
            for i,line in enumerate(old):
                if i in mapping:
                    values=line.split(b',');values[col]=(' %.17e'%mapping[i]).encode('ascii');line=b','.join(values)
                new.write(line)
        report.update(nodeIDs=ids.tolist(),newOmega=omega.tolist(),exactFloorNodeIDs=np.flatnonzero(a[:,col]==1e-4).tolist(),
            sourceSHA256=sha(source/'restart_flow.csv'),
            meshSHA256=sha(folder/'mesh.su2'),initialGuessOnly=True,physicsQualified=False)
        (folder/'pocket_initialization.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    locks=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json','case_manifest.json',
        'restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    planRelative='docs/validation/v128/POCKET_PLAN.md' if variant=='control' else 'docs/validation/v128/POCKET_ROUNDOFF.md'
    init=dict(kind=KIND,studyID=STUDY,variant=variant,changedField='Omega-at-146-interior-floor-nodes' if variant=='repaired' else 'none',
        sourceDirectory=str(source),sourceLocalIteration=149,sourceWasDiagnostic=False,localIterationStartsAt=0,
        sourceStateHasNativePrimitive=True,sameMesh=True,samePhysicsAndNumerics=True,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False,diagnosticOnly=True,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(source/k):sha(source/k) for k in locks},
        planRelativePath=planRelative,planSHA256=sha(ROOT/planRelative))
    if variant=='repaired':init['pocketReportSHA256']=sha(folder/'pocket_initialization.json')
    validate_source(folder,config(folder/'case.cfg'),init)
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n',encoding='utf-8')
    contract=copy.deepcopy(read(source/'sst_contract.json'));contract.pop('numericalStudy',None)
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        initializationStudy={k:init[k] for k in ['studyID','variant','changedField']})
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    m=copy.deepcopy(read(source/'case_manifest.json'));m['flight'].update(iterations=10,outputEvery=10,namedCheckpoints=True)
    m['sourceRunControlsFromConfig']={k:config(source/'case.cfg')[k] for k in ['ITER','OUTPUT_WRT_FREQ']}
    m.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    req=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=180,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTPocketStudy.py','docs/validation/v128/POCKET_PLAN.md',planRelative])
    (folder/'engine_request.json').write_text(json.dumps(req,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(variant=variant,folder=str(folder),newLocalLastIteration=9)))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('variant',choices=['control','repaired']);prepare(parser.parse_args().variant)
