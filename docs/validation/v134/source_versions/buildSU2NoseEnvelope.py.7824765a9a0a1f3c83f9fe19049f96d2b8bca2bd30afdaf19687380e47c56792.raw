"""保留原壁面及近壁列，用有界凸性约束修整鼻部外壳；独立保存几何准备结果。"""
from pathlib import Path
import argparse,json,os,subprocess,sys,time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,vstack,eye,hstack
from meshSU2HybridVolume import wedge_jacobian_min
from prepareSU2CheckpointContinuation import ROOT,read,sha

def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')

def edge_stencils(F):
    """相邻三角形的四点模板；原闭合面每条边恰有两邻面。"""
    E=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]])
    opposite=np.concatenate([F[:,2],F[:,0],F[:,1]])
    keys=np.sort(E,axis=1);order=np.lexsort(keys.T[::-1])
    pair=E[order].reshape(-1,2,2);other=opposite[order].reshape(-1,2)
    assert np.array_equal(np.sort(pair[:,0],axis=1),np.sort(pair[:,1],axis=1))
    return np.column_stack([pair[:,0],other])

def bounded_envelope(V,F,S,H):
    """仅移动外层：拟合参考球面后，用HiGHS保证保留连接的局部凸性。"""
    moving=V[:,0]>.03;ids=np.flatnonzero(moving);nv=len(V)
    A=np.column_stack([2*V[moving],np.ones(sum(moving))]);b=np.sum(V[moving]**2,axis=1)
    fit=np.linalg.lstsq(A,b,rcond=None)[0];center=fit[:3];radius=np.sqrt(fit[3]+center@center)
    residual=np.linalg.norm(V[moving]-center,axis=1)-radius
    u=np.clip((V[:,0]-.03)/.01,0,1);blend=u*u*(3-2*u)
    ideal=center+(radius+H)*(V-center)/np.linalg.norm(V-center,axis=1)[:,None]
    target=S+blend[:,None]*(ideal-S)
    st=edge_stencils(F);st=st[np.any(moving[st],axis=1)]
    # 在YZ投影平面求对侧点相对第一个三角形的仿射坐标，X约束保持线性。
    y=target[st,1:];mat=np.stack([y[:,1]-y[:,0],y[:,2]-y[:,0]],axis=2)
    det=np.linalg.det(mat);assert np.all(det>0),'鼻部三角形YZ投影方向非法'
    ab=np.linalg.solve(mat,y[:,3]-y[:,0]);weights=np.column_stack([1-ab.sum(axis=1),ab])
    coef=np.column_stack([-weights,np.ones(len(st))])
    # 正安全量为每条边长的百万分之一，不按结果调节。
    margin=1e-6*np.linalg.norm(target[st[:,1]]-target[st[:,0]],axis=1)
    reference=np.sum(coef*target[st,0],axis=1)
    index=np.full(nv,-1,int);index[ids]=np.arange(len(ids))
    rr=np.repeat(np.arange(len(st)),4);cc=index[st.ravel()];valid=cc>=0
    C=coo_matrix((coef.ravel()[valid],(rr[valid],cc[valid])),shape=(len(st),len(ids))).tocsr()
    I=eye(len(ids),format='csr');ones=-np.ones((len(ids),1))
    Aub=vstack([hstack([C,np.zeros((len(st),1))]),hstack([I,ones]),hstack([-I,ones])],format='csr')
    bub=np.r_[(-margin-reference)/.001,np.zeros(2*len(ids))]
    dyz=np.sum((target[ids,1:]-S[ids,1:])**2,axis=1)
    assert np.all(dyz<.001**2),'单节点YZ位移已超1mm'
    allowance=np.sqrt(.001**2-dyz)
    bounds=list(zip((S[ids,0]-target[ids,0]-allowance)/.001,(S[ids,0]-target[ids,0]+allowance)/.001))+[(0,None)]
    cost=np.r_[np.zeros(len(ids)),1.]
    result=linprog(cost,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs',options={'time_limit':60.})
    report=dict(method='bounded-spherical-reference-convex-lifting-highs',movingNodes=len(ids),constraints=len(st),
        sphereCenterM=center.tolist(),sphereRadiusM=float(radius),fitRadialRMSM=float(np.sqrt(np.mean(residual**2))),
        fitMaximumRadialErrorM=float(max(abs(residual))),minimumYZProjectedDeterminant=float(min(det)),
        displacementBoundM=.001,minimumMarginEdgeRatio=1e-6,optimizerStatus=int(result.status),
        optimizerMessage=result.message,optimizerSuccess=bool(result.success),wallGeometryChanged=False)
    if not result.success:return None,report
    candidate=target.copy();candidate[ids,0]+=.001*result.x[:-1]
    violation=np.sum(coef*candidate[st,0],axis=1)+margin
    displacement=np.linalg.norm(candidate-S,axis=1)
    report.update(maximumConstraintViolationM=float(max(violation)),maximumDisplacementM=float(max(displacement)),
        minimaxXDeviationFromReferenceM=float(.001*result.fun),
        maximumOriginalConcavityM=float(max(np.sum(coef*S[st,0],axis=1))),
        maximumCandidateConcavityM=float(max(violation-margin)))
    assert max(violation)<=1e-11 and max(displacement)<=.001+1e-12
    assert np.array_equal(candidate[~moving],S[~moving])
    return candidate,report

def build(source,folder):
    process=read(source/'mesh_process.json');assert process['exitCode']==0 and process['inputStable'] and process['sourceStable']
    for n in ['hybrid_volume.mat','mesh_report.json']:assert sha(source/n)==process['outputSHA256'][n]
    old=loadmat(str(source/'hybrid_volume.mat'));r=read(source/'mesh_report.json')
    V=old['sourceVertices'];F=old['sourceFaces'].astype(int)-1;P=old['vertices'].copy();nv=len(V)
    h=np.array(r['normalHeightsM']);nl=len(h)-1;S=P[nl*nv:(nl+1)*nv].copy()
    candidate,detail=bounded_envelope(V,F,S,h[-1]);write(folder/'envelope_report.json',detail)
    assert candidate is not None,'凸性约束不可行，不扫描参数'
    delta=candidate-S
    for j in range(nl+1):
        u=max(0.,(h[j]-.001)/(h[-1]-.001));weight=u*u*(3-2*u)
        if weight>0:P[j*nv:(j+1)*nv]+=weight*delta
    jac=wedge_jacobian_min(P[old['wedges'].astype(int)-1])
    detail.update(minimumPrismJacobian=float(min(jac)),allPrismJacobiansPositive=bool(np.all(jac>0)),
        protectedHeightsM=h[h<=.001].tolist(),protectedLayerCount=int(sum(h<=.001)-1))
    write(folder/'envelope_report.json',detail)
    assert min(jac)>0,'过渡棱柱翻转，不裁剪/删除单元'
    assert np.array_equal(P[:sum(h<=.001)*nv],old['vertices'][:sum(h<=.001)*nv])
    # 只保存新的棱柱域和原远场；旧外四面体不再有效，不导出成可用完整网格。
    savemat(str(folder/'collar_prepared.mat'),dict(vertices=P,wedges=old['wedges'],farfield=old['farfield'],
        sourceVertices=V,sourceFaces=F+1,oldOuterShell=S,newOuterShell=candidate),do_compression=True)
    r.update(noseEnvelope=detail,noseEnvelopePrepared=True,exteriorRequiresRebuild=True,
        physicsQualified=False,newQualifiedArea=0)
    write(folder/'collar_report.json',r)
    savemat(str(folder/'mesh_input.mat'),loadmat(str(source/'mesh_input.mat')),do_compression=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--target',type=Path,required=True)
    p.add_argument('--worker',action='store_true');a=p.parse_args();source=a.source.resolve();folder=a.target.resolve()
    source.relative_to(ROOT/'outputs');folder.relative_to(ROOT/'outputs')
    if a.worker:build(source,folder);return
    assert not folder.exists();folder.mkdir(parents=True)
    files=[source/n for n in ['mesh_process.json','mesh_report.json','hybrid_volume.mat','mesh_input.mat']]+[
        Path(__file__).resolve(),ROOT/'scripts/meshSU2HybridVolume.py',ROOT/'docs/validation/v134/PLAN.md']
    locks={str(p):sha(p) for p in files};start=time.monotonic();receipt=folder/'prepare_process.json'
    record=dict(status='running',maximumSeconds=180,sourceHashes=locks,sourceMeshDirectory=str(source),
        gmshRun=False,nativeFlowRuns=0,supervisorPID=os.getpid())
    write(receipt,record)
    with (folder/'prepare.log').open('wb') as log:
        proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker','--source',str(source),'--target',str(folder)],
            stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,
            env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
        record['processID']=proc.pid;write(receipt,record);timed=False
        try:code=proc.wait(timeout=180)
        except subprocess.TimeoutExpired:timed=True;proc.kill();code=proc.wait()
    stable={p:sha(Path(p))==value for p,value in locks.items()}
    record.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start,
        sourceStableAfter=stable,sourceStable=all(stable.values()),inputStable=all(stable.values()),
        outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p!=receipt})
    write(receipt,record);print(json.dumps({k:record[k] for k in ['status','exitCode','elapsedSeconds']}))
    raise SystemExit(124 if timed else code)

if __name__=='__main__':main()
