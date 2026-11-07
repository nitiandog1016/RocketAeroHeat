"""四个内部自由点的一次成熟SLSQP优化：保全部原面与连接，不改CFD内核。"""
from pathlib import Path
import json,time,shutil
import numpy as np
import scipy
from scipy.optimize import minimize
from prepareOpenFoamControl import ROOT,sha,save,header
from reviewSU2ExplicitQuality import quality


def main():
    start=time.monotonic();base=ROOT/'outputs/v200_validation';source=base/'rocket02';out=base/'rocket_relocation01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');data=np.load(source/'exact_mesh.npz');g={k:data[k] for k in data.files}
    review=json.loads((base/'rocket_geometry01/report.json').read_text());assert review['passed'] and len(review['badCells'])==19
    free=np.array([16695,16709,16721,17081]);p=g['points'];c=g['cells'];ni=len(g['neighbor']);boundary=set(tuple(sorted(f)) for f in g['faces'][ni:])
    assert not set(free)&set(g['faces'][ni:].reshape(-1));affected=np.flatnonzero(np.isin(c,free).any(axis=1));assert len(affected)==145
    local=c[affected];xyz=p[local];basepoints=p[free].copy();orders=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    mask=np.array([[tuple(sorted(row[o])) not in boundary for o in orders] for row in local])
    scales=np.array([min(np.linalg.norm(p[n]-p[v]) for n in np.unique(c[np.any(c==v,axis=1)]) if n!=v) for v in free])
    det0=np.linalg.det((xyz[:,1:]-xyz[:,:1]).transpose(0,2,1));assert np.all(det0>0)
    def metrics(x):
        q=xyz.copy();new=basepoints+x.reshape(4,3)*scales[:,None]
        for i,v in enumerate(free):q[local==v]=new[i]
        faces=q[:,orders];sf=np.cross(faces[:,:,1]-faces[:,:,0],faces[:,:,2]-faces[:,:,0])/2
        avg=np.sum(np.linalg.norm(sf,axis=2)*mask,axis=1)/mask.sum(axis=1)
        tensor=np.sum(sf[:,:,:,None]*sf[:,:,None,:]*mask[:,:,None,None],axis=1)/avg[:,None,None]**2
        det=abs(np.linalg.det(tensor));ratio=np.linalg.det((q[:,1:]-q[:,:1]).transpose(0,2,1))/det0
        return det,quality(q),ratio,new
    opt_start=time.monotonic();evaluations=0
    def constraints(x):
        nonlocal evaluations
        evaluations+=1
        if time.monotonic()-opt_start>20:raise TimeoutError('唯一局部优化20秒预算')
        d,q,v,_=metrics(x);return np.r_[d/.00105-1,q/.1-1,v/.01-1]
    failure=None;result=None
    try:
        result=minimize(lambda x:float(x@x),np.zeros(12),jac=lambda x:2*x,method='SLSQP',bounds=[(-1,1)]*12,
            constraints={'type':'ineq','fun':constraints},options={'maxiter':60,'ftol':1e-10})
    except Exception as error:failure=repr(error)
    accepted=False
    if result is not None:
        d,q,v,new=metrics(result.x);accepted=bool(result.success and d.min()>=.001 and q.min()>=.1 and v.min()>=.01)
    report=dict(passed=accepted,optimizer='scipy.optimize.minimize SLSQP',scipyVersion=scipy.__version__,
        freeVerticesZeroBased=free.tolist(),affectedCellsZeroBased=affected.tolist(),fixedBoundaryExact=True,
        evaluations=evaluations,optimizerSeconds=time.monotonic()-opt_start,engineBuilds=0,CFDRuns=0,physicalQualification=False,
        sourceSHA256={str(f):sha(f) for f in [source/'exact_mesh.npz',base/'rocket_geometry01/report.json',ROOT/'docs/validation/v200/LOCAL_MESH_PLAN.md']})
    if result is not None:report.update(optimizerSuccess=bool(result.success),message=str(result.message),iterations=int(result.nit),
        minimumLocalDeterminant=float(d.min()),minimumLocalSICN=float(q.min()),minimumVolumeRatio=float(v.min()),
        originalPointsM=basepoints.tolist(),newPointsM=new.tolist(),maximumDisplacementM=float(np.linalg.norm(new-basepoints,axis=1).max()))
    if failure:report['error']=failure
    if accepted:
        # 模块二：固定连接的新候选；全部原资产保持，保存采样质心的真实变化。
        case=base/'rocket03';shutil.copytree(source,case);p=p.copy();p[free]=new;g['points']=p
        np.savez(case/'exact_mesh.npz',**g)
        (case/'constant/polyMesh/points').write_text(header('points','vectorField')+str(len(p))+'\n(\n'+'\n'.join('('+' '.join(format(x,'.17g') for x in row)+')' for row in p)+'\n)\n')
        mapping=json.loads((case/'wall_mapping.json').read_text());ids=np.array(mapping['nativeFaceIDsZeroBased']);owners=np.array(mapping['nativeOwnerCellsZeroBased'])
        samples=p[c[owners]].mean(axis=1);center=np.array(mapping['wallCentersM']);tri=p[g['faces'][ids]]
        sf=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm=-sf/np.linalg.norm(sf,axis=1)[:,None];delta=samples-center;distance=np.einsum('ij,ij->i',delta,norm)
        assert np.all(distance>0);mapping.update(samplePositionsM=samples.tolist(),normalDistanceM=distance.tolist(),
            tangentialOffsetM=np.linalg.norm(delta-distance[:,None]*norm,axis=1).tolist())
        save(case/'wall_mapping.json',mapping)
        report['inputSHA256']={str(f):sha(f) for f in case.rglob('*') if f.is_file()}
    report['elapsedSeconds']=time.monotonic()-start;save(out/'report.json',report)
    print(json.dumps({k:report.get(k) for k in ['passed','message','iterations','evaluations','minimumLocalDeterminant','minimumLocalSICN','maximumDisplacementM','optimizerSeconds','elapsedSeconds','error']}))


if __name__=='__main__':main()
