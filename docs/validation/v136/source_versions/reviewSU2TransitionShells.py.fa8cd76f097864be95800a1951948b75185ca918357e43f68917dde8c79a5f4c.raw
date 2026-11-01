"""只读真实原壁、既有外壳和终止接口的全表面自交，不重生成网格。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from runCGALSurfaceAudit import audit_surface
from meshSU2HybridTransition import prepare,boundary
from prepareSU2CheckpointContinuation import ROOT,sha,read

def quad_fan(P,T,Q):
    centers=P[Q].mean(axis=1);apex=np.arange(len(P),len(P)+len(Q));tri=[T]
    for j in range(4):tri.append(np.c_[Q[:,j],Q[:,(j+1)%4],apex])
    return np.r_[P,centers],np.concatenate(tri)

def main():
    start=time.monotonic();source=ROOT/'outputs/v133_validation/fan_workbench';out=ROOT/'outputs/v136_validation/shells'
    assert not out.exists();out.mkdir(parents=True)
    files=[source/n for n in ['hybrid_volume.mat','mesh_report.json','mesh_process.json','mesh_input.mat']]+[Path(__file__).resolve(),ROOT/'scripts/meshSU2HybridTransition.py']
    hashes={str(p):sha(p) for p in files};d=loadmat(str(source/'hybrid_volume.mat'));P=d['vertices'];F=d['sourceFaces'].astype(int)-1
    W=d['wedges'].astype(int)-1;T,Q=boundary(W,F);assert len(Q)==0
    meta=dict(sourceSHA256=hashes,physicalQualification=False)
    reports={}
    reports['original_wall']=audit_surface(d['sourceVertices'],F,out/'original_wall',meta)
    reports['fan_outer']=audit_surface(P,T,out/'fan_outer',meta)
    P,W,V,F,far,T,Q,info,options,keep=prepare(source,False)
    combined,faces=quad_fan(P,T,Q)
    reports['terminated_interface']=audit_surface(combined,faces,out/'terminated_interface',dict(meta,quadrilateralRepresentation='four-face-center-triangles',originalTriangularFaces=len(T),originalQuadFaces=len(Q)))
    stable={name:sha(Path(name))==h for name,h in hashes.items()};assert all(stable.values())
    summary=dict(reports=reports,sourceStableAfter=stable,seconds=time.monotonic()-start,
        allSelfIntersectionFree=all(r['selfIntersectionFree'] for r in reports.values()),native3DMeshRuns=0,nativeFlowRuns=0,newQualifiedArea=0)
    (out/'report.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8')
    print(json.dumps({name:dict(faces=r['faces'],pairs=r['pairCount'],seconds=r['seconds']) for name,r in reports.items()}))
if __name__=='__main__':main()
