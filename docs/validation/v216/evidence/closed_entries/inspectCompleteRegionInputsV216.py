"""模块职责：只读记录完整母网格与现有原生差分输入，不生成新网格或流动场。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat,whosmat
from eilmerNKSAResourceBudget import ROOT,VOLUME_SHA,sha,save
BASE=ROOT/'outputs/v216_validation'
def main():
    start=time.monotonic();p=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(p)==VOLUME_SHA
    d=loadmat(str(p),variable_names=['vertices','sourceVertices','sourceFaces','wedges','pyramids','tetrahedra','farfield'])
    v=d['sourceVertices'];f=d['sourceFaces'].astype(np.int64)-1;pnt=d['vertices']
    far=d['farfield'].astype(np.int64)-1;cross=np.cross(pnt[far[:,1]]-pnt[far[:,0]],pnt[far[:,2]]-pnt[far[:,0]])
    normal=cross/np.linalg.norm(cross,axis=1)[:,None];center=pnt[far].mean(axis=1)
    report=dict(passed=True,fields=[dict(name=n,shape=list(s),dtype=t) for n,s,t in whosmat(str(p))],
        sourceBounds=[v.min(axis=0).tolist(),v.max(axis=0).tolist()],domainBounds=[pnt.min(axis=0).tolist(),pnt.max(axis=0).tolist()],
        wallCount=len(f),counts={k:len(d.get(k,[])) for k in ['wedges','pyramids','tetrahedra','farfield']},
        outerNormalBounds=[normal.min(axis=0).tolist(),normal.max(axis=0).tolist()],
        originalRingsXM=np.unique(v[:,0]).tolist(),motherSHA256=VOLUME_SHA,seconds=time.monotonic()-start)
    save(BASE/'inventory.json',report);print(json.dumps({k:report[k] for k in ['passed','sourceBounds','domainBounds','counts','seconds']}))
if __name__=='__main__':main()
