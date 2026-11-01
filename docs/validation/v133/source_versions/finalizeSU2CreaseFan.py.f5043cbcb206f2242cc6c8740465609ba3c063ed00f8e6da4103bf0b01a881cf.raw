"""只读保存的网格：归档外域失败质量、补齐工作台输入元数据，不重跑网格器。"""
from pathlib import Path
import json,os,time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.spatial import cKDTree
from optimizeSU2FanExterior import sicn
from prepareSU2CheckpointContinuation import ROOT,read,sha
BASE=ROOT/'outputs/v133_validation'

def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')

def main():
    target=BASE/'outer_result.json';assert not target.exists()
    source=BASE/'fan';rawpath=BASE/'optimized_recovered/optimized_raw.mat'
    raw=loadmat(str(rawpath),squeeze_me=True);old=loadmat(str(source/'hybrid_volume.mat'))
    tags=raw['nodeTags'].astype(int);remap=np.full(int(max(tags))+1,-1,int);remap[tags]=np.arange(len(tags))
    T=remap[raw['tetrahedronNodeTags'].astype(int)];P=raw['vertices'];q=sicn(P,T)
    oldP=old['vertices'];oldT=old['tetrahedra'].astype(int)-1;oldQ=sicn(oldP,oldT)
    bad=oldT[oldQ<.1];inner=np.unique(old['sourceFaces'].astype(int)-1+32*9878)
    near,ids=cKDTree(P).query(oldP[bad].reshape(-1,3));assert max(near)<1e-14
    existing={tuple(row) for row in np.sort(T[q<.1],axis=1)}
    retained=[tuple(row) in existing for row in np.sort(ids.reshape(-1,4),axis=1)]
    report=dict(nativeOptimizationCalls=2,nativeTerminalExitCodes=[1,1],secondFailure='declared 900000 cell budget exceeded before final export',
        firstOutputNotPersisted=True,secondRawResultPreserved=True,sourceCells=len(oldT)+len(old['wedges']),
        candidateCells=len(T)+len(old['wedges']),declaredMaximumCells=900000,
        originalMinimumSICN=float(min(oldQ)),candidateMinimumSICN=float(min(q)),
        originalBelow01=int(sum(oldQ<.1)),candidateBelow01=int(sum(q<.1)),
        oldLowQualityCellsExactlyRetained=retained,
        originalBadTets=bad.tolist(),badTetCenters=oldP[bad].mean(axis=1).tolist(),
        badTetInnerInterfaceVertexCounts=[int(sum(np.isin(row,inner))) for row in bad],
        qualityGatePassed=bool(min(q)>=.1 and not np.any(q<.1)),cellBudgetPassed=len(T)+len(old['wedges'])<=900000,
        candidateAdopted=False,defaultChanged=False,nativeFlowRuns=0,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/optimizeSU2FanExterior.py',rawpath,source/'hybrid_volume.mat',
            BASE/'optimized/mesh_process.json',BASE/'optimized_recovered/mesh_process.json']})
    write(target,report)
    # 首次几何输入未包含主工作台要求的optionsJSON，只重导出元数据到独立目录。
    folder=BASE/'fan_workbench';assert not folder.exists();folder.mkdir();start=time.monotonic()
    upstream=ROOT/'outputs/v120_validation/mesh32/mesh_input.mat'
    files=[source/'mesh_process.json',source/'geometry_audit.json',source/'hybrid_volume.mat',source/'mesh_report.json',upstream,Path(__file__).resolve()]
    locks={str(p):sha(p) for p in files}
    assert read(source/'geometry_audit.json')['passed']
    os.link(str(source/'hybrid_volume.mat'),str(folder/'hybrid_volume.mat'))
    m=read(source/'mesh_report.json');m['geometrySourceDirectory']=str(source)
    m['metadataOnlyReexport']=True;m['outerOptimizationAdopted']=False
    write(folder/'mesh_report.json',m)
    data=loadmat(str(upstream))
    savemat(str(folder/'mesh_input.mat'),dict(V=old['sourceVertices'],F=old['sourceFaces'],optionsJSON=data['optionsJSON'],
        constructionOperation=m['constructionOperation'],geometrySourceSHA256=sha(source/'hybrid_volume.mat')),do_compression=True)
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    receipt=dict(status='terminal',exitCode=0,timedOut=False,inputStable=True,sourceStable=True,
        operation='metadata-only-existing-crease-fan-reexport',geometryRecomputed=False,sourceHashes=locks,sourceStableAfter=stable,
        input=str(folder/'mesh_input.mat'),inputSHA256=sha(folder/'mesh_input.mat'),elapsedSeconds=time.monotonic()-start,
        outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file()})
    write(folder/'mesh_process.json',receipt)
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','originalBadTets','badTetCenters','oldLowQualityCellsExactlyRetained']},indent=2))

if __name__=='__main__':main()
