"""完整真实棱柱扇的独立拓扑/体积及局部对偶几何检查，只读不求流场。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from auditSU2HybridTopology import audit
from reviewSU2SSTLocalDiffusion import local_geometry,neighbors_of
from prepareSU2CheckpointContinuation import ROOT,read,sha

def main():
    start=time.monotonic();folder=ROOT/'outputs/v133_validation/fan';out=folder/'geometry_audit.json';assert not out.exists()
    receipt=read(folder/'mesh_process.json');assert receipt['exitCode']==0 and receipt['sourceStable'] and receipt['inputStable']
    for name in ['hybrid_volume.mat','mesh_report.json']:assert sha(folder/name)==receipt['outputSHA256'][name]
    d=loadmat(str(folder/'hybrid_volume.mat'));old=loadmat(str(ROOT/'outputs/v120_validation/mesh32/hybrid_volume.mat'))
    r=read(folder/'mesh_report.json');report=audit(d)
    assert np.array_equal(d['sourceVertices'],old['sourceVertices']) and np.array_equal(d['sourceFaces'],old['sourceFaces'])
    assert np.array_equal(d['tetrahedra'],old['tetrahedra']) and np.array_equal(d['farfield'],old['farfield'])
    T=d['tetrahedra'].astype(int)-1
    assert np.array_equal(d['vertices'][T],old['vertices'][T])
    baseline=read(ROOT/'outputs/v132_validation/local_diffusion.json')['geometry']
    rim=r['rimNodeIDs'];sid=np.array(r['sideNodeIDs']);bid=np.array(r['baseNodeIDs'])
    correspond=[]
    for key,value in baseline.items():
        j=rim.index(value['wallParentID']);k=value['normalLayer']
        for side,ids in [('side',sid),('base',bid)]:correspond.append((key,side,int(ids[k,j])))
    focus=set(row[2] for row in correspond);P=d['vertices'];cells=[d['wedges'].astype(int)-1,T]
    adj=neighbors_of(cells,focus,len(P));core=focus|set(j for v in adj.values() for j in v)
    assert min(core)>=len(d['sourceVertices'])
    g=local_geometry(P,cells,core)
    assert all(v['closure']<=1e-10 and v['linearError']<=1e-10 for v in g.values())
    result=[]
    for key,side,i in correspond:
        v=g[i];c=np.sum(v['edges']*v['normals'],axis=1)/(np.linalg.norm(v['edges'],axis=1)*np.linalg.norm(v['normals'],axis=1))
        result.append(dict(sourceNodeID=int(key),branch=side,candidateNodeID=i,
            originalGGLinearError=baseline[key]['ggLinearGeometryError'],newGGLinearError=v['ggLinearError'],
            originalWLSCondition=baseline[key]['wlsCondition'],newWLSCondition=v['condition'],
            newWLSLinearError=v['linearError'],newClosureRelative=v['closure'],newMinimumAlignmentCosine=float(min(c))))
    report.update(localGeometry=result,gradientNodes=len(core),originalSurfaceExactlyPreserved=True,
        originalExteriorExactlyPreserved=True,outerTetsBelowSICN01=r['tetsBelowSICN01'],
        sideFirstNormalAlignment=r['sideFirstNormalAlignment'],baseFirstNormalAlignment=r['baseFirstNormalAlignment'],
        minimumPrismJacobian=r['minimumPrismJacobian'],elapsedSeconds=time.monotonic()-start,
        localComparison=dict(originalMaximumGGLinearError=max(v['ggLinearGeometryError'] for v in baseline.values()),
            candidateMaximumGGLinearError=max(row['newGGLinearError'] for row in result),
            candidateMaximumWLSCondition=max(row['newWLSCondition'] for row in result)),
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/auditSU2HybridTopology.py',
            ROOT/'scripts/reviewSU2SSTLocalDiffusion.py',folder/'hybrid_volume.mat',folder/'mesh_process.json',folder/'mesh_report.json',
            ROOT/'outputs/v132_validation/local_diffusion.json']})
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['localGeometry','sourceSHA256']},indent=2))

if __name__=='__main__':main()
