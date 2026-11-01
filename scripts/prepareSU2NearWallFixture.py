"""固定把旧小例六个首层节点改为2微米；保护其他几何，不调用网格器。"""
from pathlib import Path
import json,shutil
import numpy as np
from scipy.io import loadmat,savemat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha
from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from reviewSU2ExplicitQuality import quality

def main():
    src=ROOT/'outputs/v135_validation/fixture';out=ROOT/'outputs/v145_validation/near_wall_fixture';assert not out.exists()
    d=loadmat(str(src/'hybrid_volume.mat'));V=d['vertices'];wall=d['sourceVertices'];dist,ids=cKDTree(wall).query(V)
    moved=np.flatnonzero(abs(dist-.1)<1e-12);assert len(moved)==6
    before=V.copy();V[moved]=wall[ids[moved]]+(V[moved]-wall[ids[moved]])*(2e-6/.1)
    assert np.array_equal(V[:len(wall)],wall) and np.array_equal(V[d['farfield']-1],before[d['farfield']-1])
    g=audit(d);wj=float(min(wedge_jacobian_min(V[d['wedges']-1])));assert wj>0
    q=quality(V[d['tetrahedra']-1]);report=read(src/'mesh_report.json')
    report.update(minimumPrismJacobian=wj,minimumTetVolume=g['minimumTetVolume'],minimumTetSICN=float(min(q)),
        tetsBelowSICN01=int(sum(q<.1)),qualityGatePassed=bool(min(q)>=.1),noseRetainedHeightM=2e-6,
        normalHeightsM=[0,2e-6,.3,.6],initializationFixtureOnly=True)
    out.mkdir(parents=True);savemat(str(out/'hybrid_volume.mat'),{k:v for k,v in d.items() if not k.startswith('__')},do_compression=True)
    shutil.copyfile(str(src/'mesh_input.mat'),str(out/'mesh_input.mat'))
    (out/'mesh_report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'geometry_audit.json').write_text(json.dumps(g,indent=2)+'\n')
    r=dict(status='terminal',exitCode=0,timedOut=False,inputStable=True,sourceStable=True,input=str(out/'mesh_input.mat'),
        inputSHA256=sha(out/'mesh_input.mat'),outputSHA256={n:sha(out/n) for n in ['hybrid_volume.mat','mesh_report.json']},
        movedNodeIDsOneBased=(moved+1).tolist(),firstDistanceM=2e-6,geometryRegenerated=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),src/'hybrid_volume.mat',src/'mesh_report.json']})
    (out/'mesh_process.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
