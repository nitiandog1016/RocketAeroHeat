"""主工作台的小型解析壳网格。固定连接与尺寸，仅验证新增三维流程。"""
from pathlib import Path
import json
import time

import numpy as np
from scipy.io import loadmat, savemat

from auditSU2HybridTopology import audit
from meshSU2HybridVolume import wedge_jacobian_min
from reviewSU2ExplicitQuality import quality
from su2ThermalRuntime import ROOT, sha
from runEmbeddedSU2 import atomic_write_json


def prepare():
    """复用原八面体8面顺序，完整外层按全局顶点序三分，避免接口悬挂。"""
    folder = ROOT/'outputs/v157_validation/fixture'
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source = ROOT/'outputs/v145_validation/near_wall_fixture/mesh_input.mat'
    original = loadmat(str(source))
    wall = original['V']
    faces = original['F'].astype(int)-1
    assert len(wall) == 6 and len(faces) == 8
    radii = np.array([1., 1.000002, 1.3, 2.])
    points = np.concatenate([wall*r for r in radii])
    wedges = np.concatenate([np.c_[faces+6*j, faces+6*(j+1)] for j in range(2)])
    tetrahedra = []
    for face in faces:
        a, b, c = sorted(face+12)
        tetrahedra += [[a,b,c,c+6], [a,b,b+6,c+6], [a,a+6,b+6,c+6]]
    tetrahedra = np.array(tetrahedra)
    p = points[tetrahedra]
    sign = np.linalg.det(p[:,1:]-p[:,:1])
    tetrahedra[sign<0] = tetrahedra[sign<0][:,[0,2,1,3]]
    data = dict(vertices=points, wedges=wedges+1, tetrahedra=tetrahedra+1,
                pyramids=np.empty((0,5), int), farfield=faces+19,
                sourceVertices=wall, sourceFaces=faces+1)
    topology = audit(data)
    assert abs(topology['totalFluidVolume']-28./3) < 2e-12
    jacobian = float(min(wedge_jacobian_min(points[wedges])))
    sicn = quality(points[tetrahedra])
    assert jacobian > 0 and min(sicn) >= .1
    report = dict(nodes=24, prisms=16, pyramids=0, tetrahedra=24, originalFaces=8,
        originalVertices=6, outerBoundaryFaces=8, normalLayers=2,
        normalHeightsM=(radii[:3]-1).tolist(),
        method='analytic-octahedron-prisms-and-staircase-tetrahedra',
        constructionOperation='analytic-octahedron-prisms-and-staircase-tetrahedra',
        minimumTetSICN=float(min(sicn)), tetsBelowSICN01=0,
        minimumPrismJacobian=jacobian, minimumTetVolume=topology['minimumTetVolume'],
        bodyArea=float(4*np.sqrt(3)), bodyVolume=4./3, farfieldVolume=32./3,
        originalFacesUnchanged=True, interfaceExactlyPreserved=True,
        qualityGatePassed=True, nearWallResolutionValidated=False,
        farfieldAdequacyValidated=False, physicsQualified=False, newQualifiedArea=0)
    savemat(str(folder/'hybrid_volume.mat'), data, do_compression=True)
    savemat(str(folder/'mesh_input.mat'), dict(V=wall, F=faces+1, optionsJSON=json.dumps({})))
    atomic_write_json(folder/'mesh_report.json', report)
    atomic_write_json(folder/'geometry_audit.json', topology)
    receipt = dict(status='terminal', exitCode=0, timedOut=False, inputStable=True, sourceStable=True,
        input=str(folder/'mesh_input.mat'), inputSHA256=sha(folder/'mesh_input.mat'),
        outputSHA256={n:sha(folder/n) for n in ('hybrid_volume.mat','mesh_report.json')},
        sourceSHA256={str(p):sha(p) for p in (Path(__file__).resolve(), source)},
        elapsedSeconds=time.monotonic()-started, analyticFluidVolume=28./3,
        nativeMesherRuns=0, physicsQualified=False, newQualifiedArea=0)
    atomic_write_json(folder/'mesh_process.json', receipt)
    print(json.dumps(dict(elapsedSeconds=receipt['elapsedSeconds'], minimumTetSICN=min(sicn),
                          volumeError=abs(topology['totalFluidVolume']-28./3))))


if __name__ == '__main__':
    prepare()
