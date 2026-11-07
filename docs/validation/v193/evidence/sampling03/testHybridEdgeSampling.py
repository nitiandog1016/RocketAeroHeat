# -*- coding: utf-8 -*-
"""既存原生细场与真实原壁几何的采样合同；明确后者不是原生流场验证。"""
from pathlib import Path
import argparse
import copy
import json
import shutil
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hybridEdgeSampling import candidate,map_original_wall
from readHybridEulerNative import read_wall_frame
from readEilmerNative import read_grid,read_field,compare_source,sha

ROOT=Path(__file__).resolve().parents[1]


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__),ROOT/'scripts/hybridEdgeSampling.py']
    for p in sources:shutil.copyfile(p,out/p.name)
    base=ROOT/'outputs/v192_validation';case=base/'fine01'
    prep=base/'fine_prepare01/report.json';run=base/'fine_run01/report.json'
    frame,info=read_wall_frame(case,8,prep,run);surface=json.loads((case/'surface.json').read_text())
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,_=compare_source(grid,case/'mesh.su2')
    raw=read_field(case,8,info['cells']);record=candidate(frame,surface,grid,mapping,raw)
    # 独立规律：120x2个底面，各拆2个三角面；每列40个单元。
    expected=np.repeat(np.arange(240)*40,2)
    assert np.array_equal(record['sampling']['nativeOwnerCellsZeroBased'],expected)
    positions=np.column_stack([raw['pos.'+a][expected] for a in 'xyz'])
    assert np.array_equal(record['sampling']['samplePositionsM'],positions)
    assert not record['engineeringHeatReady'] and all(v is None for v in record['wallHeatFluxWm2'])
    (out/'control_candidate.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError):rejected.append(name)
        else:raise AssertionError('未拒绝 '+name)
    for name,change in [('changed-coordinate',lambda d:d['vertices'][0].__setitem__(0,d['vertices'][0][0]+1e-3)),
                        ('missing-face',lambda d:d['facesOneBased'].pop()),
                        ('duplicate-face',lambda d:d['facesOneBased'].__setitem__(0,d['facesOneBased'][1])),
                        ('changed-id',lambda d:d['faceIDs'].__setitem__(0,2))]:
        bad=copy.deepcopy(surface);change(bad)
        reject(name,lambda bad=bad:candidate(frame,bad,grid,mapping,raw))
    for name,change in [('wrong-pressure',lambda d:d['fields']['pressurePa'].__setitem__(0,1.)),
                       ('wrong-temperature',lambda d:d['fields']['temperatureK'].__setitem__(0,1.)),
                       ('wrong-velocity',lambda d:d['fields']['velocityMps'][0].__setitem__(0,1.)),
                       ('no-eos-audit',lambda d:d.update(algebraicInterfacePassed=False))]:
        bad=copy.deepcopy(frame);change(bad)
        reject(name,lambda bad=bad:candidate(bad,surface,grid,mapping,raw))
    bad=copy.deepcopy(mapping);bad['wall']['nativeRows'][0]=bad['wall']['nativeRows'][1]
    reject('duplicate-native-map',lambda:candidate(frame,surface,grid,bad,raw))
    # 模块二：真实原壁三角输入走同一映射器；仅源网格，不冒称Eilmer已准备。
    meshcase=ROOT/'outputs/v193_validation/mesh01';mesh=np.load(meshcase/'mesh.npz')
    P,T,F,owner=[mesh[k] for k in ['vertices','tetrahedra','wall','wallOwners']]
    tri=F[:,[0,2,1]];n=len(F)
    g=dict(points=P,cells=T,faces=tri,boundaries={'wall':dict(ids=np.arange(n),signs=np.ones(n,int),owners=owner)})
    surf=json.loads((meshcase/'surface.json').read_text())
    mapped=map_original_wall(surf,g,dict(originalNodes=tri.tolist(),nativeRows=np.arange(n)),P[T].mean(1))
    csv=np.loadtxt(meshcase/'wall_sampling_geometry.csv',delimiter=',',skiprows=1)
    assert np.array_equal(mapped['nativeOwnerCellsZeroBased'],csv[:,1])
    error=float(np.max(abs(np.array(mapped['normalDistanceM'])-csv[:,3])));assert error<1e-14
    badpos=P[T].mean(1);badpos[owner[0]]=P[F[0]].mean(0)
    reject('not-fluid-side',lambda:map_original_wall(surf,g,dict(originalNodes=tri.tolist(),nativeRows=np.arange(n)),badpos))
    report=dict(passed=True,nativeControlFaces=480,nativeOwnerExact=True,nativePositionExact=True,
        actualRocketGeometryOnlyFaces=19752,actualRocketNativePrepared=False,actualRocketSamplingMaxDifferenceM=error,
        rejections=rejected,engineeringHeatReady=False,newPhysicalQualification=False,
        sourceSHA256={str(p):sha(p) for p in sources+[meshcase/'mesh.npz',meshcase/'surface.json']},
        elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);test(p.parse_args().out)
