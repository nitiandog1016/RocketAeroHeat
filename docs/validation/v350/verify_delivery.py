"""独立回读压缩包、STL 和实际工作台导出表，按原面核对数值、单位与功率。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np


def verify(ui_folder,require_trace):
    root=Path(__file__).resolve().parents[3]
    case=root/'cases/simple_launcher'
    digest=lambda data:hashlib.sha256(data).hexdigest()
    definition=json.loads((case/'geometry.json').read_bytes())
    raw=(case/'launcher-m.stl').read_bytes()
    assert digest(raw)==definition['fileSHA256']
    count=int.from_bytes(raw[80:84],'little')
    assert len(raw)==84+50*count
    dtype=np.dtype([('normal','<f4',3),('corners','<f4',(3,3)),('attribute','<u2')])
    vertices=np.frombuffer(raw,offset=84,dtype=dtype)['corners'].astype(float)
    area=np.linalg.norm(np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0]),axis=1)/2
    assert np.all(area>0)
    with ZipFile(case/'wind_tunnel.zip') as z:
        summary=json.loads(z.read('summary.json'));data=z.read('wall.csv')
        assert digest(data)==summary['wallSHA256']
        with z.open('wall.csv') as f:stored=np.loadtxt(f,skiprows=1,delimiter=',')
    with (ui_folder/'wall_pair.csv').open(encoding='utf-8-sig',newline='') as f:display=list(csv.DictReader(f))
    assert len(display)==count and list(map(int,(r['face_id'] for r in display)))==list(range(1,count+1))
    assert np.array_equal(stored[:,0],np.arange(1,count+1))
    np.testing.assert_allclose(area,stored[:,1],rtol=0,atol=1e-15)
    for r in display:
        assert (r['pressure_unit'],r['heat_flux_unit'])==('kPa','W/m²')
        assert r['pressure_method'] and r['heat_flux_method'] and int(r['synthetic'])==0
        assert r['pressure_availability']==r['heat_flux_availability']=='research' and int(r['joint_verified'])==0
        if require_trace:assert int(r['source_traceable'])==1
    p=np.array([float(r['pressure'])*1000 for r in display])
    q=np.array([float(r['heat_flux']) for r in display])
    np.testing.assert_allclose(p,stored[:,2],rtol=0,atol=1e-8)
    np.testing.assert_allclose(q,stored[:,3],rtol=0,atol=1e-13)
    powers={}
    for g in definition['regions']:
        s=slice(g['firstFace']-1,g['lastFace'])
        powers[g['name']]=float(-area[s]@q[s])
        assert abs(powers[g['name']]-summary['regions'][g['name']]['outwardW'])<1e-12
    assert abs(sum(powers.values())-summary['outwardW'])<1e-12
    return dict(verifiedRows=count,verifiedRegions=len(powers),outwardW=sum(powers.values()),
                maxDisplayPressureDifferencePa=float(max(abs(p-stored[:,2]))),
                maxDisplayFluxDifferenceWm2=float(max(abs(q-stored[:,3]))),
                resultSHA256=summary['wallSHA256'],sourceTraceable=all(int(r['source_traceable'])==1 for r in display))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ui_folder',type=Path)
    parser.add_argument('--require-trace',action='store_true')
    args=parser.parse_args();print(json.dumps(verify(args.ui_folder,args.require_trace),ensure_ascii=False,indent=2))
