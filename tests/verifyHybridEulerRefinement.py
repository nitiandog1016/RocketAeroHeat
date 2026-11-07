# -*- coding: utf-8 -*-
"""细网格独立六四面体/壁面积/均匀初态核验；不复用原生体积作为真值。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from readEilmerNative import read_grid,read_field,sha
from hybridEulerEvolution import audit_evolution_inputs


def verify(case,prep,out):
    case,prep,out=map(Path,(case,prep,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verify.py')
    info,g,m=audit_evolution_inputs(case,prep);v=g['points'];c=g['cells'];total=0.;minimum=float('inf')
    for tet in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
        a,b,d,e=[v[c[:,i]] for i in tet];part=np.einsum('ij,ij->i',b-a,np.cross(d-a,e-a))/6
        minimum=min(minimum,float(part.min()));total+=float(part.sum())
    exact=.04*(.6-.24*np.tan(np.radians(10)))
    assert minimum>0 and abs(total-exact)<2e-14
    surface=json.loads((case/'surface.json').read_text());xyz=v[np.array(surface['facesOneBased'])-1]
    area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
    expected=.04*(.6+.4/np.cos(np.radians(10)))
    assert len(area)==480 and abs(area.sum()-expected)<1e-14
    f=read_field(case,0,9600);inp=json.loads((case/'native-input-state.json').read_text())
    for name,expected_value in [('p',5000.),('T',300.),('rho',inp['rho']),('e',inp['e']),('vel.x',inp['velx']),('vel.y',0.),('vel.z',0.)]:
        assert np.all(f[name]==expected_value)
    result=dict(passed=True,minimumTetrahedronVolumeM3=minimum,volumeM3=total,
        analyticVolumeDifferenceM3=total-exact,wallAreaM2=float(area.sum()),wallAreaDifferenceM2=float(area.sum()-expected),
        cells=9600,wallTriangles=480,initialStateExact=True,physicalQualification=False,
        inputSHA256={str(case/n):sha(case/n) for n in ['mesh.su2','surface.json','rah-air.gas','lmrsim/config','lmrsim/control']})
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','preparation','out']:p.add_argument('--'+k,required=True)
    a=p.parse_args();verify(a.case,a.preparation,a.out)
