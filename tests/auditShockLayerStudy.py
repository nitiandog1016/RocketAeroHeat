"""逐值核对工作台与冻结MAT：法向列、原网格、算法身份和未知热量不得改变。"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
from scipy.io import loadmat

from auditPressureDrivenSurfaceFlow import rows

ROOT=Path(__file__).resolve().parents[1]


def sha(file):
    return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def audit(folder,target):
    assert not target.exists()
    data=json.loads((folder/'shock-layer-data.json').read_text(encoding='utf8'))
    assert data['kind']=='shock-layer-study-v1' and data['productionQualified'] is False
    assert data['newProductionSolvedAreaM2']==0
    saved={}
    for name,source in data['sources'].items():
        file=ROOT/source['path'];assert sha(file)==source['sha256'];saved[name]=loadmat(file,simplify_cells=True)
    p,c,t=[saved[key] for key in ('parabolic','conic','translated')]
    groups=[[rows(p['cases'])[1]['column']],[rows(p['cases'])[0]['column']],
            rows(t['bodyFit']['result']['columns']),rows(t['heldOut']['columns'])]
    source_fields={'eta':'eta','position':'position','normalDistanceM':'normalDistanceM','jacobianRatio':'jacobianRatio'}
    state_fields={'pressurePa':'pressure','temperatureK':'temperature','densityKgM3':'density',
                  'speedMS':'speed','entropyJPerKgK':'entropy','sourceShockCoordinateM':'sourceShockCoordinate'}
    stations=0
    for k,(case,columns) in enumerate(zip(data['cases'],groups)):
        assert case['productionQualified'] is False and len(case['columns'])==len(columns)
        for j,(shown,col) in enumerate(zip(case['columns'],columns)):
            for shown_name,old in source_fields.items():
                np.testing.assert_array_equal(shown[shown_name],col[old])
            for shown_name,old in state_fields.items():
                np.testing.assert_array_equal(shown[shown_name],col['state'][old])
            for key in ('heatFluxWPerM2','wallTemperatureInputK','wallTemperatureSolvedK','edgeTemperatureK','recoveryTemperatureK'):
                assert shown[key]==[],key
            for key,old in [('pressure','pressureMethod'),('temperature','temperatureMethod'),('speed','speedMethod')]:
                assert shown['methods'][key]==col['state'][old]
            np.testing.assert_array_equal(shown['pickNormals'],np.tile(-col['shock']['normal'],(len(shown['eta']),1)))
            assert shown['productionQualified'] is False
            if k>=2:
                r=t['bodyFit']['result'] if k==2 else t['heldOut'];m=r['geometricMapping']
                assert shown['bodyResidualM']==r['residualM'][j]
                assert shown['originalFaceID']==m['originalFaceID'][j]
                np.testing.assert_array_equal(shown['barycentric'],m['barycentric'][j])
            else:
                assert shown['bodyResidualM']==[] and shown['originalFaceID']==[]
            stations+=len(shown['eta'])
    geometry=data['geometry'];V=t['geometry']['V'];F=t['geometry']['F'].astype(int)-1
    np.testing.assert_array_equal(geometry['vertices'],V);np.testing.assert_array_equal(geometry['faces'],F)
    focus=geometry['focus'];ids=np.asarray(focus['originalFaceIDs'],dtype=int)-1
    np.testing.assert_array_equal(np.asarray(focus['vertices'])[focus['faces']],V[F[ids]])
    assert geometry['closed'] is True and focus['closed'] is False
    assert data['fitSummary']['maximumM'][-1]==t['report']['maximumHeldOutResidualM']
    assert data['fitSummary']['toleranceM']==1e-5
    report=dict(status='passed-exact-frozen-MAT-to-workbench-geometry-columns-state-roles-and-method-audit',
                cases=len(groups),columns=sum(map(len,groups)),stations=stations,originalFaces=len(F),
                focusFaces=len(ids),dataSHA256=sha(folder/'shock-layer-data.json'),
                unknownHeatAndTemperaturesPreserved=True,productionQualified=False,
                newProductionSolvedAreaM2=0,browserVisualVerified=False)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path);parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();audit(args.folder,args.report)
