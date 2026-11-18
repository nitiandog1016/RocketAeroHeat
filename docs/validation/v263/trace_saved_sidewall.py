"""核对同一原面上的已存量与上游缺失记录。只读取冻结的工程结果。"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from check_saved_surfaces import FILES, ROOT


def trace(output):
    snapshot = ROOT / 'outputs/v176_validation/main02/engineering_data.json'
    native = ROOT / 'outputs/v104_validation/baseline02/current_case.mat'
    for path in (snapshot, native):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == FILES[path.relative_to(ROOT).as_posix()]
    data = json.loads(snapshot.read_bytes())
    mesh = loadmat(native, variable_names=['current'], simplify_cells=True)['current']
    ledger = data['evidenceLedger']
    assert ledger['sourceFrame']['sha256'] == FILES[native.relative_to(ROOT).as_posix()]
    assert np.array_equal(data['vertices'], mesh['V']) and np.array_equal(data['faces'], mesh['F']-1)
    area = mesh['area']
    assert np.array_equal(data['faceArea'], area) and len(area) == 19752
    fields = {f['key']: f for f in data['fields']}
    evidence = {f['key']: f for f in ledger['fields']}
    coverage = mesh['heatingCoverage']
    reason = '上游外缘热力状态未解/域外占位影响'
    reason_index = list(coverage['reasonNames']).index(reason)
    selected = (mesh['pressure']['fieldStatus'] == 1) & (coverage['reasonCode'] == reason_index)
    json_selected = ((np.asarray(fields['pressure']['status']) == 1) &
                     (np.asarray(fields['heatFlux']['status']) == 0) & (np.asarray(data['missingReason']) == reason))
    assert np.array_equal(selected, json_selected)
    ids = np.flatnonzero(selected)
    assert len(ids) == 936 and abs(math.fsum(area[ids])-0.09570781817031267) < 1e-12
    assert np.all(np.asarray(data['region'])[ids] == 'sidewall')

    # 逐量比较原始数据、前端字段身份及独立验证状态。给定壁温单列。
    statuses = {'pressure': (1, 1), 'edgeTemperature': (1, 1), 'recoveryTemperature': (0, 0),
                'heatFlux': (0, 0), 'wallInput': (2, 3), 'wallSolved': (0, 0)}
    for key, (display, available) in statuses.items():
        assert np.all(np.asarray(fields[key]['status'])[ids] == display)
        assert np.all(np.asarray(evidence[key]['availability'])[ids] == available)
        assert np.all(np.asarray(evidence[key]['qualification'])[ids] == 0)
        if display == 0:
            assert all(fields[key]['values'][i] is None for i in ids)
    edge = mesh['boundaryLayerEdge']
    heat = mesh['heating']
    for key, raw in [('pressure',mesh['pressure']['staticPressure']/1000), ('edgeTemperature',edge['temperature']),
                     ('wallInput',heat['wallTemperatureByFace']), ('heatFlux',heat['heatFlux']/1000),
                     ('recoveryTemperature',heat['recoveryTemperature'])]:
        np.testing.assert_allclose(np.asarray(fields[key]['values'], dtype=float)[ids],raw[ids],
                                   rtol=1e-14,atol=0,equal_nan=True)
    assert np.all(heat['wallTemperatureByFace'][ids] == 300)
    path = mesh['surfaceStreamlines']
    regime = mesh['continuumRegime']
    assert path['velocitySource'] == edge['method']
    assert np.all(regime['downstreamStateUnresolvedFace'][ids]) and np.all(edge['thermalStateValidFace'][ids])
    direct = regime['directStateUnresolvedFace']
    rows = []
    for face in ids:
        node = face
        hops = 0
        while not direct[node]:
            parent_id = int(path['predecessorFace'][node])
            assert 1 <= parent_id <= len(area) and hops < len(area)
            parent = parent_id-1
            assert path['connectedFace'][node] and path['connectedFace'][parent]
            assert path['pathCost'][parent] < path['pathCost'][node]
            node = parent
            hops += 1
        rows.append(dict(faceID=int(face)+1,area_m2=area[face],
                         wallPressure_kPa=fields['pressure']['values'][face],
                         edgeTemperature_K=edge['temperature'][face],
                         edgeMethodIndex=fields['edgeTemperature']['methodIndex'][face],
                         wallInput_K=heat['wallTemperatureByFace'][face],recoveryTemperature_K=None,
                         heatIntoWall_kW_m2=None,solvedSolidTw_K=None,
                         nearestDirectGapFaceID=int(node)+1,upstreamLinks=hops))

    gap_ids = np.unique([r['nearestDirectGapFaceID'] for r in rows])
    hops = [r['upstreamLinks'] for r in rows]
    summary = dict(faceIDsStartAt=1,methodIndicesStartAt=0,totalOriginalFaces=len(area),
                   totalOriginalAreaM2=math.fsum(area),selectedFaces=len(ids),selectedAreaM2=math.fsum(area[ids]),
                   uniqueNearestDirectGaps=len(gap_ids),minLinks=min(hops),maxLinks=max(hops),
                   nearestGapsFlaggedAsPlaceholders=int(np.count_nonzero(edge['diagnosticPlaceholderFace'][gap_ids-1])),
                   reason=reason,edgeMethodNames=fields['edgeTemperature']['methodNames'],
                   sourceFrame=ledger['sourceFrame'],snapshotSHA256=FILES[snapshot.relative_to(ROOT).as_posix()],
                   newJointAreaM2=0,qualifiedAreaM2=0)
    output.mkdir(exist_ok=False)
    with (output/'faces.csv').open('w',encoding='utf-8',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(len(ids),len(gap_ids),min(hops),max(hops),summary['nearestGapsFlaggedAsPlaceholders'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('new_output_folder',type=Path)
    trace(parser.parse_args().new_output_folder)
