# -*- coding: utf-8 -*-
"""真实原三角面初场接口/来源拒绝和完整场守卫；不重复CFD小控制。"""
from pathlib import Path
import argparse
import copy
import json
import shutil
import sys
import time
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hybridRocketEuler import check_inputs,check_admission,read_wall,SavedFieldGuard,BUDGET,sha,save


def verify(case, preparation, out):
    start=time.monotonic();case,preparation,out=map(Path,(case,preparation,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verify.py')
    r=check_inputs(case,preparation);frame=read_wall(case,0,preparation)
    assert len(frame['faceIDs'])==19752 and frame['stateIdentity']=='uniform-initial-not-solved'
    assert all(v is None for v in frame['wallHeatFluxWm2']) and not frame['engineeringHeatReady']
    original=json.loads((case/'native-input-state.json').read_text())
    for key,name in [('pressurePa','p'),('temperatureK','T'),('densityKgm3','rho'),('internalEnergyJkg','e')]:
        np.testing.assert_allclose(frame['fields'][key],original[name],rtol=2e-13,atol=2e-13)
    np.testing.assert_allclose(frame['fields']['velocityMps'],np.tile([original[k] for k in ['velx','vely','velz']],(19752,1)),rtol=0,atol=0)
    rejected=[]
    for n in [-1,1]:
        try:read_wall(case,n,preparation)
        except ValueError:rejected.append('snapshot_'+str(n))
        else:raise AssertionError('缺运行来源被放行')
    a=dict(schema='RAH.HybridRocketEulerAdmission.v1',profile='full-wall-euler-400-v1',budget=BUDGET,
           temperatureRangeK=[200.,6000.],evidenceSHA256={str(preparation):sha(preparation)})
    check_admission(a)
    for key,value in [('profile','coarse'),('budget',dict(BUDGET,maxStep=401)),('temperatureRangeK',[80.,6000.]),
                      ('evidenceSHA256',{str(preparation):'0'*64})]:
        bad=copy.deepcopy(a);bad[key]=value
        try:check_admission(bad)
        except ValueError:rejected.append(key)
        else:raise AssertionError('错误准入被放行')
    # 真实配置文件保持不变，只替换读取结果验证完整配置/控制比较。
    original_read=Path.read_text
    for name,key,value in [('config','viscous',True),('config','flux_calculator','roe'),('control','max_step',401)]:
        target=case/'lmrsim'/name;changed=json.loads(target.read_text());changed[key]=value
        def altered(path,*args,**kwargs):
            return json.dumps(changed) if path.resolve()==target.resolve() else original_read(path,*args,**kwargs)
        with patch.object(Path,'read_text',altered):
            try:check_inputs(case,preparation)
            except ValueError:rejected.append('config_'+key)
            else:raise AssertionError('实际配置变动被放行')
    guard=SavedFieldGuard(case,64982,out/'initial_guard.json');assert guard() is None
    assert len(guard.rows)==1 and guard.rows[0]['invalidCells']==0
    # 完整/不完整/低温/非正状态的合成守卫检查，明确不作为物理解。
    synth=out/'guard_fixture';snap=synth/'lmrsim/snapshots/0001';snap.mkdir(parents=True)
    shutil.copyfile(case/'lmrsim/snapshots/fluid.metadata',snap.parent/'fluid.metadata')
    g=SavedFieldGuard(synth,2,out/'fixture_guard.json')
    data=np.ones((len(g.names),2));data[g.names.index('T')]=[300,199]
    raw=data.astype('<f8').tobytes();f=snap/'fluid-0000';f.write_bytes(raw[:10]);assert g() is None and not g.rows
    f.write_bytes(raw);bad=g();assert bad['invalidCells']==1 and bad['firstInvalidCellZeroBased']==1
    lowtemp=bad
    nextdir=snap.parent/'0002';nextdir.mkdir();data[g.names.index('T')]=[300,300];data[g.names.index('rho')][0]=-1
    (nextdir/'fluid-0000').write_bytes(data.astype('<f8').tobytes());bad=g();assert bad['invalidCells']==1 and bad['outsideTemperatureCount']==0
    result=dict(passed=True,originalFaces=19752,originalAreaM2=frame['originalAreaM2'],rejected=rejected,
        initialEOS=frame['fullVolumeEOSErrors'],lowTemperatureGuard=lowtemp,nonpositiveGuard=bad,
        scope='actual original-wall initial mapping; manifest/config substitutions; synthetic complete-file guard; no CFD',
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,
        sourceSHA256={str(p.resolve()):sha(p) for p in [Path(__file__),preparation,Path(__file__).parents[1]/'scripts/hybridRocketEuler.py']})
    save(out/'report.json',result)
    print(json.dumps(dict(passed=True,rejections=len(rejected),elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['case','preparation','out']:p.add_argument('--'+key,required=True)
    a=p.parse_args();verify(a.case,a.preparation,a.out)
