# -*- coding: utf-8 -*-
"""用唯一原始优化输出回放最终批次选择，并核对原生外边界实际状态；不重跑优化。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from optimizeHybridRocketEulerMesh import appended_cohort
from readEilmerNative import sha


def verify(out):
    start=time.monotonic();root=Path(__file__).resolve().parents[1];base=root/'outputs/v194_validation'
    out=Path(out);out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_verify.py')
    old=np.load(root/'outputs/v193_validation/mesh01/mesh.npz');raw=np.load(base/'netgen01/optimized_raw.npz')
    recovered=np.load(base/'netgen_recovered02/new_cohort.npz')
    args=[raw['nodeTags'],raw['vertices'],raw['tetrahedronTags'],raw['tetrahedronNodeTags'],old['vertices'],old['tetrahedra'],26325]
    tags,conn=appended_cohort(*args)
    assert np.array_equal(tags,recovered['sourceElementTags'])
    assert np.array_equal(conn,recovered['sourceNodeTags'][recovered['tetrahedra']])
    rejected=[]
    for name,index,row in [('duplicate_element',2,0),('duplicate_node',0,0),('wrong_old_connectivity',3,0),
                           ('moved_old_node',1,0),('missing_new_sequence',2,-1)]:
        changed=list(args);changed[index]=args[index].copy()
        if name in ('duplicate_element','duplicate_node'):changed[index][row]=changed[index][1]
        elif name=='wrong_old_connectivity':changed[index][0,0]=changed[index][0,1]
        elif name=='moved_old_node':changed[index][0,0]+=.01
        else:changed[index][-1]+=2
        try:appended_cohort(*changed)
        except AssertionError:rejected.append(name)
        else:raise AssertionError('未拒绝错误批次: '+name)
    # 模块二：直接核对实际配置的边界动作和来流，不把标签正确等同实际边界正确。
    case=base/'rocket01';cfg=json.loads((case/'lmrsim/config').read_text())
    native=json.loads((case/'native-input-state.json').read_text());b=cfg['block_0']
    assert cfg['gas_model_file']=='rah-air.gas' and b['nboundaries']==2
    assert b['boundary_0']['pre_recon_action']==[dict(type='internal_copy_then_reflect')]
    action=b['boundary_1']['pre_recon_action'];assert len(action)==1 and action[0]['type']=='from_upwind_copy'
    fs=action[0]['flowstate']
    for k in ['p','T','velx','vely','velz']:assert fs[k]==native[k]
    assert fs['massf']==[1.] and fs['mu_t']==0 and fs['k_t']==0
    paths=[Path(__file__),root/'scripts/optimizeHybridRocketEulerMesh.py',base/'netgen01/optimized_raw.npz',
           base/'netgen_recovered02/new_cohort.npz',case/'lmrsim/config',case/'native-input-state.json']
    result=dict(passed=True,allNewCellsExact=True,newCells=len(tags),rejected=rejected,
        scope='saved raw cohort replay and native boundary configuration only; no optimizer or CFD rerun',
        nativeBoundaryActionsPassed=True,sourceSHA256={str(p):sha(p) for p in paths},
        physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=True,newCells=len(tags),rejections=len(rejected),elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
