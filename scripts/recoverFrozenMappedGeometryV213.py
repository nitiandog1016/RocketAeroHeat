"""模块职责：仅恢复已证的棱柱节点格式，并监督唯一原生控制；不改几何。"""
from pathlib import Path
import json,shutil,sys,os,time
import numpy as np
from runEilmerProcess import run,sha
from auditFrozenMappedGeometryV213 import PERM,volume
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v1'
CASE=BASE/'control02';PLAN=ROOT/'docs/validation/v213/GEOMETRY_RECOVERY.md'

def convert():
    start=time.monotonic();old=BASE/'control01';assert json.loads((BASE/'geometry_order_audit.json').read_text())['passed']
    CASE.mkdir(exist_ok=False)
    for n in ['job.lua','phase-air.lua','grid.lua','mapping.npz']:shutil.copyfile(old/n,CASE/n)
    lines=(old/'mesh.su2').read_text().splitlines();n=int(lines[1].split('=')[1])
    z=np.load(str(CASE/'mapping.npz'));assert n==len(z['cells'])==416
    for i in range(n):
        row=lines[i+2].split();assert row[0]=='13' and list(map(int,row[1:7]))==z['cells'][i].tolist()
        row[1:7]=[row[1+j] for j in PERM];lines[i+2]=' '.join(row)
    (CASE/'mesh.su2').write_text('\n'.join(lines)+'\n')
    d=json.loads((old/'donors.json').read_text())
    for r in d['donors']:r['vertices']=[r['vertices'][j] for j in PERM]
    d['vertexConvention']='native-VTK; official Gmsh permutation [1,0,2,4,3,5]'
    (CASE/'donors.json').write_text(json.dumps(d,indent=2)+'\n')
    np.savez(CASE/'native_connectivity.npz',cells=z['cells'][:,PERM])
    assert np.all(volume(z['points'][z['cells'][:,PERM]])>0)
    assert np.all(volume(np.array([r['vertices'] for r in d['donors']]))>0)
    report=dict(passed=True,seconds=time.monotonic()-start,permutation=PERM.tolist(),
        oldMappingSHA256=sha(old/'mapping.npz'),newMappingSHA256=sha(CASE/'mapping.npz'),
        originalCaseUnchanged=True,physicalGeometryChanged=False)
    assert report['oldMappingSHA256']==report['newMappingSHA256']
    (CASE/'conversion.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

def native(stage):
    build=json.loads((BASE/'build_result.json').read_text());binary=PACKAGE/'bin/probe'
    assert build['passed'] and sha(binary)==build['binarySHA256']
    if stage=='grid':name='prep_grid';budget=30;cmd=['lmr','prep-grids','--job=grid.lua'];record='prep_grid02'
    elif stage=='prepare':name='prep_sim';budget=30;cmd=[str(binary),'prepare'];record='prep_sim03'
    elif stage=='query':
        name='query';budget=60;cmd=['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary)];record='query02'
        assert json.loads((BASE/'prep_sim03/report.json').read_text())['passed']
    else:raise ValueError('stage')
    spent=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob(name+'*/report.json'))
    assert 0<budget-spent
    locks=[p for p in CASE.rglob('*') if p.is_file()]+[PLAN,Path(__file__).resolve(),binary,PACKAGE/'lib/prepsim.lua']
    r=run(cmd,CASE,BASE/record,budget-spent,32,locks,2048,2048)
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'],exit=r['nativeExitCode'],cumulativeSeconds=spent+r['elapsedSeconds'])))
    return r['passed']

if __name__=='__main__':
    raise ValueError('v213节点格式恢复已结束；复用control02及保留失败，禁止重复查询。')
    if sys.argv[1]=='convert':convert()
    else:raise SystemExit(0 if native(sys.argv[1]) else 1)
