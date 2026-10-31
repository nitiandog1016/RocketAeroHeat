"""用已确认的同网格/同配置原生检查点继续稳态研究，显式保存累计更新身份。"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTROL = {'ITER','OUTPUT_WRT_FREQ','WRT_RESTART_OVERWRITE','WRT_SURFACE_OVERWRITE',
           'OUTPUT_FILES','VOLUME_OUTPUT','HISTORY_OUTPUT','SCREEN_WRT_FREQ_INNER',
           'RESTART_SOL','SOLUTION_FILENAME','READ_BINARY_RESTART'}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def config(path):
    pairs=re.findall(r'^\s*([A-Z0-9_]+)\s*=\s*([^%\r\n]+)',path.read_text(encoding='utf-8-sig'),re.M)
    assert len(pairs)==len(dict(pairs)), '配置键重复'
    return {k:v.strip() for k,v in pairs}


def prepare(source,target):
    source=Path(source).resolve(); target=Path(target).resolve()
    source.relative_to(ROOT/'outputs'); target.relative_to(ROOT/'outputs')
    assert source!=target and not (target/'engine_process.json').exists()
    assert not (target/'initial_flow.csv').exists() and not (target/'continuation.json').exists()
    receipt=read(source/'engine_process.json')
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut'] and not receipt['cancelled'], '仅接受正常终态源任务'
    assert all(receipt['inputsStable'].values()) and receipt['runtimeFlavor']=='official'
    a=config(source/'case.cfg'); b=config(target/'case.cfg')
    assert a['SOLVER']=='NAVIER_STOKES' and a['KIND_TURB_MODEL']=='NONE' and a['FLUID_MODEL']=='IDEAL_GAS', '仅支持现有定比热层流NS检查点'
    assert a.get('TIME_DOMAIN','NO')=='NO' and a.get('MULTIZONE_MESH','NO')=='NO'
    assert {k:v for k,v in a.items() if k not in CONTROL}=={k:v for k,v in b.items() if k not in CONTROL}, '物理/数值配置不一致，拒绝继续'
    for name in ['case.cfg','mesh.su2']:
        assert sha(source/name)==receipt['inputSHA256'][str(source/name)],name
    assert sha(source/'mesh.su2')==sha(target/'mesh.su2'), '体网格不一致'
    for name in ['restart_flow.csv','surface_flow.vtk','history.csv']:
        assert sha(source/name)==receipt['outputSHA256'][name],name
    assert a['WRT_RESTART_OVERWRITE']=='NO' and a['WRT_SURFACE_OVERWRITE']=='NO', '缺少具名源时刻'
    matched=[]
    for p in source.glob('restart_flow_*.csv'):
        m=re.fullmatch(r'restart_flow_(\d{6,})\.csv',p.name)
        if not m:continue
        wall=source/('surface_flow_'+m[1]+'.vtk')
        if wall.exists():
            assert sha(p)==receipt['outputSHA256'][p.name] and sha(wall)==receipt['outputSHA256'][wall.name]
            if sha(p)==receipt['outputSHA256']['restart_flow.csv'] and sha(wall)==receipt['outputSHA256']['surface_flow.vtk']: matched.append(int(m[1]))
    assert matched, '没有确认的原生文件对'
    step=max(matched)
    hist=np.atleast_1d(np.genfromtxt(str(source/'history.csv'),names=True,delimiter=','))
    assert np.array_equal(hist['Inner_Iter'],np.arange(len(hist))) and step==int(hist['Inner_Iter'][-1])
    state=np.loadtxt(str(source/'restart_flow.csv'),delimiter=',',skiprows=1)
    assert state.shape[1]==9 and np.array_equal(state[:,0],np.arange(len(state))) and np.isfinite(state).all()
    rho=state[:,4]; internal=state[:,8]/rho-.5*np.sum((state[:,5:8]/rho[:,None])**2,axis=1)
    assert min(rho)>0 and min(internal)>0, '源场含非正状态'
    before=0; priorSeconds=float(receipt['elapsedSeconds']); locks=[source/'case.cfg',source/'mesh.su2',source/'engine_process.json',source/'restart_flow.csv',source/'surface_flow.vtk',source/'history.csv']
    if (source/'continuation.json').exists():
        old=read(source/'continuation.json');assert sha(source/'continuation.json')==receipt['inputSHA256'][str(source/'continuation.json')]
        before=int(old['completedUpdatesBeforeRun']);priorSeconds+=old['priorSuccessfulEngineSeconds'];locks.append(source/'continuation.json')
    manifest=read(target/'case_manifest.json')
    # 初猜说明必须更换为真实源场；不保留“C1初猜”冒充本轮输入。
    lineage=dict(schemaVersion=1,sourceDirectory=str(source),sourceSavedLocalIteration=step,
                 sourceCumulativeSavedIteration=before+step,completedUpdatesBeforeRun=before+step+1,
                 priorSuccessfulEngineSeconds=priorSeconds,sourceStateSHA256=sha(source/'restart_flow.csv'),
                 inputSourceHashes={str(p):sha(p) for p in locks},samePhysicsAndNumerics=True,sameMesh=True,
                 pseudoTimeHistoryRestarts=True,physicalTimeContinuation=False,physicsQualified=False,
                 method='exact native conservative checkpoint; no interpolation or clipping')
    # 所有检查完成后才写新目录；共享输入绝不能被输出覆写。
    os.link(str(source/'restart_flow.csv'),str(target/'initial_flow.csv'))
    text=(target/'case.cfg').read_text(encoding='utf-8-sig')
    for k,value in [('RESTART_SOL','YES'),('SOLUTION_FILENAME','initial_flow'),('READ_BINARY_RESTART','NO')]:
        text,n=re.subn(r'^\s*'+k+r'\s*=.*$',k+'= '+value,text,flags=re.M)
        if not n:text+='\n'+k+'= '+value+'\n'
        assert n<=1
    (target/'case.cfg').write_text(text,encoding='utf-8')
    manifest['initialization']=lineage
    manifest['continuation']=lineage
    (target/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    (target/'continuation.json').write_text(json.dumps(lineage,indent=2)+'\n',encoding='utf-8')
    return lineage


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--target',required=True);a=p.parse_args()
    r=prepare(a.source,a.target);print(json.dumps(dict(prepared=True,sourceLocalIteration=r['sourceSavedLocalIteration'],completedUpdatesBeforeRun=r['completedUpdatesBeforeRun'])))
