"""固定原壁面/外域时沿原法向网格列凸组合守恒变量；不宣称全域积分守恒。"""
from pathlib import Path
import argparse,hashlib,json,re
import numpy as np
from scipy.io import loadmat

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def transfer(oldFolder,newFolder,source,output):
    assert not output.exists()
    # 源迭代来自真实终态、具名文件及已锁定谱系，不能写死为某次研究步数。
    receiptPath=source.parent/'engine_process.json'
    receipt=json.loads(receiptPath.read_text(encoding='utf-8-sig'))
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut'] and not receipt['cancelled']
    assert all(receipt['inputsStable'].values()) and receipt['runtimeFlavor']=='official'
    match=re.fullmatch(r'restart_flow_(\d{6,})\.csv',source.name);assert match
    step=int(match[1]);assert sha(source)==receipt['outputSHA256'][source.name]
    assert sha(source)==receipt['outputSHA256']['restart_flow.csv'], '仅支持源任务最终具名场'
    lineagePath=source.parent/'continuation.json';before=0;locks=[receiptPath,source]
    if lineagePath.exists():
        assert sha(lineagePath)==receipt['inputSHA256'][str(lineagePath)]
        lineage=json.loads(lineagePath.read_text(encoding='utf-8-sig'))
        assert lineage['sameMesh'] and lineage['samePhysicsAndNumerics']
        before=int(lineage['completedUpdatesBeforeRun']);locks.append(lineagePath)
    for folder in (oldFolder,newFolder):
        rp=folder/'mesh_process.json';mr=json.loads(rp.read_text(encoding='utf-8-sig'))
        assert mr['status']=='terminal' and mr['exitCode']==0 and mr['sourceStable'] and mr['inputStable']
        for name in ('hybrid_volume.mat','mesh_report.json'):
            assert sha(folder/name)==mr['outputSHA256'][name];locks.append(folder/name)
        locks.append(rp)
    a=loadmat(str(oldFolder/'hybrid_volume.mat'));b=loadmat(str(newFolder/'hybrid_volume.mat'))
    oldReceipt=json.loads((oldFolder/'mesh_process.json').read_text());inp=loadmat(oldReceipt['input'])
    oldH=np.asarray(json.loads(str(inp['optionsJSON'][0]))['layerHeights'])
    newReport=json.loads((newFolder/'mesh_report.json').read_text());newH=np.asarray(newReport['normalHeightsM'])
    V=a['sourceVertices'];n=len(V);oldN=len(oldH)-1;newN=len(newH)-1
    assert np.array_equal(V,b['sourceVertices']) and np.array_equal(a['sourceFaces'],b['sourceFaces']) and oldH[-1]==newH[-1]
    state=np.loadtxt(str(source),delimiter=',',skiprows=1)
    assert state.shape==(len(a['vertices']),9) and np.array_equal(state[:,0],np.arange(len(state)))
    # 原生ASCII约15位有效数字；实测最大舍入差8.88e-16m，容差固定1e-14m。
    coordinateError=float(np.max(abs(state[:,1:4]-a['vertices'])))
    assert coordinateError<=1e-14, '源节点编号/坐标不匹配'
    result=np.empty((len(b['vertices']),5));brackets=[]
    for j,h in enumerate(newH):
        k=min(np.searchsorted(oldH,h,side='right')-1,oldN-1);k=max(k,0)
        weight=(h-oldH[k])/(oldH[k+1]-oldH[k]);assert 0<=weight<=1
        lo=np.arange(k*n,(k+1)*n);hi=lo+n;ids=np.arange(j*n,(j+1)*n)
        assert np.max(abs((1-weight)*a['vertices'][lo]+weight*a['vertices'][hi]-b['vertices'][ids]))<2e-14
        result[ids]=(1-weight)*state[lo,4:9]+weight*state[hi,4:9]
        brackets.append(dict(newLayer=j,oldLowerLayer=int(k),weight=float(weight)))
    result[(newN+1)*n:]=state[(oldN+1)*n:,4:9]
    assert np.array_equal(result[:n],state[:n,4:9])
    assert np.array_equal(result[newN*n:],state[oldN*n:,4:9])
    rho=result[:,0];e=result[:,4]/rho-.5*np.sum((result[:,1:4]/rho[:,None])**2,axis=1)
    assert np.isfinite(result).all() and min(rho)>0 and min(e)>0
    values=np.column_stack([np.arange(len(result)),b['vertices'],result])
    np.savetxt(str(output),values,delimiter=',',fmt=['%d']+['%.17e']*8,comments='',
               header='"PointID","x","y","z","Density","Momentum_x","Momentum_y","Momentum_z","Energy"')
    report=dict(kind='normal-column-state-transfer',sourceState=str(source),sourceSavedLocalIteration=step,
                sourceCumulativeIteration=before+step,sameMesh=False,newMeshLocalIterationStartsAt=0,
                coordinateToleranceM=1e-14,sourceCoordinateRoundingErrorM=coordinateError,
                oldLayers=oldN,newLayers=newN,newNodes=len(result),exactWallAndExterior=True,
                brackets=brackets,minimumDensity=float(min(rho)),minimumInternalEnergy=float(min(e)),
                strictlyConservativeRemap=False,initialGuessOnly=True,physicsQualified=False,
                sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve()]+locks},
                outputSHA256=sha(output))
    output.with_suffix('.transfer.json').write_text(json.dumps(report,indent=2)+'\n');return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old-folder',type=Path,required=True);p.add_argument('--new-folder',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=transfer(a.old_folder.resolve(),a.new_folder.resolve(),a.source.resolve(),a.output.resolve());print(json.dumps({k:r[k] for k in ['newNodes','minimumDensity','minimumInternalEnergy']}))
