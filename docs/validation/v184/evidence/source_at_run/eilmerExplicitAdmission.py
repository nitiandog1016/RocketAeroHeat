"""完整域显式发展准入：准确候选目录、物理/网格、已有证据和固定资源必须闭合。"""
from pathlib import Path
import copy,hashlib,json
from eilmerExplicitProfile import DEVELOPMENT_PROFILE,validate_explicit


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def physical_hash(condition):
    c=copy.deepcopy(condition);c.pop('explicit',None);c.pop('development',None);c.pop('maxNewtonSteps',None)
    return hashlib.sha256(json.dumps(c,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def need(ok,message):
    if not ok:raise ValueError(message)


def audit_admission(manifest,case):
    """只读小证据；不重跑旧小例或遍历大场。返回实际锁定的登记文件。"""
    c=manifest['condition'];d=validate_explicit(c['explicit']);need(d['profile']==DEVELOPMENT_PROFILE,'非完整域方案')
    path=Path(d['candidate']);need(path.is_file() and sha(path)==d['candidateSHA256'],'候选登记缺失/改变')
    a=json.loads(path.read_text());need(a['schema']=='RAH.EilmerExplicitCandidate.v1','候选版本错误')
    need(case is not None and Path(a['case']).resolve()==Path(case).resolve(),'候选目录不符')
    need(a['profile']==DEVELOPMENT_PROFILE and a['maxSteps']==200 and a['checkpointStep']==100,'候选步数不符')
    need(a['physicalConditionSHA256']==physical_hash(c),'候选物理工况改变')
    mesh=manifest['mesh']
    need(mesh['sourceSHA256']==a['sourceMeshSHA256'] and mesh['cells']==705490 and mesh['nodes']==335726 and mesh['wallFaces']==19752,'候选原网格不符')
    need(abs(mesh['wallAreaM2']-1.9646292565715551)<2e-13,'原面积改变')
    expected=dict(nativeWallClockSeconds=600,runEntrySeconds=660,parentSeconds=720,prepareSeconds=120,
                  cpuCount=2,memoryMiB=8192,outputMiB=1024,automaticContinuation=0)
    need(a['budget']==expected,'候选预算变化')
    locks=[path]
    for entry in a['evidence'].values():
        p=Path(entry['path']);need(p.is_file() and sha(p)==entry['sha256'],'候选证据改变')
        locks.append(p)
    need(set(a['evidence'])=={'sourceManifest','sourcePreparation','sourceGeometry','explicitInterface','explicitRun'},'候选证据集合不完整')
    docs={k:json.loads(Path(v['path']).read_text()) for k,v in a['evidence'].items()}
    need(docs['sourcePreparation']['passed'] and docs['sourceGeometry']['passed'] and docs['sourceGeometry']['originalCellOrderAndBoundaryBijection'],'原准备/几何证据失败')
    need(docs['explicitInterface']['passed'] and docs['explicitInterface']['independentSignedHeatCSVScaledError']<=2e-12,'显式接口未过')
    need(docs['explicitRun']['passed'] and docs['explicitRun']['nativeCompletion']['finalStep']==6,'缺少六步原生证据')
    need(docs['sourceManifest']['mesh']['sourceSHA256']==a['sourceMeshSHA256'],'来源清单网格不符')
    need(docs['sourcePreparation']['caseManifestSHA256']==a['evidence']['sourceManifest']['sha256'],'原准备清单绑定错误')
    return a,locks
