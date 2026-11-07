"""补充完整域必要载荷下界与真实调用分派，不把部分布局当总内存。"""
from pathlib import Path
import json,re,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import PREFIX,sha
from eilmerWorkspaceRuntime import audit_workspace_execution
start=time.monotonic();base=Path(__file__).resolve().parent;out=base/'capacity01';out.mkdir(exist_ok=False)
scope=json.loads((base/'scope03/report.json').read_text());cap=scope['capacity']
cfgpath=ROOT/'outputs/v186_validation/real01/lmrsim/config';cfg=json.loads(cfgpath.read_text())
block=cfg['block_0'];assert all(block['boundary_'+str(i)]['ghost_cell_data_available'] for i in range(2))
assert block['ncells']==705490 and block['nfaces']==1741600
n=cap['cellsIncludingGhostCapacity'];nf=cap['faces'];sizes=cap['sizes']
extra=dict(faceFlowAndGradients=nf*(sizes['FlowState']+sizes['FlowGradients']+24+sizes['WLSQGradWorkspace']),
    faceNamedCoreArrays=nf*(8+24+8+24+24+24+24+1+1+7*8+16+8+8),
    perCellStageSliceAndUdfPayload=n*(2*4*16+7*8),cellSingleGridTimeGeometry=n*(24+8+8+24))
total=cap['necessaryPartialSumBytes']+sum(extra.values())
assert total<8*1024**3
source=PREFIX/'source/src/lmr'
names=['commands/runsim.d','timemarching.d','simcore_gasdynamic_step.d','fluidfvcell.d','fluidblock.d','ufluidblock.d','lsqinterp.d','flowgradients.d','newtonkrylovsolver.d']
hashes={str(source/name):sha(source/name) for name in names};chains=[]
for name in ['commands/runsim.d','timemarching.d']:
    lines=(source/name).read_text().splitlines()
    for i,line in enumerate(lines):
        if re.search(r'\b(initTimeMarchingSimulation|integrateInTime|gasdynamic_explicit_increment_with_fixed_grid|performNewtonKrylovUpdates)\s*\(',line):
            chains.append(dict(path=name,line=i+1,text='\n'.join('%d: %s'%(j+1,lines[j]) for j in range(max(0,i-6),min(len(lines),i+8)))))
assert any('gasdynamic_explicit_increment_with_fixed_grid' in item['text'] for item in chains)
(out/'dispatch.json').write_text(json.dumps(chains,indent=2)+'\n')
result=dict(passed=True,sourceSHA256=hashes,configSHA256=sha(cfgpath),allBoundaryGhostsConstructed=True,
    additionalNecessaryPayloadBytes=extra,partialNecessaryLowerBoundBytes=total,limitBytes=8*1024**3,
    knownRemovalPayloadBytes=cap['removablePayloadSumBytes'],fullProcessFitProven=False,
    excluded='vertex/connectivity and cloud allocations, extra gas-state buffers, array capacity slack, temporary buffers, allocator/GC, runtime and address-space mappings',
    dispatchScope='fixed transient real explicit run only; check-jacobian/adjoint/steady/dual-time not admitted by this package interface',
    elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
# 相容回执绑定真实小例与独立比较，不修改原包的pending标记或原运行源。
runpath=base/'small01_run/report.json';verify=base/'verify01/report.json';rp=ROOT/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
r=json.loads(runpath.read_text());v=json.loads(verify.read_text());runtime=json.loads(rp.read_text())
m=json.loads((base/'small01/case_manifest.json').read_text())
assert r['passed'] and v['passed'] and v['maxScaledError']==0 and len(v['rejected'])==18
assert all(v['fieldAndWallBytesIdentical'].values()) and audit_workspace_execution(m['condition'],r)
compat=dict(schema='RAH.EilmerWorkspaceCompatibility.v1',passed=True,runtimePackageSHA256=sha(rp),
    binarySHA256=runtime['filesSHA256'][runtime['binary']],runReceiptSHA256=sha(runpath),verificationSHA256=sha(verify),
    caseManifestSHA256=sha(base/'small01/case_manifest.json'),physicalQualification=False,newQualifiedAreaM2=0)
(base/'compatibility.json').write_text(json.dumps(compat,indent=2)+'\n')
print(json.dumps(result))
