"""有限核对显式模式仍分配的NK工作区，生成未构建草稿；不改共享源码或运行包。"""
from pathlib import Path
import difflib
import json
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import PREFIX,sha

start=time.monotonic();base=Path(__file__).resolve().parent;out=base/'workspace01';out.mkdir(exist_ok=False)
original=PREFIX/'source/src/lmr';pkg=ROOT/'tools/runtime/eilmer-explicit-dense-v1'
names=r'\b(?:saved_gradients|saved_lsqgradients|saved_source_terms|cell_jacobians|grad_save|gradients_save|Q_save|invBlockDiag|dRdU)\b'
accesses=[]
for path in sorted(original.rglob('*.d')):
    if not path.is_file():continue
    for line,text in enumerate(path.read_text().splitlines(),1):
        if re.search(names,text):accesses.append(dict(path=str(path.relative_to(original)),line=line,text=text))
assert {x['path'] for x in accesses}<={'coredata.d','fluidblock.d','fluidfvcell.d','ufluidblock.d','newtonkrylovsolver.d'}
(out/'workspace_accesses.json').write_text(json.dumps(accesses,indent=2)+'\n')
spans={'fluidblock.d':[(289,301),(1480,1518),(1565,1593),(1785,1809)],
       'fluidfvcell.d':[(237,255),(1114,1133)],'ufluidblock.d':[(327,339)],
       'coredata.d':[(42,55)],'newtonkrylovsolver.d':[(2616,2651)]}
for name,ranges in spans.items():
    lines=(original/name).read_text().splitlines()
    (out/(name+'.txt')).write_text('\n'.join(f'{i}: {lines[i-1]}' for lo,hi in ranges for i in range(lo,hi+1))+'\n')
assert 'fs.mu_t = fmin(fs.mu_t, factor * fs.gas.mu);' in (original/'fluidfvcell.d').read_text()
patches=[]
for name,path,old in [
    ('fluidblock.d',pkg/'overlay/src/lmr/fluidblock.d','''        version(newton_krylov) {
            celldata.cell_jacobians.length = (ncells+nghost)*neq*neq;
            celldata.saved_gradients.reserve(ncells + nghost);
            foreach (n; 0 .. ncells+nghost) celldata.saved_gradients ~= FlowGradients(myConfig);
            celldata.saved_source_terms.length = (ncells + nghost)*neq;
            celldata.doNotPerturb.length = (ncells+nghost);
        }'''),
    ('ufluidblock.d',original/'ufluidblock.d','''        celldata.saved_lsqgradients.reserve(ncells + nghost);
        foreach (n; 0 .. ncells+nghost) celldata.saved_lsqgradients ~= LSQInterpGradients(myConfig.n_species, myConfig.n_modes, myConfig.turb_model.nturb);''')]:
    before=path.read_text();assert before.count(old)==1
    if name=='fluidblock.d':
        body=old.splitlines()[1:-1]
        new='        version(newton_krylov) {\n            // 保留NK编译物理分派；显式流动不分配其矩阵/备份工作区。\n            if (myConfig.solverMode != SolverMode.transient) {\n'+'\n'.join('    '+line for line in body)+'\n            }\n        }'
    else:
        new='        // 显式路径不使用Jacobian/伴随备份梯度，其他模式保留原分配。\n        if (myConfig.solverMode != SolverMode.transient) {\n'+'\n'.join('    '+line for line in old.splitlines())+'\n        }'
    after=before.replace(old,new)
    patches.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='a/src/lmr/'+name,tofile='b/src/lmr/'+name))
patch=ROOT/'patches/eilmer/explicit-workspace-v1.patch';assert not patch.exists();patch.write_text(''.join(patches))
gcroot=PREFIX/'ldc2-1.41.0-linux-x86_64/import'
gcconfig=gcroot/'core/gc/config.d';gcimpl=gcroot/'core/internal/gc/impl/conservative/gc.d'
for path,lo,hi in [(gcconfig,17,35),(gcimpl,2114,2141)]:
    lines=path.read_text().splitlines()
    (out/('druntime-'+path.name+'.txt')).write_text('\n'.join(f'{i}: {lines[i-1]}' for i in range(lo,hi+1))+'\n')
n=705490+26324
report=dict(passed=True,status='unbuilt-workspace-only-draft',sourceCommit='f53f4609a0331d48efee69a4e4f3c3598378cc03',
    packageSHA256=sha(pkg/'runtime.json'),sourceSHA256={str(original/name):sha(original/name) for name in spans},
    runtimeSourceSHA256={str(p):sha(p) for p in [gcconfig,gcimpl]},directWorkspaceAccesses=len(accesses),
    denseJacobianPayloadBytes=n*7*7*8,savedSourcePayloadBytes=n*7*8,
    minimumNamedDensePayloadBytes=n*(49+7)*8,
    additionalSavedGradientsAndPerCellMatricesNotQuantified=True,
    NKCompileFlagMustRemain=True,reason='NK=0 changes turbulence_viscosity_limit and is not a storage-only build option.',
    GCHeapSizeFactorNotProvenCause=True,reasonGC='heapSizeFactor sets collection thresholds; RSS alone does not establish the failed request or mapped virtual capacity.',
    patch=str(patch.relative_to(ROOT)),patchSHA256=sha(patch),engineBuilds=0,nativeRuns=0,
    fullDynamicCoveragePending=True,failedAllocationLineKnown=False,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
