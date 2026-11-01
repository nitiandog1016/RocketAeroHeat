"""归档MMG本轮小型证据及准确执行源，保留真实大网格和本地运行时。"""
from pathlib import Path
import json,shutil
from prepareCGALRemeshingRuntime import ROOT,sha
from prepareMMGLocalRuntime import DEST

def main():
    out=ROOT/'docs/validation/v138';actual=ROOT/'outputs/v138_validation/actual';small=out/'hollow_test02'
    rr=json.loads((actual/'mesh_process.json').read_text());sr=json.loads((small/'report.json').read_text())
    assert rr['exitCode']==1 and rr['sourceStable'] and sr['passed']
    for n,h in rr['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in rr['outputSHA256'].items():assert sha(actual/n)==h
    for n,h in sr['filesSHA256'].items():assert sha(small/n)==h
    # 用Python整数/显式邻接独立核对247条边，避免Windows默认int32乘法编码溢出。
    import numpy as np
    from scipy.io import loadmat
    seed=json.loads((actual/'interior_seed_review.json').read_text());d=loadmat(str(actual/'hybrid_volume.mat'));T=d['tetrahedra'].astype(np.int64)-1
    for edge in seed['edges']:
        a,b=[v-1 for v in edge['vertexIDsOneBased']];rows=[v-1 for v in edge['incidentTetRowsOneBased']]
        expected=np.flatnonzero(np.any(T==a,axis=1)&np.any(T==b,axis=1));assert np.array_equal(expected,rows)
        adjacency={}
        for row in rows:
            c,e=[int(v) for v in T[row] if v not in [a,b]]
            adjacency.setdefault(c,set()).add(e);adjacency.setdefault(e,set()).add(c)
        assert all(len(v)==2 for v in adjacency.values())
        seen=set();stack=[next(iter(adjacency))]
        while stack:
            v=stack.pop()
            if v in seen:continue
            seen.add(v);stack.extend(adjacency[v]-seen)
        assert seen==set(adjacency)
    for name in ['candidate_report.json','geometry_audit.json','mesh_report.json','mesh_process.json','independent_review.json',
                 'independent_review02.json','interior_seed_review.json','failure.log',
                 'native/interface.json','native/process/process.json','native/process/stdout.log']:
        target=out/'actual'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(actual/name),str(target))
    versions=out/'source_versions/final';versions.mkdir(parents=True,exist_ok=True)
    sources=['scripts/prepareMMGLocalRuntime.py','scripts/buildMMGLocalRuntime.py','scripts/repairMMGLocalWrapper.py',
        'scripts/mmgLocalInterface.py','scripts/remeshSU2WithMMG.py','scripts/reviewMMGLocalCandidate.py','scripts/reviewMMGInteriorSeeds.py','scripts/archiveMMGV138Evidence.py',
        'tests/testMMGLocalRemeshing.py','third_party/mmg/local_tetra_remesh.cpp','scripts/buildCGALRemeshingRuntime.py',
        'scripts/prepareCGALRemeshingRuntime.py','scripts/cgalTetrahedralContract.py','tests/testCGALLocalRemeshing.py',
        'scripts/reviewSU2ExplicitQuality.py','scripts/optimizeSU2FanExterior.py','scripts/meshSU2HybridVolume.py',
        'scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py']
    for n in sources:
        target=versions/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/n),str(target))
    for name in ['COPYING','COPYING.LESSER','CMakeLists.txt','src/mmg3d/API_functions_3d.c','src/mmg3d/mmg3d1.c',
                 'src/mmg3d/libmmg3d.c','src/mmg3d/isosiz_3d.c','src/mmg3d/libmmg3d.h']:
        target=out/'official'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(DEST/'source'/name),str(target))
    runtime=json.loads((DEST/'runtime.json').read_text());assert sha(Path(runtime['executable']))==runtime['executableSHA256']
    old=json.loads((out/'runtime/runtime_before_interface_fix.json').read_text());assert sha(DEST/'versions/first/local_tetra_remesh.exe')==old['executableSHA256']
    failure=dict(preconfigureFailure='external D:/matlab Ninja path could not be made relative to E:/RocketAeroHeat',
        preconfigureExitCode=1,preconfigureSource='runtime/build_before_path_fix.py',
        firstHollowExitCode=3,firstHollowSeconds=.391,firstHollowReason='renum parameter guarded by USE_SCOTCH; library built without SCOTCH',
        firstHollowNumericalRemeshingStarted=False,firstHollowSources='source_versions/before_interface_fix',
        priorRuntime='runtime/runtime_before_interface_fix.json',priorBinaryLocal=str(DEST/'versions/first/local_tetra_remesh.exe'),
        priorBinarySHA256=old['executableSHA256'],correctedWrapperOnly=True,upstreamNumericalKernelModified=False,
        originalLibraryBuildRepeated=False,actualRocketCandidateCount=1,
        readonlySeedReviewFailure='int32 edge key overflow caught by closed-link assertion; no numeric mesh changed',
        correctedConstraintReport='actual/independent_review02.json',supersededConstraintReport='actual/independent_review.json',
        correctedEdgeReview='actual/interior_seed_review.json',independentExactIntegerConnectedLinkCheck=True)
    (out/'failures_and_correction.json').write_text(json.dumps(failure,indent=2)+'\n')
    outputFiles={p.relative_to(actual).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in actual.rglob('*') if p.is_file()}
    runtimeFiles=[p for p in DEST.rglob('*') if p.is_file()]
    summary=dict(passed=True,currentSourceSHA256={n:sha(ROOT/n) for n in sources},
        outputsFileLengthBytes=sum(r['bytes'] for r in outputFiles.values()),outputFiles=outputFiles,
        runtimeFileLengthBytes=sum(p.stat().st_size for p in runtimeFiles),runtimeFileCount=len(runtimeFiles),
        actualGeometryAccepted=False,physicsQualified=False,newQualifiedArea=0,noNewSU2OrMATLAB=True,
        independentConnectedInteriorEdgeStars=len(seed['edges']),
        runtimeAndLargeMATUploaded=False,deletedUniqueDataFiles=0,
        scope='file lengths only, not allocated disk; no whole historical outputs rescan')
    assert summary['outputsFileLengthBytes']<100*1024**2 and summary['runtimeFileLengthBytes']<350*1024**2
    (out/'storage_and_sources.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['outputFiles','currentSourceSHA256']},indent=2))
if __name__=='__main__':main()
