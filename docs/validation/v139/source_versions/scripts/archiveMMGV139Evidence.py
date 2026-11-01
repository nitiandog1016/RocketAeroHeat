"""归档v139实际输入、准确执行源码和只读复核；大型网格与原生运行时仅本地。"""
from pathlib import Path
import json,shutil,sys
import numpy as np
import scipy
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    out=ROOT/'docs/validation/v139';actual=ROOT/'outputs/v139_validation/actual';runtime=ROOT/'tools/runtime/mmg/5.8.0-seeds'
    r=json.loads((actual/'mesh_process.json').read_text());s=json.loads((out/'seed_test/report.json').read_text())
    review=json.loads((actual/'independent_review.json').read_text());assert r['exitCode']==1 and r['sourceStable'] and s['passed'] and review['passed']
    for receipt in [r,review,json.loads((actual/'native/process/process.json').read_text())]:
        for n,h in receipt['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in r['outputSHA256'].items():assert sha(actual/n)==h
    for n,h in s['filesSHA256'].items():assert sha(out/'seed_test'/n)==h
    meta=json.loads((runtime/'runtime.json').read_text());assert sha(Path(meta['executable']))==meta['executableSHA256']
    for n,h in meta['sourceSHA256'].items():assert sha(Path(n))==h
    targets=['candidate_report.json','geometry_audit.json','seed_geometry_audit.json','split_report.json','mesh_report.json',
        'mesh_process.json','independent_review.json','failure.log','native/interface.json','native/process/process.json','native/process/stdout.log']
    for name in targets:
        dst=out/'actual'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(actual/name),str(dst))
    sources=['scripts/buildMMGSeedRuntime.py','scripts/splitClosedTetEdges.py','scripts/mmgLocalInterface.py',
        'scripts/remeshSU2WithSeeds.py','scripts/reviewMMGSeedCandidate.py','scripts/archiveMMGV139Evidence.py',
        'tests/testMMGSeedOptimization.py','third_party/mmg/local_tetra_remesh.cpp','scripts/buildCGALRemeshingRuntime.py',
        'scripts/prepareCGALRemeshingRuntime.py','scripts/prepareMMGLocalRuntime.py','scripts/cgalTetrahedralContract.py',
        'scripts/reviewMMGLocalCandidate.py','scripts/reviewSU2ExplicitQuality.py','scripts/optimizeSU2FanExterior.py',
        'scripts/meshSU2HybridVolume.py','scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py']
    for name in sources:
        dst=out/'source_versions'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/name),str(dst));assert sha(dst)==sha(ROOT/name)
    outputs={p.relative_to(actual).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in actual.rglob('*') if p.is_file()}
    runtimeFiles={p.relative_to(runtime).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in runtime.rglob('*') if p.is_file()}
    report=dict(passed=True,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        actualCandidateCount=1,smallNumericalTestCount=1,noNewMATLABOrSU2=True,noUIChangesOrReruns=True,
        sourceSHA256={name:sha(ROOT/name) for name in sources},outputFiles=outputs,runtimeFiles=runtimeFiles,
        outputsFileLengthBytes=sum(v['bytes'] for v in outputs.values()),runtimeFileLengthBytes=sum(v['bytes'] for v in runtimeFiles.values()),
        candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,deletedUniqueDataFiles=0,
        largeMATAndRuntimeUploaded=False,scope='this phase file lengths only, not NTFS allocation; no historical storage rescan')
    assert report['outputsFileLengthBytes']<100*1024**2 and report['runtimeFileLengthBytes']<15*1024**2
    (out/'storage_and_sources.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','outputFiles','runtimeFiles']},indent=2))

if __name__=='__main__':main()
