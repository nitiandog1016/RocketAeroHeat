"""归档底面外层重分布的已验收几何、独立复核和准确执行源，不重跑计算。"""
from pathlib import Path
import json,shutil
from prepareTetGenRuntime import ROOT,sha

def main():
    out=ROOT/'docs/validation/v144';actual=ROOT/'outputs/v144_validation/actual'
    native=actual/'native/process/process.json';mesh=json.loads((actual/'mesh_process.json').read_text())
    candidate=json.loads((actual/'candidate_report.json').read_text());review=json.loads((actual/'independent_review.json').read_text())
    small=json.loads((out/'fixture/report.json').read_text());assert mesh['exitCode']==0 and mesh['sourceStable'] and candidate['candidateAccepted'] and review['passed'] and small['passed']
    for r in [mesh,review,json.loads(native.read_text()),json.loads((out/'actual_supervisor/process.json').read_text()),
        json.loads((out/'fixture_process/process.json').read_text()),json.loads((out/'review_process/process.json').read_text())]:
        assert all(r['sourceStableAfter'].values())
        for n,h in r['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in small['sourceSHA256'].items():assert sha(Path(n))==h
    assert sha(out/'fixture/fixture.npz')==small['fixtureSHA256']
    for n,h in mesh['outputSHA256'].items():assert sha(actual/n)==h
    cg=json.loads((actual/'boundary_audit/process.json').read_text());assert cg['exitCode']==0 and cg['runtimeStableAfter']
    for n,h in cg['outputSHA256'].items():assert sha(actual/'boundary_audit'/n)==h
    for p in actual.rglob('*'):
        if p.is_file() and p.suffix in ['.json','.log']:
            dst=out/'actual'/p.relative_to(actual);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(dst))
    sources=['scripts/recenterBaseOuterColumn.py','scripts/meshSU2BaseCapEnvelope.py','scripts/reviewBaseCapEnvelope.py',
        'scripts/archiveBaseCapEnvelopeV144.py','tests/testBaseOuterColumn.py','scripts/buildSU2CreaseFan.py',
        'scripts/tetgenPLCInterface.py','scripts/tetgenFixedInterface.py','scripts/reviewTetGenFixedRepair.py',
        'scripts/optimizeSU2FanExterior.py','scripts/reviewSU2ExplicitQuality.py','scripts/cgalTetrahedralContract.py',
        'scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py','scripts/meshSU2HybridVolume.py',
        'scripts/runCGALSurfaceAudit.py','scripts/cgalSurfaceIntersections.py','scripts/buildCGALRemeshingRuntime.py',
        'docs/validation/v144/plan.md']
    for n in sources:
        dst=out/'source_versions'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/n),str(dst));assert sha(dst)==sha(ROOT/n)
    files={p.relative_to(actual).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in actual.rglob('*') if p.is_file()}
    evidence={p.relative_to(out).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    size=sum(r['bytes'] for r in files.values());assert size<100*1024**2
    result=dict(passed=True,sourceSHA256={n:sha(ROOT/n) for n in sources},outputFiles=files,outputsFileLengthBytes=size,
        evidenceFiles=evidence,evidenceFileLengthBytes=sum(r['bytes'] for r in evidence.values()),newRuntimeBytes=0,
        actualCandidateCount=1,actualCGALSurfaceCount=1,actualPLCCount=1,smallFixtureCount=1,
        geometryCandidateAccepted=True,physicsQualified=False,newQualifiedArea=0,
        noNewMATLABOrSU2OrUITests=True,largeOutputsUploaded=False,deletedUniqueDataFiles=0,
        scope='current phase file lengths only, not allocated disk; receipt itself and later prose excluded from evidence tally')
    (out/'storage_and_sources.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['sourceSHA256','outputFiles','evidenceFiles']},indent=2))

if __name__=='__main__':main()
