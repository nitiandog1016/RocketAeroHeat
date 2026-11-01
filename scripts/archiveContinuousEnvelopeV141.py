"""归档连续包络的实际失败与独立审核，保存准确执行源和小型证据。"""
from pathlib import Path
import json,shutil
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    out=ROOT/'docs/validation/v141';actual=ROOT/'outputs/v141_validation/actual'
    r=json.loads((actual/'mesh_process.json').read_text());review=json.loads((actual/'independent_review.json').read_text())
    small=json.loads((out/'fixture/report.json').read_text());assert r['exitCode']==1 and r['sourceStable'] and review['passed'] and small['passed']
    for receipt in [r,review,json.loads((out/'fixture_process/process.json').read_text()),json.loads((actual/'exterior_process/process.json').read_text())]:
        for n,h in receipt['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in r['outputSHA256'].items():assert sha(actual/n)==h
    assert sha(out/'fixture/fixture.npz')==small['fixtureSHA256']
    cg=json.loads((actual/'boundary_audit/process.json').read_text());assert cg['exitCode']==0 and cg['runtimeStableAfter'] and all(cg['sourceStableAfter'].values())
    for n,h in cg['outputSHA256'].items():assert sha(actual/'boundary_audit'/n)==h
    names=['mesh_process.json','candidate_report.json','geometry_audit.json','preparation_report.json','shell_geometry_audit.json',
        'independent_review.json','failure.log','boundary_audit/process.json','boundary_audit/intersections.json','boundary_audit/stdout.log',
        'exterior_process/process.json','exterior_process/stdout.log']
    for n in names:
        dst=out/'actual'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(actual/n),str(dst))
    sources=['scripts/redistributePrismEnvelope.py','scripts/meshSU2ContinuousEnvelope.py','scripts/reviewContinuousEnvelope.py',
        'scripts/archiveContinuousEnvelopeV141.py','tests/testContinuousPrismEnvelope.py','scripts/terminatePrismLayers.py',
        'scripts/meshSU2HybridVolume.py','scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py',
        'scripts/reviewSU2ExplicitQuality.py','scripts/optimizeSU2FanExterior.py','scripts/mmgLocalInterface.py','scripts/cgalTetrahedralContract.py',
        'scripts/runCGALSurfaceAudit.py','scripts/cgalSurfaceIntersections.py','scripts/buildCGALRemeshingRuntime.py','scripts/prepareCGALRemeshingRuntime.py']
    for n in sources:
        dst=out/'source_versions'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/n),str(dst));assert sha(dst)==sha(ROOT/n)
    files={p.relative_to(actual).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in actual.rglob('*') if p.is_file()}
    candidate=json.loads((actual/'candidate_report.json').read_text());cut=json.loads((actual/'preparation_report.json').read_text())['cutSourceXM']
    centers=candidate['badTetCentersM'];summary=dict(passed=True,sourceSHA256={n:sha(ROOT/n) for n in sources},outputFiles=files,
        outputsFileLengthBytes=sum(v['bytes'] for v in files.values()),newRuntimeBytes=0,deletedUniqueDataFiles=0,
        smallCaseCount=1,actualPreparationCount=1,actualCGALSurfaceCount=1,actualExteriorMeshingCount=1,
        badTetCentersAtOrAboveNoseCut=sum(p[0]>=cut for p in centers),badTetCentersBelowNoseCut=sum(p[0]<cut for p in centers),
        noNewMATLABOrSU2OrUITests=True,TetGenInstalledOrTested=False,candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,
        largeMATUploaded=False,scope='this phase file lengths only, not allocated disk; no historical storage rescan')
    assert summary['outputsFileLengthBytes']<100*1024**2
    (out/'storage_and_sources.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['sourceSHA256','outputFiles']},indent=2))

if __name__=='__main__':main()
