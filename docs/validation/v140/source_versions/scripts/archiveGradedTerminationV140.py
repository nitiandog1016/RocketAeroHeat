"""归档v140准确执行源、小例和预检拒绝；保留独有局部网格而不复制完整工程。"""
from pathlib import Path
import json,shutil
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    out=ROOT/'docs/validation/v140';actual=ROOT/'outputs/v140_validation/actual'
    r=json.loads((actual/'mesh_process.json').read_text());v=json.loads((actual/'independent_review.json').read_text())
    s=json.loads((out/'fixture/report.json').read_text());assert r['exitCode']==1 and not r['gmshStarted'] and r['sourceStable'] and v['passed'] and s['passed']
    for receipt in [r,v,json.loads((out/'fixture_process/process.json').read_text())]:
        for n,h in receipt['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in r['outputSHA256'].items():assert sha(actual/n)==h
    assert sha(out/'fixture/fixture.npz')==s['fixtureSHA256']
    for name in ['mesh_process.json','preparation_report.json','shell_geometry_audit.json','independent_review.json','failure.log']:
        dst=out/'actual'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(actual/name),str(dst))
    sources=['scripts/terminatePrismLayers.py','scripts/meshSU2GradedTermination.py','scripts/reviewGradedTermination.py',
        'scripts/archiveGradedTerminationV140.py','tests/testGradedPrismTermination.py','scripts/meshSU2HybridVolume.py',
        'scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py','scripts/reviewSU2ExplicitQuality.py',
        'scripts/optimizeSU2FanExterior.py','scripts/runCGALSurfaceAudit.py','scripts/cgalSurfaceIntersections.py',
        'scripts/buildCGALRemeshingRuntime.py','scripts/prepareCGALRemeshingRuntime.py']
    for name in sources:
        dst=out/'source_versions'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/name),str(dst));assert sha(dst)==sha(ROOT/name)
    files={p.relative_to(actual).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in actual.rglob('*') if p.is_file()}
    report=dict(passed=True,sourceSHA256={n:sha(ROOT/n) for n in sources},outputFiles=files,
        outputsFileLengthBytes=sum(v['bytes'] for v in files.values()),newRuntimeBytes=0,deletedUniqueDataFiles=0,
        smallCaseCount=1,actualPreparationCount=1,actualExteriorMeshingCount=0,noNewMATLABOrSU2OrUITests=True,
        noNewCGALSurfaceAudit=True,candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,
        unexecutedBranch='CGAL and Gmsh branch after rejected quality gate; no claim of validated complete pipeline',
        lockScope='small process locks entry and template; shared audit dependencies hashed at delivery only; actual process locks listed inputs and audit sources',
        largeMATUploaded=False,scope='this phase file lengths only, not allocated disk; no historical storage rescan')
    assert report['outputsFileLengthBytes']<100*1024**2
    (out/'storage_and_sources.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','outputFiles']},indent=2))

if __name__=='__main__':main()
