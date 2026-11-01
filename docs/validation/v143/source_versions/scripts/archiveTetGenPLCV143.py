"""保存PLC真实拒绝、独立复核与准确执行源，不重跑数值或改历史输出。"""
from pathlib import Path
import json,shutil
from prepareTetGenRuntime import ROOT,sha
from buildTetGenPLCRuntime import PLC

def main():
    out=ROOT/'docs/validation/v143';actual=ROOT/'outputs/v143_validation/actual'
    receipt=json.loads((actual/'mesh_process.json').read_text());review=json.loads((actual/'independent_review.json').read_text())
    small=json.loads((out/'hollow_test/report.json').read_text());assert receipt['exitCode']==1 and receipt['sourceStable'] and small['passed'] and review['passed']
    receipts=[receipt,review,small,json.loads((actual/'rear_interface.json').read_text()),
        json.loads((actual/'native/process/process.json').read_text()),json.loads((out/'hollow_test/native/process/process.json').read_text()),
        json.loads((out/'build/process.json').read_text()),json.loads((out/'actual_supervisor/process.json').read_text()),
        json.loads((out/'review_process/process.json').read_text())]
    for r in receipts:
        assert all(r['sourceStableAfter'].values())
        for n,h in r['sourceSHA256'].items():assert sha(Path(n))==h
    for n,h in receipt['outputSHA256'].items():assert sha(actual/n)==h
    for p in actual.rglob('*'):
        if p.is_file() and p.suffix in ['.json','.log']:
            dst=out/'actual'/p.relative_to(actual);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(dst))
    sources=['scripts/buildTetGenPLCRuntime.py','scripts/tetgenPLCInterface.py','scripts/meshSU2TetGenPLC.py',
        'scripts/reviewTetGenPLC.py','scripts/locateTetGenRearInterface.py','scripts/archiveTetGenPLCV143.py',
        'third_party/tetgen/fixed_surface_plc.cpp','tests/testTetGenPLC.py','scripts/tetgenFixedInterface.py',
        'scripts/reviewTetGenFixedRepair.py','scripts/optimizeSU2FanExterior.py','scripts/cgalTetrahedralContract.py',
        'scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py','scripts/meshSU2HybridVolume.py',
        'scripts/reviewSU2ExplicitQuality.py','scripts/buildCGALRemeshingRuntime.py','scripts/prepareTetGenRuntime.py',
        'docs/validation/v143/plan.md']
    for n in sources:
        dst=out/'source_versions'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/n),str(dst));assert sha(dst)==sha(ROOT/n)
    files={p.relative_to(actual).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in actual.rglob('*') if p.is_file()}
    runtime={p.relative_to(PLC).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in PLC.rglob('*') if p.is_file()}
    evidence={p.relative_to(out).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    total=sum(r['bytes'] for r in files.values());rsize=sum(r['bytes'] for r in runtime.values());assert total<100*1024**2 and rsize<10*1024**2
    report=dict(passed=True,sourceSHA256={n:sha(ROOT/n) for n in sources},outputFiles=files,outputsFileLengthBytes=total,
        runtimeFiles=runtime,newRuntimeFileLengthBytes=rsize,evidenceFiles=evidence,
        evidenceFileLengthBytes=sum(r['bytes'] for r in evidence.values()),smallNativeRunCount=1,actualNativeRunCount=1,
        newMATLABOrSU2OrUITests=False,candidateAccepted=False,physicsQualified=False,newQualifiedArea=0,
        deletedUniqueDataFiles=0,largeOutputsAndRuntimeUploaded=False,
        scope='current phase file lengths only, not allocated disk; receipt itself and later prose excluded from evidence tally')
    (out/'storage_and_sources.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','outputFiles','runtimeFiles','evidenceFiles']},indent=2))

if __name__=='__main__':main()
