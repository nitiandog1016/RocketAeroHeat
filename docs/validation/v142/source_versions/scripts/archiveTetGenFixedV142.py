"""归档TetGen固定接口的实际失败、只读复核与准确执行源码；不重算数值。"""
from pathlib import Path
import json,shutil
from prepareTetGenRuntime import ROOT,DEST,sha

def main():
    out=ROOT/'docs/validation/v142';test=json.loads((out/'hollow_test/report.json').read_text())
    review=json.loads((out/'hollow_test/independent_review.json').read_text());runtime=json.loads((DEST/'runtime.json').read_text())
    assert not test['passed'] and review['passed'] and not review['realCandidateExecuted']
    receipts=[test,review,json.loads((out/'hollow_test/native/process/process.json').read_text())]+runtime['receipts']
    for r in receipts:
        for n,h in r['sourceSHA256'].items():assert sha(Path(n))==h
        assert all(r['sourceStableAfter'].values())
    sources=['scripts/prepareTetGenRuntime.py','scripts/buildTetGenRuntime.py','scripts/tetgenFixedInterface.py',
        'scripts/reviewTetGenFixedRepair.py','scripts/archiveTetGenFixedV142.py','third_party/tetgen/fixed_tetra_repair.cpp',
        'tests/testTetGenFixedRepair.py','tests/testCGALLocalRemeshing.py','scripts/reviewSU2ExplicitQuality.py',
        'scripts/optimizeSU2FanExterior.py','scripts/buildCGALRemeshingRuntime.py','scripts/prepareCGALRemeshingRuntime.py',
        'scripts/cgalTetrahedralContract.py','scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py',
        'docs/validation/v142/plan.md']
    for n in sources:
        dst=out/'source_versions'/n;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(ROOT/n),str(dst));assert sha(dst)==sha(ROOT/n)
    # 许可原字节与有限源码行证据；官方完整源码可按固定URL/blob验证重建。
    for n in ['LICENSE','README','makefile']:
        dst=out/'runtime'/('upstream_'+n);shutil.copyfile(str(DEST/'source'/n),str(dst))
    core=DEST/'source/tetgen.cxx';lines=core.read_text().splitlines()
    excerpts=[]
    for label,lo,hi in [('reconstruct_input',23336,23360),('Y_disables_input_free_classification',23841,23865),
        ('smoothing_requires_free_vertex',27440,27453),('original_point_markers',29815,29835),
        ('PLC_recover_and_carve',32210,32250),('explicit_hole_seed',22537,22555)]:
        excerpts.append(dict(label=label,firstLine=lo,lastLine=hi,text='\n'.join(lines[lo-1:hi])))
    (out/'runtime/source_contract.json').write_text(json.dumps(dict(sourceSHA256=sha(core),version='1.5.1',excerpts=excerpts),indent=2)+'\n')
    files={p.relative_to(DEST).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in DEST.rglob('*') if p.is_file()}
    evidence={p.relative_to(out).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    size=sum(r['bytes'] for r in files.values());assert size<=30*1024**2
    smallsize=sum(p.stat().st_size for p in (out/'hollow_test').rglob('*') if p.is_file());assert smallsize<100*1024**2
    result=dict(passed=True,runtimeFiles=files,runtimeFileLengthBytes=size,
        evidenceFiles=evidence,evidenceFileLengthBytes=sum(r['bytes'] for r in evidence.values()),
        sourceSHA256={n:sha(ROOT/n) for n in sources},smallFixtureFileLengthBytes=smallsize,
        smallNativeRunCount=1,realCandidateCount=0,newActualMATBytes=0,deletedUniqueDataFiles=0,
        noNewMATLABOrSU2OrUITests=True,runtimeUploaded=False,physicsQualified=False,newQualifiedArea=0,
        scope='this phase file lengths, not allocated space; receipt itself excluded from evidence tally; no historical scan')
    (out/'storage_and_sources.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['sourceSHA256','runtimeFiles','evidenceFiles']},indent=2))

if __name__=='__main__':main()
