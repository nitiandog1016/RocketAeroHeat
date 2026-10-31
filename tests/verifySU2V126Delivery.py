"""独立核对v126交付证据、真实终态、限制、存储与Git暂存字节。"""
from pathlib import Path
import hashlib,json,subprocess,re,os
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v126'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def main():
    m=read(DEST/'delivery_manifest.json')
    for row in m['files']+m['sourceBindings']+m['directEvidence']:assert sha(ROOT/row['archived'])==row['sha256']
    assert [x['exitCode'] for x in m['matlab']]==[1,0] and all(x['sourceStable'] for x in m['matlab'])
    assert m['nativeEngineRuns']==3 and all(x['exitCode']==0 for x in m['nativeRuns'])
    assert m['fullDomainGeneralSSTInterfacePassed'] and m['mainFiveStepStartTested'] and m['nativeHighPrecisionTemperatureWritten']
    a=read(DEST/'evidence/actual_surface_audit.json');p=read(DEST/'native_main/workflow/pipeline_process.json')
    assert a['interfacePassed'] and not a['diagnosticOnly'] and a['checkpointLocalIteration']==49
    assert a['fullPrimitiveEOSRelativeError']<2e-12 and a['nativePrecisionHeatGradientError']<2e-6
    assert not a['originalResidualThresholdPassed'] and a['newQualifiedArea']==0
    assert p['status']=='terminal' and p['exitCode']==0 and not p['cancelled'] and not p['timedOut']
    assert read(DEST/'evidence/export_check.json')['passed']
    ui=read(DEST/'evidence/main02/report.json');assert ui['sameConditionReopened'] and ui['webgl']['webglError']==0
    assert ui['originalFaces']==19752 and ui['sectionSegments']==92 and not ui['physicsQualified']
    d=read(DEST/'evidence/multistep_diagnosis.json')
    assert d['energyResidualReboundLog']>.9 and not d['actualCFLHistoryAvailable'] and d['solverRuns']==0
    assert not m['physicsQualified'] and m['newQualifiedArea']==0 and not m['largeResultsBackedUpToGitHub']
    storage=read(DEST/'evidence/storage_receipt.json')
    for p,h in storage['originalSHA256'].items():assert sha(Path(p))==h
    for row in storage['plan']:assert os.path.samefile(row['source'],row['destination'])
    assert storage['uniqueResultsDeleted']==0 and storage['directoryLogicalBytesBefore']<1200*1024**2
    for p in [ROOT/'docs/PROJECT_GOALS.md',ROOT/'docs/SU2_SST_MULTISTEP_V126.md',ROOT/'docs/reviews/TECHNICAL_ROUTE_REVIEW_2026-11-01.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link:assert (p.parent/link.split('#')[0]).exists(),link
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v126'],cwd=ROOT).decode().splitlines();assert names
    indexed=set(subprocess.check_output(['git','ls-files','--','docs/validation/v126'],cwd=ROOT).decode().splitlines())
    allfiles={p.relative_to(ROOT).as_posix() for p in DEST.rglob('*') if p.is_file()};assert allfiles<=indexed,allfiles-indexed
    raw=subprocess.check_output(['git','cat-file','--batch'],cwd=ROOT,input=''.join(':'+p+'\n' for p in names).encode());cursor=0
    for name in names:
        end=raw.index(b'\n',cursor);header=raw[cursor:end].decode().split();assert header[1]=='blob';n=int(header[2])
        assert hashlib.sha256(raw[end+1:end+1+n]).hexdigest()==sha(ROOT/name),name;cursor=end+1+n+1
    assert cursor==len(raw)
    for item in m['unmergedWorkbenchResources']:
        assert not item['committedByThisVersion'] and sha(ROOT/item['path'])==item['taskSHA256']
    staged=subprocess.check_output(['git','diff','--cached','--name-only'],cwd=ROOT).decode().splitlines()
    assert not any(p in staged for p in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py'])
    assert not any(p.startswith('outputs/') for p in staged)
    report=dict(passed=True,evidenceFiles=len(m['files']),sourceBindings=len(m['sourceBindings']),
        indexedEvidenceFiles=len(names),allDeliveryFilesIndexed=True,nativeTerminal=3,matlabTerminal=2,matlabExitCodes=[1,0],
        storageSHA256Rechecked=storage['checkedFiles'],unmergedResourcesPreserved=True,
        generalSSTInterfacePassed=True,mainFiveStepStartTested=True,physicsQualified=False,newQualifiedArea=0,
        verifierSHA256=sha(Path(__file__).resolve()))
    target=DEST/'delivery_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
