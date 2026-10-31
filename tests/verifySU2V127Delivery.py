"""独立核对v127方案否决、主工作台证据、原始字节和受保护存储。"""
from pathlib import Path
import hashlib,json,subprocess,re,os
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v127'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def main():
    m=read(DEST/'delivery_manifest.json')
    for row in m['files']+m['sourceBindings']+m['directEvidence']:assert sha(ROOT/row['archived'])==row['sha256']
    assert [x['exitCode'] for x in m['matlab']]==[0] and all(x['sourceStable'] for x in m['matlab'])
    assert m['nativeEngineRuns']==2 and all(x['exitCode']==0 for x in m['nativeRuns'])
    c=read(DEST/'evidence/comparison.json')
    assert not c['allScreeningGatesPassed'] and not m['allScreeningGatesPassed']
    assert [k for k,v in c['screeningGates'].items() if not v]==['densityMeanImproved03','energyMeanImproved03']
    assert .08<c['meanResidualGainsLog']['rms[Rho]']<.1 and .09<c['meanResidualGainsLog']['rms[RhoE]']<.1
    assert c['pairedFaceHeatRelativeL1Difference']>.1 and not c['defaultChanged'] and not c['automaticContinuationAllowed']
    for name in ['control','cap05']:
        a=read(DEST/'evidence'/name/'sst_state_audit.json');r=read(DEST/'evidence'/name/'engine_process.json')
        assert a['interfacePassed'] and a['checkpointLocalIteration']==149 and not a['originalResidualThresholdPassed']
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert a['fullPrimitiveEOSRelativeError']<2e-12 and a['nativePrecisionHeatGradientError']<2e-6
    assert read(DEST/'evidence/export_check.json')['passed'] and read(DEST/'evidence/manifest_control_check.json')['passed']
    ui=read(DEST/'evidence/workbench01/report.json')
    assert ui['wrongWallBoundaryRejected'] and ui['oldSSTRestartCompatible'] and ui['solverRuns']==0
    for v in ui['variants'].values():assert v['faces']==19752 and v['sectionSegments']==92 and v['actualRunControls']['ITER']==150
    extrema=read(DEST/'evidence/state_extrema.json')
    assert extrema['solverRuns']==0 and extrema['cases']['control149']['h0Above110ReferenceNodeCount']==1
    assert not m['physicsQualified'] and m['newQualifiedArea']==0 and not m['largeResultsBackedUpToGitHub']
    assert m['mainWorkbenchReadOnlyTested'] and not m['newMainSolveLaunchTested'] and not m['newNativeScreenshots']
    storage=read(DEST/'evidence/storage_receipt.json')
    for p,h in storage['originalSHA256'].items():assert sha(Path(p))==h
    for row in storage['plan']:assert os.path.samefile(row['source'],row['destination'])
    assert storage['uniqueResultsDeleted']==0 and storage['directoryLogicalBytesBefore']<1700*1024**2
    for p in [ROOT/'docs/PROJECT_GOALS.md',ROOT/'docs/SU2_SST_CFL_COMPARISON_V127.md',ROOT/'docs/reviews/TECHNICAL_ROUTE_REVIEW_2026-11-01.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link:assert (p.parent/link.split('#')[0]).exists(),link
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v127'],cwd=ROOT).decode().splitlines();assert names
    indexed=set(subprocess.check_output(['git','ls-files','--','docs/validation/v127'],cwd=ROOT).decode().splitlines())
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
        indexedEvidenceFiles=len(names),allDeliveryFilesIndexed=True,nativeTerminal=2,matlabTerminal=1,matlabExitCodes=[0],
        storageSHA256Rechecked=storage['checkedFiles'],unmergedResourcesPreserved=True,
        lowerCFLCandidateAccepted=False,mainReadOnlyTested=True,physicsQualified=False,newQualifiedArea=0,
        verifierSHA256=sha(Path(__file__).resolve()))
    target=DEST/'delivery_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
