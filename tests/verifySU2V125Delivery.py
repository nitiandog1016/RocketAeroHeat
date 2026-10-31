"""独立核对v125证据字节/终态/导出/存储与Git暂存，保留原失败和范围限制。"""
from pathlib import Path
import hashlib,json,subprocess,re,os
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v125'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))

def main():
    m=read(DEST/'delivery_manifest.json')
    for row in m['files']+m['sourceBindings']:assert sha(ROOT/row['archived'])==row['sha256']
    assert [x['exitCode'] for x in m['matlab']]==[1,0] and all(x['sourceStable'] for x in m['matlab'])
    assert m['nativeEngineRuns']==1 and m['nativeEngineExitCodes']==[0]
    assert m['phaseDiagnosticPassed'] and not m['nativeHighPrecisionTemperatureWritten']
    assert not m['fullDomainGeneralSSTInterfacePassed'] and not m['mainFiveStepStartTested']
    d=read(DEST/'evidence/phase_diagnosis.json');a=read(DEST/'evidence/phase01/sst_state_audit.json')
    assert d['phaseProvedForThisUpdate'] and d['inputKHeatGradientError']<2e-6<d['updatedKHeatGradientError']
    assert a['interfacePassed'] and a['diagnosticOnly'] and a['checkpointLocalIteration']==0
    assert a['last25HeatPowerRelativeChange'] is None and not a['originalResidualThresholdPassed']
    assert read(DEST/'evidence/export_check.json')['passed']
    assert read(DEST/'evidence/view02/report.json')['oldFailedFieldStillRejected']
    assert not m['physicsQualified'] and m['newQualifiedArea']==0
    storage=read(DEST/'evidence/storage_receipt.json')
    for p,h in storage['originalSHA256'].items():assert sha(Path(p))==h
    for row in storage['plan']:assert os.path.samefile(row['source'],row['destination'])
    assert storage['uniqueResultsDeleted']==0 and storage['distinctFileLengthBefore']<600*1024**2
    for p in [ROOT/'docs/PROJECT_GOALS.md',ROOT/'docs/SU2_SST_PHASE_V125.md']:
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link:assert (p.parent/link.split('#')[0]).exists(),link
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v125'],cwd=ROOT).decode().splitlines();assert names
    indexed=set(subprocess.check_output(['git','ls-files','--','docs/validation/v125'],cwd=ROOT).decode().splitlines())
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
    report=dict(passed=True,evidenceFiles=len(m['files']),sourceBindings=len(m['sourceBindings']),
        indexedEvidenceFiles=len(names),allDeliveryFilesIndexed=True,nativeTerminal=1,matlabTerminal=2,matlabExitCodes=[1,0],
        storageSHA256Rechecked=storage['checkedFiles'],unmergedResourcesPreserved=True,
        generalSSTInterfacePassed=False,mainFiveStepStartTested=False,physicsQualified=False,newQualifiedArea=0,
        verifierSHA256=sha(Path(__file__).resolve()))
    target=DEST/'delivery_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
