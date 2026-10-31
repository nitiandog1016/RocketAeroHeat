"""独立核对本轮终态、准确归档、后验拒绝、存储身份及Git暂存；不授予物理资格。"""
from pathlib import Path
import hashlib,json,subprocess,os,re
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v124'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    m=read(DEST/'delivery_manifest.json')
    for item in m['files']:assert sha(ROOT/item['archived'])==item['sha256']
    for item in m['sourceBindings']:assert sha(ROOT/item['archived'])==item['sha256']
    assert [x['exitCode'] for x in m['matlab']]==[0,1,1,1,0] and all(x['sourceStable'] for x in m['matlab'])
    assert m['nativeEngineRuns']==2 and m['nativeEngineExitCodes']==[0,0] and m['fullDomainPipelineExitCode']==1
    assert not m['fullDomainResultAccepted'] and not m['physicsQualified'] and m['newQualifiedArea']==0
    a=read(DEST/'evidence/split_phase_diagnosis.json');assert not a['gradientGatePassed'] and a['maximumGradientRelativeError']>a['originalGradientThreshold']
    assert a['yPlusRelativeError']<2e-6 and a['last25HeatPowerRelativeChange']>.02
    assert read(DEST/'evidence/controls04/report.json')['existingHEOSReadOnlyPassed']
    storage=read(DEST/'storage_receipt.json')
    for p,h in storage['originalSHA256'].items():assert sha(Path(p))==h
    for row in storage['plan']:assert sha(Path(row['source']))==row['sha256'] and os.path.samefile(row['source'],row['destination'])
    # 当前目标中所有本地Markdown目标均存在，旧历史记录不用覆盖当前第4B起点。
    for file in [ROOT/'docs/PROJECT_GOALS.md',ROOT/'docs/SU2_SST_INTERFACE_V124.md',ROOT/'docs/reviews/TECHNICAL_ROUTE_REVIEW_2026-11-01.md']:
        for link in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
            if '://' not in link:assert (file.parent/link.split('#')[0]).exists(),link
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v124'],cwd=ROOT).decode().splitlines();assert names
    indexed=set(subprocess.check_output(['git','ls-files','--','docs/validation/v124'],cwd=ROOT).decode().splitlines())
    required={p.relative_to(ROOT).as_posix() for p in DEST.rglob('*') if p.is_file()}
    assert required<=indexed,sorted(required-indexed)
    raw=subprocess.check_output(['git','cat-file','--batch'],cwd=ROOT,input=''.join(':'+p+'\n' for p in names).encode());cursor=0
    for name in names:
        end=raw.index(b'\n',cursor);head=raw[cursor:end].decode().split();assert head[1]=='blob';length=int(head[2])
        assert hashlib.sha256(raw[end+1:end+1+length]).hexdigest()==sha(ROOT/name),name
        cursor=end+1+length+1
    assert cursor==len(raw)
    report=dict(passed=True,evidenceFiles=len(m['files']),sourceVersions=len(list((DEST/'source_versions').glob('*.raw'))),
        exactIndexedEvidenceFiles=len(names),allDeliveryFilesIndexed=True,storageFilesRechecked=storage['checkedFiles'],
        matlabTerminal=5,matlabExitCodes=[0,1,1,1,0],nativeTerminal=2,nativeExitCodes=[0,0],fullPipelineExitCode=1,
        fullDomainResultAccepted=False,physicsQualified=False,newQualifiedArea=0,verifierSHA256=sha(Path(__file__).resolve()))
    target=DEST/'delivery_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
