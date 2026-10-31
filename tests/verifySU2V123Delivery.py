"""核对交付证据的原字节、实际终态与Git暂存版本；不把接口成功当物理资格。"""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v123'
def digest(raw):return hashlib.sha256(raw).hexdigest()
def main():
    m=json.loads((DEST/'delivery_manifest.json').read_text())
    for q in m['files']:assert digest((ROOT/q['archived']).read_bytes())==q['sha256']
    rawhash={digest(p.read_bytes()) for p in (DEST/'source_versions').glob('*.raw')}
    assert all(q['sha256'] in rawhash for q in m['sourceBindings'])
    assert [r['exitCode'] for r in m['matlab']]==[1,0] and all(r['sourceStable'] for r in m['matlab'])
    assert m['nativeEngineRuns']==1 and sorted(r['exitCode'] for r in m['pipelines'])==[0,1,125]
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v123'],cwd=ROOT).decode().splitlines();assert names
    indexed=set(subprocess.check_output(['git','ls-files','--','docs/validation/v123'],cwd=ROOT).decode().splitlines())
    required={q['archived'].replace('\\','/') for q in m['files']}
    required.update(p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw'))
    assert required <= indexed, '交付清单仍有未暂存文件：'+str(sorted(required-indexed))
    raw=subprocess.check_output(['git','cat-file','--batch'],cwd=ROOT,input=''.join(':'+p+'\n' for p in names).encode())
    cursor=0
    for name in names:
        end=raw.index(b'\n',cursor);head=raw[cursor:end].decode().split();assert head[1]=='blob';length=int(head[2])
        assert digest(raw[end+1:end+1+length])==digest((ROOT/name).read_bytes()),name
        cursor=end+1+length+1
    assert cursor==len(raw)
    report=dict(passed=True,evidenceFiles=len(m['files']),sourceVersions=len(rawhash),stagedEvidenceFiles=len(names),
        allManifestFilesIndexed=True,allStagedEvidenceByteExact=True,matlabTerminal=2,nativeEngineTerminal=1,pipelineTerminal=3,
        syntheticSupervisionChildren=3,newQualifiedArea=0,physicsQualified=False,verifierSHA256=digest(Path(__file__).read_bytes()))
    file=DEST/'staged_evidence_check_final.json';assert not file.exists();file.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
