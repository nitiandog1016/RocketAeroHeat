"""独立核对小型交付与Git暂存原字节，防止证据因换行规则再次改变。"""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/validation/v122'
def digest(raw):return hashlib.sha256(raw).hexdigest()
def main():
    manifest=json.loads((DEST/'delivery_manifest.json').read_text());count=0
    for row in manifest['files']:
        assert digest((ROOT/row['archived']).read_bytes())==row['sha256'];count+=1
    names=subprocess.check_output(['git','diff','--cached','--name-only','--','docs/validation/v122'],cwd=ROOT).decode().splitlines()
    assert names
    raw=subprocess.check_output(['git','cat-file','--batch'],cwd=ROOT,input=''.join(':'+n+'\n' for n in names).encode())
    offset=0;staged=[]
    for name in names:
        end=raw.index(b'\n',offset);head=raw[offset:end].decode().split();assert head[1]=='blob'
        length=int(head[2]);start=end+1;body=raw[start:start+length];offset=start+length+1
        assert digest(body)==digest((ROOT/name).read_bytes()),name
        staged.append((name,digest(body)))
    assert offset==len(raw)
    engines=manifest['engines'];assert len(engines)==7 and engines['actual04']['exitCode']==0
    assert sum(x['exitCode']!=0 for x in engines.values())==3
    assert manifest['matlab']['review01_process']['exitCode']==1 and manifest['matlab']['review02_process']['exitCode']==0
    report=dict(passed=True,manifestFiles=count,stagedEvidenceFiles=len(names),allStagedEvidenceByteExact=True,
        stagedNameAndDigestSHA256=digest(json.dumps(staged).encode()),enginesTerminal=7,matlabTerminal=2,
        physicsQualified=False,newQualifiedArea=0,verifierSHA256=digest(Path(__file__).read_bytes()))
    target=DEST/'staged_evidence_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
