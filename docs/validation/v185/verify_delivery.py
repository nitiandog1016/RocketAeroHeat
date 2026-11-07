"""交付前核对Git原始blob与证据/源码SHA；并发修改不混入本轮提交。"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'docs/validation/v185'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git']+list(args),cwd=str(ROOT))


def verify():
    destination=BASE/'delivery_receipt.json';assert not destination.exists()
    manifest_path=BASE/'evidence_manifest.json'
    manifest=json.loads(manifest_path.read_text('utf-8'))
    # 模块一：每份证据的索引原字节与运行时保存的SHA相同。
    for record in manifest['files']:
        data=git('show',':'+record['path'])
        assert len(data)==record['bytes'] and sha(data)==record['sha256'],record['path']
    for name,digest in manifest['finalSourcesSHA256'].items():
        assert sha(git('show',':'+name))==digest,name
    assert sha(git('show',':docs/validation/v185/evidence_manifest.json'))==sha(manifest_path.read_bytes())
    # 模块二：独有大构建物不进Git，29项并发文件只读保持。
    names=git('diff','--cached','--name-only','-z').decode('utf-8').strip('\0').split('\0')
    allowed=set(manifest['finalSourcesSHA256'])|{'.gitattributes','AGENTS.md','docs/PROJECT_GOALS.md','docs/EILMER_DENSE_STORAGE_V185.md'}
    assert all(n in allowed or n.startswith('docs/validation/v185/') for n in names)
    assert len(names)>=len(manifest['files'])+len(allowed)
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    assert len(protected)==29
    for name,digest in protected.items():
        assert name.replace('\\','/') not in names
        assert sha((ROOT/name.replace('\\','/')).read_bytes())==digest,name
    git('diff','--cached','--check')
    result=dict(evidenceIndexRawSHA256Passed=True,evidenceFiles=len(manifest['files']),
        finalSourceIndexRawSHA256Passed=True,finalSources=len(manifest['finalSourcesSHA256']),
        protectedFilesUnchanged=29,stagedScopeChecked=True,stagedCountBeforeReceipt=len(names),
        parentCommit=git('rev-parse','HEAD').decode().strip(),manifestSHA256=sha(manifest_path.read_bytes()),
        physicalQualification=False,newQualifiedAreaM2=0,uniqueFilesDeleted=0)
    destination.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':verify()
