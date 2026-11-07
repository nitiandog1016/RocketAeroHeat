"""交付前核对Git原始blob、证据清单与并发保护；不触发计算或清理。"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'docs/validation/v184'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git'] + list(args), cwd=str(ROOT))


def verify():
    """模块一：逐项检查暂存原字节，避免换行转换破坏已有SHA证据。"""
    destination = BASE / 'delivery_receipt.json'
    assert not destination.exists(), '已有交付回执不覆盖'
    manifest_path = BASE / 'evidence_manifest.json'
    manifest = json.loads(manifest_path.read_text('utf-8'))
    for record in manifest['files']:
        data = git('show', ':' + record['path'])
        assert len(data) == record['bytes'] and sha(data) == record['sha256'], record['path']
    for name, digest in manifest['finalSourcesSHA256'].items():
        assert sha(git('show', ':' + name)) == digest, name
    assert sha(git('show', ':docs/validation/v184/evidence_manifest.json')) == sha(manifest_path.read_bytes())

    # 模块二：限定本轮暂存范围，保护原有未提交源码与研究文件。
    names = git('diff', '--cached', '--name-only', '-z').decode('utf-8').strip('\0').split('\0')
    allowed = set(manifest['finalSourcesSHA256']) | {
        '.gitattributes', 'AGENTS.md', 'docs/PROJECT_GOALS.md', 'docs/EILMER_EXPLICIT_FULL_DOMAIN_V184.md',
    }
    assert all(name in allowed or name.startswith('docs/validation/v184/') for name in names)
    assert len(names) >= len(manifest['files']) + len(allowed)
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    assert len(protected) == 29
    for name, digest in protected.items():
        assert name.replace('\\', '/') not in names
        assert sha((ROOT / name.replace('\\', '/')).read_bytes()) == digest, name
    git('diff', '--cached', '--check')
    result = dict(
        evidenceIndexRawSHA256Passed=True, evidenceFiles=len(manifest['files']),
        finalSourceIndexRawSHA256Passed=True, finalSources=len(manifest['finalSourcesSHA256']),
        protectedFilesUnchanged=len(protected), stagedScopeChecked=True,
        stagedCountBeforeReceipt=len(names), parentCommit=git('rev-parse', 'HEAD').decode().strip(),
        manifestSHA256=sha(manifest_path.read_bytes()), physicalQualification=False,
        newQualifiedAreaM2=0, uniqueFilesDeleted=0)
    destination.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    verify()
