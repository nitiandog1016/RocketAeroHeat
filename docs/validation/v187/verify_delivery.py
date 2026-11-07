"""交付前核对原始Git内容、源码身份和精确暂存范围，保护并发工作。"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'docs/validation/v187'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git'] + list(args), cwd=str(ROOT))


def verify():
    destination = BASE / 'delivery_receipt.json'
    assert not destination.exists(), '已有回执不得覆盖'
    manifest_path = BASE / 'evidence_manifest.json'
    manifest = json.loads(manifest_path.read_text('utf-8'))
    # 模块一：证据和最终源码须以运行记录的原始字节进入Git。
    assert len(manifest['files']) == manifest['filesCount'] == 231
    assert len(manifest['finalSourcesSHA256']) == 13
    assert sum(r['bytes'] for r in manifest['files']) == manifest['totalSmallBytes']
    for record in manifest['files']:
        data = git('show', ':' + record['path'])
        assert len(data) == record['bytes'] and sha(data) == record['sha256'], record['path']
        assert sha((ROOT / record['path']).read_bytes()) == record['sha256'], record['path']
    for name, digest in manifest['finalSourcesSHA256'].items():
        assert sha(git('show', ':' + name)) == digest, name
        assert sha((ROOT / name).read_bytes()) == digest, name
    assert sha(git('show', ':docs/validation/v187/evidence_manifest.json')) == sha(manifest_path.read_bytes())
    # 模块二：精确列举交付文件，拒绝混入大场、运行包或其他并发修改。
    names = git('diff', '--cached', '--name-only', '-z').decode('utf-8').strip('\0').split('\0')
    allowed = set(manifest['finalSourcesSHA256']) | {r['path'] for r in manifest['files']}
    allowed |= {'.gitattributes', 'AGENTS.md', 'docs/PROJECT_GOALS.md', 'docs/EILMER_WORKSPACE_STORAGE_V187.md'}
    allowed |= {'docs/validation/v187/' + name for name in (
        'PLAN.md', 'FULL_DOMAIN_PLAN.md', 'NEXT_DECISION.md', 'README.md', 'collect_evidence.py',
        'evidence_manifest.json', 'budget_receipt.json', 'verify_delivery.py')}
    assert set(names) == allowed, dict(extra=sorted(set(names) - allowed), missing=sorted(allowed - set(names)))
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    assert len(protected) == 29
    for name, digest in protected.items():
        relative = name.replace('\\', '/')
        assert relative not in names
        assert sha((ROOT / relative).read_bytes()) == digest, name
    git('diff', '--cached', '--check')
    # 模块三：记录本次交付检查，不将来源/存储检查当作物理验证。
    result = dict(evidenceIndexRawSHA256Passed=True, evidenceFiles=len(manifest['files']),
        finalSourceIndexRawSHA256Passed=True, finalSources=len(manifest['finalSourcesSHA256']),
        protectedFilesUnchanged=29, stagedScopeChecked=True, stagedCountBeforeReceipt=len(names),
        parentCommit=git('rev-parse', 'HEAD').decode().strip(), manifestSHA256=sha(manifest_path.read_bytes()),
        physicalQualification=False, newQualifiedAreaM2=0, uniqueFilesDeleted=0)
    destination.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    verify()
