"""收集v182小证据，核对保护文件与准确执行源；大构建物仅列清单，不删除。"""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import time

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'outputs/v182_validation'
DEST = ROOT / 'docs/validation/v182/evidence'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(name):
    return json.loads((BASE / name).read_text('utf-8-sig'))


def collect():
    start = time.monotonic()
    assert not DEST.exists()
    # 模块一：保护并发文件，确认诊断执行源与后加说明的区别。
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks = {p: sha(ROOT / Path(p.replace('\\', '/'))) == h for p, h in protected.items()}
    assert len(checks) == 29 and all(checks.values())
    final_d = ROOT / 'scripts/native/eilmerOperatorDiagnostic.d'
    executed_d = BASE / 'build01/executed_diagnostic.d'
    final_lines = final_d.read_text('utf8').splitlines()
    comment_lines = [line for line in final_lines if line.startswith((' * v182结果限制', ' * 该列不能', ' * 原准确执行源'))]
    assert len(comment_lines) == 3
    assert [line for line in final_lines if line not in comment_lines] == executed_d.read_text('utf8').splitlines()
    for path in [ROOT / 'tests/verifyEilmerNativeOperators.py', Path(__file__)]:
        ast.parse(path.read_text('utf8'))

    # 模块二：只读核对原始失败、终态清点与预算，不重新计算原生结果。
    build = read('build01/recovery04_report.json')
    runs = read('diagnostics_report.json')
    analysis = read('analysis01/report.json')
    scope = read('scope02.json')
    assert build['passed'] and build['totalBuildSeconds'] <= 180
    assert runs['passed'] and runs['elapsedSeconds'] <= 60
    assert analysis['passed'] is False and scope['passed'] is False
    assert scope['productionCauseOfV181Proven'] is False
    assert max(analysis['records']['real-state0']['J_mixed_vs_central']) > 1e-4
    assert max(analysis['records']['real-state2']['J_mixed_vs_central']) < 1e-4
    for item in runs['results']:
        proc = item['process']
        assert proc['inputsStable'] and not proc['liveGroupAfterCleanup']

    # 模块三：复制必要的小文件，大文件只登记位置/SHA/大小。
    entries, local = [], []
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix == '.pyc':
            continue
        source = path.relative_to(ROOT).as_posix()
        record = dict(path=source, bytes=path.stat().st_size, sha256=sha(path))
        if record['bytes'] > 1024 ** 2:
            local.append(record)
            continue
        target = DEST / path.relative_to(BASE)
        if target.suffix == '.log':
            target = target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert sha(target) == record['sha256']
        entries.append(dict(path=target.relative_to(ROOT).as_posix(), source=source,
                            bytes=record['bytes'], sha256=record['sha256']))
    names = ['scripts/native/eilmerOperatorDiagnostic.d', 'tests/verifyEilmerNativeOperators.py']
    manifest = dict(files=entries, filesCount=len(entries), totalSmallBytes=sum(x['bytes'] for x in entries),
                    largeLocalOnly=local, largeLocalOnlyBytes=sum(x['bytes'] for x in local),
                    finalSourcesSHA256={p: sha(ROOT / p) for p in names},
                    exactExecutedDiagnosticSHA256=sha(executed_d), finalDOnlyAddsThreeWarningComments=True,
                    protectedFiles=checks, uniqueFilesDeleted=0, newQualifiedAreaM2=0,
                    note='Initial derivative gate failed; scope02 supersedes the overstrong scope01 interpretation. Copied old loads/history are not new CFD steps. Large artifacts are not backed up on GitHub.')
    (DEST.parent / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf8')
    budget = dict(buildProgramSeconds=build['totalBuildSeconds'], buildProgramBudgetSeconds=180,
                  diagnosticProgramSeconds=runs['elapsedSeconds'], diagnosticProgramBudgetSeconds=60,
                  analysisProgramSeconds=sum(read(p)['elapsedSeconds'] for p in
                                             ['analysis01/report.json', 'scope01/report.json', 'scope02.json']),
                  analysisProgramBudgetSeconds=60,
                  explicitSourceReviewSeconds=read('explicit_entry01/report.json')['elapsedSeconds'],
                  explicitSourceReviewParentLimitSeconds=30,
                  collectionSeconds=time.monotonic()-start,
                  productionKernelSourceCompiles=0, newCFDSteps=0,
                  originalPlanSingleParentSeconds=210, continuousSingleParentApplied=False,
                  budgetDeviation='No single parent covered manual build recoveries; per-stage parents and inherited remaining program budget applied. Editing/reading/manual coordination wall time is separate and not reported as program time.',
                  outputDirectoryLogicalBytes=manifest['totalSmallBytes']+manifest['largeLocalOnlyBytes'],
                  outputSizeScope='Only outputs/v182_validation counted; logical bytes, not new allocated disk space. No scan of all outputs.',
                  largeArtifactsBackedUpOnGitHub=False, uniqueFilesDeleted=0,
                  derivativeGatePassed=False, physicalQualification=False, newQualifiedAreaM2=0)
    (DEST.parent / 'budget_receipt.json').write_text(json.dumps(budget, indent=2) + '\n', encoding='utf8')
    print(json.dumps(dict(files=len(entries), smallBytes=manifest['totalSmallBytes'], largeCount=len(local),
                          localBytes=manifest['largeLocalOnlyBytes'], protected=len(checks), derivativeGatePassed=False)))


if __name__ == '__main__':
    collect()
