# -*- coding: utf-8 -*-
"""归档本轮小证据，校验保护文件；不拷贝大场或删除任何文件。"""
from pathlib import Path
import hashlib
import json
import shutil

root = Path(__file__).resolve().parents[2]
src = Path(__file__).resolve().parent
dest = root/'docs/validation/v189/evidence'
dest.mkdir(parents=True, exist_ok=False)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


names = ['analysis01.stdout.txt', 'analysis01_parent.json', 'analysis01_source.py',
         'analysis01_parent_source.py', 'analysis02.stdout.txt', 'analysis02_parent.json',
         'analysis02_source.py', 'analysis02_parent_source.py',
         'analysis02/report.json', 'analysis02/edge_components.json',
         'analysis02/original_face_diagnosis.csv', 'selection01.stdout.txt',
         'selection01_parent.json', 'selection01_source.py',
         'selection01/report.json', 'selection01/region_geometry.json',
         'integrity01.stdout.txt', 'integrity01_parent.json', 'integrity01/report.json',
         'integrity01/injected/report.json', 'eilmer_source_binding.json', 'deliver.py']
for name in names:
    p = src/name; target = dest/name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(p, target)
    assert sha(p) == sha(target)

# 模块一：失败目录中已写出的逐面数据与恢复结果逐字相同，只登记别名。
aliases = []
for folder in ('analysis01', 'integrity01/injected'):
    for name in ('edge_components.json', 'original_face_diagnosis.csv'):
        p = src/folder/name; canonical = src/'analysis02'/name
        assert sha(p) == sha(canonical)
        aliases.append(dict(local=str(p.relative_to(root)), publicCanonical='evidence/analysis02/'+name,
                            sha256=sha(p), bytes=p.stat().st_size))
protected = json.loads((root/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
assert all(sha(root/p) == value for p, value in protected.items())
parents = [json.loads((src/(name+'_parent.json')).read_text()) for name in ('analysis01','analysis02','selection01','integrity01')]
assert [p['exitCode'] for p in parents] == [1, 0, 0, 0]
assert all(p['elapsedSeconds'] < 60 for p in parents)
seconds = sum(p['elapsedSeconds'] for p in parents)
assert seconds < 120
geometry = json.loads((dest/'selection01/region_geometry.json').read_text())
contract = json.loads((dest.parent/'region_contract.json').read_text())
assert contract['local']['wallFaceCount'] == geometry['localWall']['faceCount']
assert abs(contract['local']['wallAreaM2']-geometry['localWall']['areaM2']) < 1e-14
assert contract['executionAuthorized'] is False
assert all(p['matchesFixedGitCommit'] for p in json.loads((dest/'eilmer_source_binding.json').read_text())['files'])

# 模块二：最终源码与执行过的准确版本分别绑定，失败不被覆盖。
source_paths = ['scripts/analyzeHybridRegionInputs.py', 'scripts/verifyHybridRegionSelection.py',
                'tests/testHybridInputIntegrity.py', 'AGENTS.md', 'docs/PROJECT_GOALS.md',
                'docs/HYBRID_ROUTE_SELECTION_V189.md', 'docs/validation/v189/PLAN.md',
                'docs/validation/v189/NEXT_DECISION.md', 'docs/validation/v189/SOURCES.json',
                'docs/validation/v189/region_contract.json']
files = list(dest.rglob('*'))
manifest = dict(schema='RAH.HybridH1Delivery.v1', passed=True,
                executedParentSeconds=seconds, budgetSeconds=120, protectedCount=len(protected),
                protectedUnchanged=True, newSolvedAreaM2=0, newPhysicalQualification=False,
                CFDRuns=0, builds=0, matlabRuns=0, UIChanges=0,
                finalSourceSHA256={p:sha(root/p) for p in source_paths},
                evidence=[dict(path=str(p.relative_to(dest.parent)).replace('\\','/'),
                               sha256=sha(p), bytes=p.stat().st_size) for p in files if p.is_file()],
                identicalLocalAliases=aliases,
                analysis02Scope='Exact executed recovery source archived; final source only adds passed=false on exceptions, separately fault-injection checked.',
                runtimeAndLargeFieldsBackedUpOnGitHub=False,
                uniqueFilesDeleted=0, cleanupsPerformed=0)
(dest.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(passed=True, parentSeconds=seconds, evidenceFiles=len(manifest['evidence']),
                     evidenceBytes=sum(p['bytes'] for p in manifest['evidence']), protectedCount=len(protected))))
