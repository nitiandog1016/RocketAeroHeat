# -*- coding: utf-8 -*-
"""用受控SHA故障验证输入变化时报告不能标为通过，不修改原输入。"""
from pathlib import Path
import importlib.util
import json
import sys
import time

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('hybrid', root/'scripts/analyzeHybridRegionInputs.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
started = time.monotonic()
real_sha = module.sha
seen = {}


def changed_summary_sha(path):
    """仅注入第二次摘要哈希读取，磁盘原数据和真实哈希函数保持。"""
    key = str(path)
    seen[key] = seen.get(key, 0)+1
    if Path(path).name == 'summary.json' and seen[key] > 1:
        return '0'*64
    return real_sha(path)


module.sha = changed_summary_sha
rejected = False
try:
    module.analyze(out/'injected')
except AssertionError:
    rejected = True
finally:
    module.sha = real_sha
report = json.loads((out/'injected/report.json').read_text(encoding='utf-8'))
assert rejected and report['passed'] is False and report['inputStable'] is False
assert report['newSolvedAreaM2'] == 0 and report['newPhysicalQualification'] is False
assert all(real_sha(path) == h for path, h in report['sourceSHA256'].items())
result = dict(passed=True, test='injected second source SHA mismatch, no real file mutation',
              rejected=rejected, qualificationStillFalse=True, realInputsUnchanged=True,
              elapsedSeconds=time.monotonic()-started)
(out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps(result))
