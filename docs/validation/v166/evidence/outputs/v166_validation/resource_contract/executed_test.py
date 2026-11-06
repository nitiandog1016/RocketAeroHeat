"""模块用途。核对具名900秒预算与普通720秒边界，全部测试不启动原生进程。"""
import copy
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from runEmbeddedSU2 import validate_resource_budget
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save


def main():
    base = ROOT/'outputs/v166_validation'
    manifest = read(base/'real01/case_manifest.json')
    request = dict(maximumSeconds=900, threads=2, resourceProfile='v166-single-development')
    assert validate_resource_budget(request, manifest) == (900., 2)
    assert validate_resource_budget(dict(maximumSeconds=720, threads=8), {}) == (720., 8)
    rejected = []
    for name in ('901-seconds', 'missing-profile', 'three-threads', 'wrong-runtime',
                 '201-updates', 'changed-prior', '1537-mib'):
        r, m = copy.deepcopy(request), copy.deepcopy(manifest)
        if name == '901-seconds': r['maximumSeconds'] = 901
        elif name == 'missing-profile': r.pop('resourceProfile')
        elif name == 'three-threads': r['threads'] = 3
        elif name == 'wrong-runtime': m['requiredRuntimeFlavor'] = 'rah-ideal-air-step-v1'
        elif name == '201-updates': m['idealAir']['checkpointContract']['maximumUpdates'] = 201
        elif name == 'changed-prior': m['initialization']['priorSuccessfulUpdates'] = 201
        else: m['idealAir']['checkpointContract']['maximumOutputBytes'] = 1537*1024**2
        try: validate_resource_budget(r, m)
        except AssertionError: rejected.append(name)
    assert len(rejected) == 7
    destination = base/'resource_contract'; destination.mkdir(exist_ok=False)
    (destination/'executed_test.py').write_bytes(Path(__file__).read_bytes())
    save(destination/'report.json', dict(passed=True, rejected=rejected, ordinaryLimitUnchanged=True,
        sourceSHA256={str(p):sha(p) for p in (Path(__file__).resolve(), ROOT/'scripts/runEmbeddedSU2.py')}, nativeRuns=0))
    print(json.dumps(dict(passed=True, rejected=rejected)))


if __name__ == '__main__': main()
