"""独立对照新存储包与v183既存小例，检查容量、数值等价和来源拒绝，不重跑控制。"""
from pathlib import Path
import argparse
import copy
import csv
import json
import shutil
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from runEilmerProcess import sha
from readEilmerNative import read_field, export_loads, audit_run_scope
from auditEilmerCompletion import audit_completion
from eilmerDenseRuntime import audit_dense_runtime, audit_dense_execution
from eilmerExplicitProfile import validate_explicit
from inspectEilmerRun import inspect


def verify(out):
    """模块一：绑定既存控制和本次实际二进制，先查真实完成再作逐字段比较。"""
    start = time.monotonic()
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out / 'executed_verifier.py')
    result = dict(passed=False, newNativeRuns=0, physicalQualification=False, newQualifiedAreaM2=0)
    try:
        base = ROOT / 'outputs/v185_validation'; case = base / 'small01'
        old = ROOT / 'outputs/v183_validation/small01'
        rp = base / 'small01_run/report.json'; r = json.loads(rp.read_text())
        m = json.loads((case / 'case_manifest.json').read_text()); c = m['condition']
        expected = json.loads((ROOT / 'outputs/v183_validation/small01_recovery/report.json').read_text())
        assert expected['passed'] and r['passed'] and r['commands'][0]['passed']
        assert not r['commands'][0]['liveGroupAfterCleanup'] and r['commands'][0]['inputsStable']
        for p, h in expected['outputSHA256'].items(): assert sha(p) == h
        for p, h in r['outputSHA256'].items(): assert sha(p) == h
        locked = {str(p): sha(p) for folder in [case, old] for p in folder.rglob('*') if p.is_file()}
        audit_dense_execution(c, r)
        log = base / 'small01_run/process00/stdout.txt'
        completion = audit_completion(case, log)
        assert completion['savedSteps'] == [0, 3, 6] and completion['wallFaces'] == 7
        assert inspect(case, rp, out / 'inspection.json')['passed']
        for name in ['lmrsim/config', 'lmrsim/control', 'native-input-state.json', 'lmrsim/grid/grid-0000', 'job.lua']:
            assert (case/name).read_bytes() == (old/name).read_bytes(), name
        compared = 0; errors = {}; byte_matches = {}

        def compare(name, a, b, initial_dt=False):
            nonlocal compared
            assert a.shape == b.shape, name
            if initial_dt:
                assert np.array_equal(np.isnan(a), np.isnan(b)) and not np.isinf(a).any() and not np.isinf(b).any()
                finite = np.isfinite(a); a, b = a[finite], b[finite]
            else:
                assert np.isfinite(a).all() and np.isfinite(b).all(), name
            scale = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
            error = float(np.max(np.abs(a-b)/scale)) if a.size else 0.0
            assert error <= 2e-12, (name, error)
            compared += a.size; errors[name] = error

        # 模块二：物理场、局部时间步、全部壁载荷和残差；wall-clock只作资源数据。
        for snapshot in (0, 1, 2):
            a, b = read_field(case, snapshot, 4), read_field(old, snapshot, 4)
            assert list(a) == list(b)
            path = 'lmrsim/snapshots/%04d/fluid-0000' % snapshot
            byte_matches[path] = (case/path).read_bytes() == (old/path).read_bytes()
            for field in a:
                compare('%d/%s' % (snapshot, field), a[field], b[field], snapshot == 0 and field == 'dt_local')
        newloads = sorted((case/'lmrsim/loads').glob('*/*-wall.dat'))
        oldloads = sorted((old/'lmrsim/loads').glob('*/*-wall.dat'))
        assert [p.relative_to(case) for p in newloads] == [p.relative_to(old) for p in oldloads]
        for path, control in zip(newloads, oldloads):
            a = np.atleast_1d(np.genfromtxt(path, names=True, deletechars=''))
            b = np.atleast_1d(np.genfromtxt(control, names=True, deletechars=''))
            assert a.dtype.names == b.dtype.names
            name = path.relative_to(case).as_posix()
            byte_matches[name] = path.read_bytes() == control.read_bytes()
            for field in a.dtype.names: compare(name+'/'+field, a[field], b[field])
        hist = np.loadtxt(case/'lmrsim/transient-residuals.dat')
        oldhist = np.loadtxt(old/'lmrsim/transient-residuals.dat')
        assert hist.shape == oldhist.shape == (6, 10)
        compare('residuals-excluding-wall-clock', np.delete(hist, 2, 1), np.delete(oldhist, 2, 1))
        meta = np.loadtxt(case/'lmrsim/loads/loads-metadata')
        oldmeta = np.loadtxt(old/'lmrsim/loads/loads-metadata')
        compare('loads-step-and-pseudo-time', meta, oldmeta)
        wall = export_loads(case, newloads[-1], out/'wall.csv', rp)
        rows = list(csv.DictReader((out/'wall.csv').open()))
        native = np.atleast_1d(np.genfromtxt(newloads[-1], names=True, deletechars=''))
        csv_error = 0.0
        for row in rows:
            j = int(row['nativeBoundaryRow'])
            q = -native['outsign'][j]*(native['q_cond'][j]+native['q_diff'][j])
            csv_error = max(csv_error, abs(float(row['heatFluxIntoWall_W_m2'])-q)/max(1,abs(q)))
            assert 'rah-explicit-dense-v1' in row['method'] and row['status'] == 'research'
            assert row['givenWallTemperature_K'] == '300.0' and row['physicsQualified'] == 'False'
        assert csv_error <= 2e-12 and wall['wallFaces'] == 7

        # 模块三：错误包、旧二进制、容量/越界文本和未登记完整域必须拒绝。
        rejected = []
        def reject(label, fn):
            try: fn()
            except ValueError: rejected.append(label)
            else: raise AssertionError('未拒绝: '+label)
        bad = copy.deepcopy(c); bad['explicit']['runtimePackageSHA256'] = '0'*64
        reject('wrong-package-SHA', lambda: audit_dense_runtime(bad))
        badpath = copy.deepcopy(c); badpath['explicit']['runtimePackage'] = str(out/'absent.json')
        reject('unregistered-package-path', lambda: audit_dense_runtime(badpath))
        badrun = copy.deepcopy(r); badrun['commands'][0]['argv'][0] = 'lmr-run'
        reject('old-installed-binary', lambda: audit_dense_execution(c,badrun))
        missing = copy.deepcopy(r); missing['commands'][0]['inputSHA256'].pop(c['explicit']['runtimePackage'])
        reject('missing-runtime-lock', lambda: audit_dense_execution(c,missing))
        badscope = copy.deepcopy(m); badscope['mesh']['cells'] = 705490
        reject('unregistered-full-domain', lambda: audit_run_scope(badscope))
        long = copy.deepcopy(c['explicit']); long['maxSteps'] = 200
        reject('unregistered-long-run', lambda: validate_explicit(long))
        text = log.read_text()
        line = next(line for line in text.splitlines() if line.startswith('RAH-STORAGE'))
        allocation = completion['storageAllocation']
        variants = {
            'missing-allocation-record': text.replace(line, ''),
            'old-stage-capacity': text.replace('payloadBytes='+str(allocation['arrayPayloadBytes']), 'payloadBytes='+str(allocation['oldArrayPayloadBytes'])),
            'wrong-ghost-count': text.replace('ghosts='+str(allocation['ghosts']), 'ghosts=0'),
            'native-bounds-error': text+'\ncore.exception.RangeError: Range violation\n',
        }
        for label, changed in variants.items():
            path = out/(label+'.txt'); path.write_text(changed)
            reject(label, lambda path=path: audit_completion(case,path))
        # v184失败与旧v183控制只读兼容，不让新分派追认旧失败。
        prior = ROOT/'outputs/v184_validation/real01'
        reject('prior-OOM-still-failed', lambda: audit_completion(prior,ROOT/'outputs/v184_validation/real01_run/process00/stdout.txt'))
        prior_manifest = json.loads((prior/'case_manifest.json').read_text())
        reject('prior-development-still-closed', lambda: audit_run_scope(prior_manifest,prior))
        assert audit_completion(old,ROOT/'outputs/v183_validation/small01_run/process00/stdout.txt')['finalStep'] == 6
        assert all(sha(p)==h for p,h in locked.items())
        result.update(passed=True, comparisonTolerance=2e-12, comparedScalars=int(compared), maxScaledError=max(errors.values()),
            fieldErrors=errors, fieldAndWallBytesIdentical=byte_matches, storageAllocation=allocation,
            residualWallClockExcluded=True, csvScaledError=csv_error,
            heatSigns=dict(positive=int(np.sum(-native['outsign']*native['q_total']>0)), negative=int(np.sum(-native['outsign']*native['q_total']<0))),
            lastOverFirstMassResidual=float(hist[-1,3]/hist[0,3]), lastOverFirstEnergyResidual=float(hist[-1,7]/hist[0,7]),
            sourceCasesUnchanged=True, rejected=rejected, actualArrayBoundsChecksEnabled=True,
            wholeDomainMemoryFitVerified=False, sameUpdateFormulas=True,
            sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__), *[ROOT/'scripts'/n for n in
                ['eilmerDenseRuntime.py','eilmerExplicitProfile.py','readEilmerNative.py','runPreparedEilmerCase.py','auditEilmerExplicitCompletion.py','inspectEilmerRun.py']]]})
    except Exception as error:
        result['error'] = repr(error)
        raise
    finally:
        result['elapsedSeconds'] = time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','comparedScalars','maxScaledError','storageAllocation','elapsedSeconds']}))


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--out',required=True); verify(p.parse_args().out)
