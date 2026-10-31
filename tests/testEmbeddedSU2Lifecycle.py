"""v107监督器实际失败/超时及历史行解析；不运行收敛物理基准。"""
from pathlib import Path
import importlib.util
import json
import shutil
import subprocess
import sys
import hashlib
import argparse

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='本程序原生箭体短迭代目录，至少有两行历史')
    parser.add_argument('--out', required=True)
    parser.add_argument('--audit-existing', action='store_true', help='只审核已有终态，不再启动引擎')
    args = parser.parse_args()
    source = Path(args.source).resolve()
    out = Path(args.out).resolve()
    out.relative_to(ROOT / 'outputs')
    if not args.audit_existing:
        out.mkdir(exist_ok=False)
    script = ROOT / 'scripts/runEmbeddedSU2.py'
    spec = importlib.util.spec_from_file_location('supervisor', str(script))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    progress = module.history_progress(source / 'history.csv')
    # 独立CSV读取全历史作对照，覆盖带空格/引号的SU2实际表头。
    import csv
    with (source / 'history.csv').open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    reference = {k.strip(): float(v) for k, v in zip(rows[0], rows[-1])}
    assert progress == {key: reference[key] for key in ('Inner_Iter','rms[Rho]','rms[RhoE]','HF')}
    partial = out / 'partial.csv'
    partial.write_text('"Inner_Iter",   "rms[Rho]"\n1,-3\n2,', encoding='utf8')
    assert module.history_progress(partial) == {'Inner_Iter': 1.0, 'rms[Rho]': -3.0}
    results = []
    for label in ('invalid_config', 'time_limit'):
        folder = out / label
        code = None
        if not args.audit_existing:
            folder.mkdir()
            shutil.copyfile(str(source / 'mesh.su2'), str(folder / 'mesh.su2'))
            import re
            cfg = (source / 'case.cfg').read_text(encoding='utf8')
            cfg = re.sub(r'^ITER=.*$', 'ITER= 30000', cfg, flags=re.M)
            cfg = re.sub(r'^OUTPUT_WRT_FREQ=.*$', 'OUTPUT_WRT_FREQ= 1000', cfg, flags=re.M)
            if label == 'invalid_config':
                cfg += 'INTENTIONALLY_INVALID_V107_OPTION= YES\n'
            (folder / 'case.cfg').write_text(cfg, encoding='utf8')
            request = dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
                           threads=2, maximumSeconds=1 if label == 'time_limit' else 30,
                           jobToken='actual-owned-v107-' + label)
            (folder / 'engine_cancel.json').write_text(json.dumps(dict(jobToken='wrong-token')), encoding='utf8')
            path = folder / 'engine_request.json'
            path.write_text(json.dumps(request), encoding='utf8')
            completed = subprocess.run([sys.executable, '-X', 'utf8', str(script), '--request', str(path)],
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       creationflags=subprocess.CREATE_NO_WINDOW)
            (folder / 'supervisor.log').write_bytes(completed.stdout)
            code = completed.returncode
        state = json.loads((folder / 'engine_process.json').read_text(encoding='utf8'))
        assert state['status'] == 'terminal' and not state['cancelled'] and all(state['inputsStable'].values())
        assert state['timedOut'] == (label == 'time_limit')
        if code is not None:
            assert code == 124 if label == 'time_limit' else code != 0
        results.append(dict(case=label, supervisorExitCodeObserved=code, terminal=state))
    report = dict(progress=progress, partialRowIgnored=True, wrongCancellationTokenIgnored=True,
                  cases=results, physicsQualified=False, reusedExistingTerminal=args.audit_existing,
                  sourceSHA256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (script, Path(__file__).resolve(), source / 'case.cfg', source / 'mesh.su2')})
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf8')
    print('Actual failure, timeout, wrong-token refusal and complete-row progress passed.')


if __name__ == '__main__':
    main()
