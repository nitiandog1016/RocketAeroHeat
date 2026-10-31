"""逐项公布独立原MAT恢复校验进度，不删除任何数据。"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[4]
DEST = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
from auditMat5HistoryArchive import audit_manifest


def main():
    plan = json.loads((DEST / 'plan.json').read_text(encoding='utf-8'))
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(audit_manifest, ROOT / g['manifest']['path'], ROOT / plan['store'])
                   for g in plan['groups']]
        for future in as_completed(futures):
            results.append(future.result())
            (DEST / 'independent_audit_progress.json').write_bytes(
                (json.dumps({'verified': len(results), 'total': 20, 'results': results},
                            ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
            print('VERIFIED %d/20 %s' % (len(results), results[-1]['original']), flush=True)
    report = {'status': 'passed-independent-exact-original-MAT-SHA', 'verifiedOriginals': len(results),
              'verifiedOriginalBytes': sum(x['bytes'] for x in results), 'results': results,
              'freeBytes': shutil.disk_usage(str(ROOT)).free, 'archiveReplayFunctionUsed': False,
              'runnerSHA256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'auditorSHA256': hashlib.sha256((ROOT / 'scripts/auditMat5HistoryArchive.py').read_bytes()).hexdigest()}
    (DEST / 'predelete_archive_audit.json').write_bytes(
        (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps({k: v for k, v in report.items() if k != 'results'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
