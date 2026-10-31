"""隐藏执行内嵌SU2，保存真实进程终态和输入散列；退出0不等于物理解合格。"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text(encoding='utf-8-sig'))
    folder = Path(request['workDirectory']).resolve()
    folder.relative_to(ROOT / 'outputs')
    cfg = folder / request['configFile']
    mesh = folder / request['meshFile']
    assert cfg.parent == folder and mesh.parent == folder
    runtime = json.loads((ROOT / 'tools/runtime/su2/8.6.0/runtime.json').read_text(encoding='utf-8'))
    exe = Path(runtime['executable']).resolve()
    runtime_dir = ROOT / 'tools/runtime/su2/8.6.0'
    exe.relative_to(runtime_dir)
    for item in runtime['files']:
        assert sha(runtime_dir / item['path']) == item['sha256'], '运行时身份改变'
    timeout = float(request['maximumSeconds'])
    threads = int(request['threads'])
    assert 1 <= threads <= 8 and 1 <= timeout <= 720
    receipt = folder / 'engine_process.json'
    assert not receipt.exists(), '已有真实进程记录，禁止覆盖或盲目重启'
    locked = {str(file): sha(file) for file in (request_path, cfg, mesh, exe, Path(__file__).resolve())}
    state = dict(status='starting', startedUTC=datetime.now(timezone.utc).isoformat(),
                 request=request, inputSHA256=locked, executableVersion=runtime['version'],
                 archiveSHA256=runtime['archiveSHA256'], physicalSolutionQualified=False)

    def save():
        receipt.write_text(json.dumps(state, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    save()
    command = [str(exe), str(cfg), '--threads', str(threads)]
    if request.get('dryRun', False):
        command.append('--dryrun')
    start = time.monotonic()
    with (folder / 'engine_stdout.log').open('wb') as stream:
        process = subprocess.Popen(command, cwd=str(folder), stdout=stream, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0,
                                   env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
        state.update(status='running', processID=process.pid, command=command)
        save()
        timed_out = False
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            code = process.wait()
    state.update(status='terminal', exitCode=code, timedOut=timed_out,
                 elapsedSeconds=time.monotonic() - start, endedUTC=datetime.now(timezone.utc).isoformat(),
                 inputsStable={path: sha(Path(path)) == digest for path, digest in locked.items()})
    state['outputSHA256'] = {file.name: sha(file) for file in sorted(folder.iterdir())
                            if file.is_file() and file != receipt and str(file) not in locked}
    save()
    print(json.dumps(dict(status='terminal', exitCode=code, timedOut=timed_out,
                         elapsedSeconds=state['elapsedSeconds'], receipt=str(receipt))))
    return 124 if timed_out else (code if code else (0 if all(state['inputsStable'].values()) else 2))


if __name__ == '__main__':
    raise SystemExit(main())
