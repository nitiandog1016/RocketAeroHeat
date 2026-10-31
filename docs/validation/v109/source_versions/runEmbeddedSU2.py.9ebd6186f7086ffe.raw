"""隐藏执行内嵌SU2，保存真实进程终态和输入散列；退出0不等于物理解合格。"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import subprocess
import time
import csv
import math

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
    for name in request.get('sourceFiles', []):
        file = (ROOT / name).resolve()
        file.relative_to(ROOT)
        assert file.is_file(), '缺少绑定的输入/源码'
        locked[str(file)] = sha(file)
    # 重启动初场是数值输入，必须和配置/网格一起绑定，输出不得覆写它。
    import re
    config_text = cfg.read_text(encoding='utf-8-sig')
    if re.search(r'^\s*RESTART_SOL\s*=\s*YES', config_text, re.M):
        match = re.search(r'^\s*SOLUTION_FILENAME\s*=\s*([^\s%]+)', config_text, re.M)
        assert match, '重启动必须显式指定SOLUTION_FILENAME'
        # SU2 GetFilename为稳态ASCII输入添加.csv；v107按字面文件名锁定并未实测重启动。
        assert re.search(r'^\s*READ_BINARY_RESTART\s*=\s*NO', config_text, re.M), '当前重启动只支持稳态ASCII'
        assert not re.search(r'^\s*TIME_DOMAIN\s*=\s*YES', config_text, re.M), '非定常重启动尚未实现'
        stem = match.group(1)
        assert Path(stem).suffix == '', 'SOLUTION_FILENAME须用无扩展名的文件基名'
        initial = (folder / (stem + '.csv')).resolve()
        assert initial.parent == folder and initial.is_file()
        target = re.search(r'^\s*RESTART_FILENAME\s*=\s*([^\s%]+)', config_text, re.M)
        assert target and target.group(1) != stem, '检查点不得覆写初始场'
        locked[str(initial)] = sha(initial)
    # 独占创建防止两个监督器同时启动；回执采用原子替换，UI不会读到半份JSON。
    with receipt.open('x', encoding='utf-8') as stream:
        stream.write('{}')
    state = dict(status='starting', startedUTC=datetime.now(timezone.utc).isoformat(),
                 request=request, inputSHA256=locked, executableVersion=runtime['version'],
                 archiveSHA256=runtime['archiveSHA256'], supervisorProcessID=os.getpid(),
                 physicalSolutionQualified=False)

    def save():
        temporary = receipt.with_suffix('.pending')
        temporary.write_text(json.dumps(state, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        os.replace(str(temporary), str(receipt))

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
        cancelled = False
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - start
                cancel = folder / 'engine_cancel.json'
                if cancel.exists() and request.get('jobToken'):
                    try:
                        cancelled = json.loads(cancel.read_text(encoding='utf-8-sig')).get('jobToken') == request['jobToken']
                    except (ValueError, OSError):
                        pass
                timed_out = elapsed >= timeout
                if timed_out or cancelled:
                    process.kill()
                    break
                state.update(elapsedSeconds=elapsed, progress=history_progress(folder / 'history.csv'))
                save()
                time.sleep(0.25)
            code = process.wait()
        finally:
            # 包括回执/磁盘异常；只回收本监督器拥有的子进程，绝不按进程名批量杀。
            if process.poll() is None:
                process.kill()
                process.wait()
    state.update(status='terminal', exitCode=code, timedOut=timed_out,
                 cancelled=cancelled, progress=history_progress(folder / 'history.csv'),
                 elapsedSeconds=time.monotonic() - start, endedUTC=datetime.now(timezone.utc).isoformat(),
                 inputsStable={path: sha(Path(path)) == digest for path, digest in locked.items()})
    state['outputSHA256'] = {file.name: sha(file) for file in sorted(folder.iterdir())
                            if file.is_file() and file != receipt and str(file) not in locked}
    save()
    print(json.dumps(dict(status='terminal', exitCode=code, timedOut=timed_out,
                         elapsedSeconds=state['elapsedSeconds'], receipt=str(receipt))))
    return 125 if cancelled else (124 if timed_out else (code if code else (0 if all(state['inputsStable'].values()) else 2)))


def history_progress(path):
    """只读最后一个完整历史行；迭代号/残差不是完成百分比或物理资格。"""
    if not path.exists():
        return {}
    try:
        with path.open('rb') as stream:
            header = stream.readline().decode('utf-8').strip()
            stream.seek(0, 2)
            size = stream.tell()
            stream.seek(max(0, size - 65536))
            lines = stream.read().decode('utf-8').splitlines(keepends=True)
        complete = [line for line in lines if line.endswith('\n')]
        keys = [key.strip() for key in next(csv.reader([header], skipinitialspace=True))]
        for line in reversed(complete):
            values = next(csv.reader([line]))
            if len(values) != len(keys):
                continue
            try:
                row = {key: float(value) for key, value in zip(keys, values)}
                return {key: value for key, value in row.items()
                        if key in ('Inner_Iter', 'rms[Rho]', 'rms[RhoE]', 'HF') and math.isfinite(value)}
            except ValueError:
                continue
    except (OSError, ValueError, UnicodeError):
        pass
    return {}


if __name__ == '__main__':
    raise SystemExit(main())
