"""隐藏运行 MATLAB，绑定源码与显式输入，只在实际进程结束后记录终态。"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('name')
    parser.add_argument('expression')
    parser.add_argument('--input', action='append', default=[])
    parser.add_argument('--folder', default='outputs/v91_validation')
    args = parser.parse_args()
    assert re.fullmatch(r'[A-Za-z0-9_]+', args.name)
    folder = (ROOT / args.folder).resolve()
    folder.relative_to(ROOT / 'outputs')
    folder.mkdir(parents=True, exist_ok=True)
    receipt = folder / (args.name + '_process.json')
    assert not receipt.exists(), '已有进程记录，不得盲目重启。'
    # 全部计算/界面/测试源码锁定；只读显式输入，避免反复读取所有历史大矩阵。
    files = list(ROOT.glob('*.m'))
    for directory, suffixes in [('internalFlow', {'.m'}), ('regionalFlow', {'.m'}),
                                ('tests', {'.m', '.py', '.js'}), ('web', {'.js', '.html', '.css'}),
                                ('scripts', {'.py', '.ps1'})]:
        files += [p for p in (ROOT / directory).rglob('*') if p.is_file() and p.suffix in suffixes]
    files += [ROOT / p for p in args.input]
    files += list((ROOT / 'tests/data/temperature_research').glob('*'))
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(set(files)) if p.is_file()}
    for path in args.input:
        assert (ROOT / path).resolve().is_file(), path
    state = dict(status='running', sourceLockScope='all-sources-and-explicit-inputs',
                 sourceSHA256=before, explicitInputs=args.input, expression=args.expression,
                 startedUTC=datetime.now(timezone.utc).isoformat(), runnerProcessID=os.getpid())

    def save():
        receipt.write_bytes((json.dumps(state, indent=2) + '\n').encode('utf8'))

    save()
    log = folder / (args.name + '_stdout.log')
    with log.open('wb') as stream:
        command = ['D:/matlab/bin/matlab.exe', '-batch',
                   "restoredefaultpath; maxNumCompThreads(1); fprintf('RAH_ENGINE_PID=%d\\n',feature('getpid')); " + args.expression,
                   '-logfile', str(folder / (args.name + '.log'))]
        process = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW,
                                   env=dict(os.environ, MATLABPATH='', LC_ALL='C', LC_CTYPE='C', LANG='C'))
        state['processID'] = process.pid
        save()
        code = process.wait()
    stable = {path: sha(ROOT / path) == expected for path, expected in before.items()}
    state.update(status='terminal', exitCode=code, sourceStableAfter=stable,
                 endedUTC=datetime.now(timezone.utc).isoformat(),
                 engineProcessIDs=sorted(set(map(int, re.findall(r'RAH_ENGINE_PID=(\d+)', log.read_text('utf8', errors='replace'))))))
    save()
    print(json.dumps(dict(name=args.name, exitCode=code, sourceStable=all(stable.values()))), flush=True)
    raise SystemExit(code if code else (0 if all(stable.values()) else 2))


if __name__ == '__main__':
    main()
