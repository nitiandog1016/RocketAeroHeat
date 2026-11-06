"""下载官方固定安装包到E盘；按范围传输，逐包SHA256校验后才交给安装程序。"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import subprocess
import time

ROOT = Path(r'E:\RocketAeroHeat\tools\runtime\setup-v177')
ASSETS = [
    ('wsl.3.0.3.0.x64.msi', 367575040,
     'ae977286974255bd26bf925f00f40f5692c6b2084155cbc0563c09aea091f211',
     'https://github.com/microsoft/WSL/releases/download/3.0.3/wsl.3.0.3.0.x64.msi'),
    ('ubuntu-24.04.5-wsl-amd64.wsl', 388975696,
     'bb415d824822c4b878125729af451a5d18fb13d1cf5cbed9a7393ad64ac6039e',
     'https://releases.ubuntu.com/24.04.5/ubuntu-24.04.5-wsl-amd64.wsl'),
    ('ldc2-1.41.0-linux-x86_64.tar.xz', 89108296,
     '4a439457f0fe59e69d02fd6b57549fc3c87ad0f55ad9fb9e42507b6f8e327c8f',
     'https://github.com/ldc-developers/ldc/releases/download/v1.41.0/ldc2-1.41.0-linux-x86_64.tar.xz'),
]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    target = ROOT / 'downloads'
    target.mkdir(parents=True, exist_ok=True)
    receipts = []
    for name, size, expected, url in ASSETS:
        started = time.monotonic()
        final = target / name
        if final.exists():
            assert final.stat().st_size == size and digest(final) == expected
            print('Already verified: ' + name, flush=True)
        else:
            parts = target / (name + '.parts')
            parts.mkdir(exist_ok=True)
            blocksize = 8 * 1024 * 1024

            def fetch(start):
                end = min(start + blocksize, size) - 1
                part = parts / str(start)
                if part.exists() and part.stat().st_size == end-start+1:
                    return part.stat().st_size
                argv = ['curl.exe', '--fail', '--location', '--silent', '--show-error',
                        '--range', '{}-{}'.format(start, end), '--max-time', '120',
                        '--retry', '1', '--retry-all-errors', '--output', str(part), url]
                result = subprocess.run(argv, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, timeout=250)
                if result.returncode or not part.exists() or part.stat().st_size != end-start+1:
                    raise RuntimeError('Range download failed {} {}: {}'.format(
                        name, start, result.stdout.decode('utf8', errors='replace')))
                return part.stat().st_size

            completed = 0
            last = time.monotonic()
            with ThreadPoolExecutor(max_workers=4) as pool:
                futures = [pool.submit(fetch, offset) for offset in range(0, size, blocksize)]
                for future in as_completed(futures):
                    completed += future.result()
                    if time.monotonic()-last > 15 or completed == size:
                        print('{} {:.1f}%'.format(name, completed*100/size), flush=True)
                        last = time.monotonic()
            assembled = target / (name + '.assembled')
            with assembled.open('wb') as out:
                for offset in range(0, size, blocksize):
                    with (parts / str(offset)).open('rb') as src:
                        for block in iter(lambda: src.read(1024*1024), b''):
                            out.write(block)
            assert assembled.stat().st_size == size and digest(assembled) == expected, name
            assembled.rename(final)
            # 只清理本次拼接的下载分块，已核对完整安装包；不触碰科研输出。
            for offset in range(0, size, blocksize):
                part = parts / str(offset)
                assert parts.resolve() in part.resolve().parents
                part.unlink()
            parts.rmdir()
        receipt = dict(name=name, bytes=size, sha256=expected, url=url,
                       verified=True, seconds=time.monotonic()-started)
        receipts.append(receipt)
        (ROOT / (name + '.verified.json')).write_text(json.dumps(receipt, indent=2)+'\n')
        (ROOT / 'downloads_verified.json').write_text(json.dumps(receipts, indent=2)+'\n')
        print('Verified: ' + name, flush=True)


if __name__ == '__main__':
    main()
