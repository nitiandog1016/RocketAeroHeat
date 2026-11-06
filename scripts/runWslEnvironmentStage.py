"""监督独立WSL环境安装步骤；保存准确输出与退出码，不运行科研CFD。"""
from pathlib import Path
import argparse
import datetime
import json
import subprocess
import time

ROOT = Path(r'E:\RocketAeroHeat\tools\runtime\setup-v177')
WSL = r'E:\RocketAeroHeat\tools\runtime\windows-wsl\wsl.exe'
DISTRO = 'RocketAeroHeat-Ubuntu-24.04'
LINUX_PROJECT = '/mnt/e/RocketAeroHeat'
PREFIX = '/opt/rocket-aeroheat/eilmer-v5.0.0'
ARCHIVE = LINUX_PROJECT + '/tools/runtime/setup-v177/downloads/ldc2-1.41.0-linux-x86_64.tar.xz'


def decode_output(data):
    # WSL自己的诊断可能是UTF16，后接Linux程序的UTF8；原字节始终保留。
    if b'\x00' not in data[:200]:
        return data.decode('utf8', errors='replace')
    position = 0
    while True:
        line_end = data.find(b'\n\x00', position)
        if line_end < 0:
            return data.decode('utf-16-le', errors='replace')
        position = line_end + 2
        tail = data[position:]
        if tail and b'\x00' not in tail[:80]:
            return (data[:position].decode('utf-16-le', errors='replace')
                    + tail.decode('utf8', errors='replace'))


def main(stage, attempt=1):
    base = [WSL, '-d', DISTRO, '-u']
    stages = {
        'version': ([WSL, '--version'], 30),
        'import': ([WSL, '--import', DISTRO,
                    str(Path(r'E:\RocketAeroHeat\tools\runtime\wsl') / DISTRO),
                    str(ROOT/'downloads/ubuntu-24.04.5-wsl-amd64.wsl'), '--version', '2'], 600),
        'dependencies': (base+['root', '--exec', '/usr/bin/timeout', '1900', '/bin/bash',
                              LINUX_PROJECT+'/scripts/setupEilmerLinuxHost.sh'], 1950),
        'build': (base+['rah', '--exec', '/usr/bin/timeout', '1900', '/bin/bash',
                       LINUX_PROJECT+'/scripts/setupEilmerRuntime.sh', PREFIX, ARCHIVE], 1950),
        'health': (base+['rah', '--exec', '/usr/bin/timeout', '120', '/bin/bash',
                        LINUX_PROJECT+'/scripts/checkEilmerRuntime.sh',
                        'installation-health-{:02d}'.format(attempt)], 150),
    }
    argv, limit = stages[stage]
    assert attempt >= 1
    label = stage if attempt == 1 else '{}{:02d}'.format(stage, attempt)
    receipt = ROOT/(label+'-result.json')
    raw_path = ROOT/(label+'-output.raw')
    error_path = ROOT/(label+'-error.raw')
    if receipt.exists() or raw_path.exists() or error_path.exists():
        raise RuntimeError('Refusing to overwrite a previous stage receipt: '+stage)
    result = {'stage':stage, 'attempt':attempt, 'argv':argv, 'limitSeconds':limit, 'status':'running',
              'startedUTC':datetime.datetime.utcnow().isoformat()+'Z', 'newCFDTasks':0}
    receipt.write_text(json.dumps(result, indent=2)+'\n')
    started = time.monotonic()
    try:
        with raw_path.open('wb') as out, error_path.open('wb') as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err)
            result['windowsProcessId'] = process.pid
            receipt.write_text(json.dumps(result, indent=2)+'\n')
            try:
                code = process.wait(timeout=limit)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
                # 安装依赖和构建各有Linux端timeout；这里只终止自身Windows前端。
                raise
        result.update(exitCode=code, status='passed' if code == 0 else 'failed')
    except Exception as error:
        result.update(status='failed', error=str(error))
    result['elapsedSeconds'] = time.monotonic()-started
    receipt.write_text(json.dumps(result, indent=2)+'\n')
    data = raw_path.read_bytes()
    output = decode_output(data)
    (ROOT/(label+'-output.txt')).write_text(output, encoding='utf8')
    errors = error_path.read_bytes()
    error_text = decode_output(errors)
    (ROOT/(label+'-error.txt')).write_text(error_text, encoding='utf8')
    print(output[-5000:])
    if error_text:
        print(error_text[-2000:])
    print(json.dumps(result))
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['version','import','dependencies','build','health'])
    parser.add_argument('--attempt', type=int, default=1)
    options = parser.parse_args()
    raise SystemExit(main(options.stage, options.attempt))
