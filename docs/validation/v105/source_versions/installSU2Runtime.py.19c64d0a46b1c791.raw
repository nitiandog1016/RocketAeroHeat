"""安装固定官方SU2运行时到项目本地；不修改系统PATH，不下载求解结果。"""

from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
VERSION = '8.6.0'
URL = 'https://github.com/su2code/SU2/releases/download/v8.6.0/SU2_v8.6.0_win64-omp.zip'
ARCHIVE_SHA = '4bc28b656def415b4d22c8274af2473dd719aa5ebfa6bc1c95c08b1bac26508c'
ARCHIVE_BYTES = 29350764


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proxy', default=None)
    args = parser.parse_args()
    target = ROOT / 'tools/runtime/su2' / VERSION
    target.mkdir(parents=True, exist_ok=True)
    manifest = target / 'runtime.json'
    if manifest.exists():
        info = json.loads(manifest.read_text(encoding='utf-8'))
        assert info['archiveSHA256'] == ARCHIVE_SHA
        for item in info['files']:
            assert sha(target / item['path']) == item['sha256']
        print(json.dumps(dict(status='existing-runtime-verified', executable=info['executable'])))
        return
    archive = target / 'official-download.zip'
    if not archive.exists():
        handlers = []
        if args.proxy:
            handlers.append(urllib.request.ProxyHandler({'http': args.proxy, 'https': args.proxy}))
        opener = urllib.request.build_opener(*handlers)
        partial = target / 'official-download.part'
        with opener.open(urllib.request.Request(URL, headers={'User-Agent': 'RocketAeroHeat-runtime-setup'}),
                         timeout=45) as response, partial.open('wb') as output:
            count = 0
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                count += len(block)
                assert count <= ARCHIVE_BYTES, '下载超过官方记录大小'
                output.write(block)
        assert partial.stat().st_size == ARCHIVE_BYTES and sha(partial) == ARCHIVE_SHA
        partial.replace(archive)
    assert archive.stat().st_size == ARCHIVE_BYTES and sha(archive) == ARCHIVE_SHA
    # 模块二：只安装CFD主程序、动态依赖及许可证；所有ZIP路径限制在本运行时目录。
    files = []
    with zipfile.ZipFile(str(archive)) as package:
        members = package.infolist()
        selected = [item for item in members if not item.is_dir() and
                    (Path(item.filename).name.lower() == 'su2_cfd.exe' or
                     Path(item.filename).suffix.lower() == '.dll' or
                     Path(item.filename).name.upper().startswith(('COPYING', 'LICENSE')))]
        assert any(Path(item.filename).name.lower() == 'su2_cfd.exe' for item in selected)
        assert sum(item.file_size for item in selected) < 250 * 1024 * 1024
        for item in selected:
            path = (target / item.filename).resolve()
            assert target.resolve() in path.parents and not Path(item.filename).is_absolute()
            assert not (item.external_attr >> 16) & 0o170000 == 0o120000
            path.parent.mkdir(parents=True, exist_ok=True)
            raw = package.read(item)
            if path.exists():
                assert path.read_bytes() == raw
            else:
                with path.open('xb') as output:
                    output.write(raw)
            files.append(dict(path=path.relative_to(target).as_posix(), sha256=sha(path), bytes=len(raw)))
    executable = next(item['path'] for item in files if item['path'].lower().endswith('su2_cfd.exe'))
    info = dict(status='installed-hash-verified-not-physics-validated', version=VERSION,
                archiveURL=URL, archiveSHA256=ARCHIVE_SHA, archiveBytes=ARCHIVE_BYTES,
                executable=str((target / executable).resolve()), files=files,
                license='LGPL-2.1; preserve upstream license and notices',
                sourceURL='https://github.com/su2code/SU2/tree/v8.6.0',
                systemPATHModified=False, solverResultsDownloaded=False)
    manifest.write_text(json.dumps(info, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=info['status'], executable=info['executable'],
                         runtimeBytes=sum(item['bytes'] for item in files), files=files)))


if __name__ == '__main__':
    main()
